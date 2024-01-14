/* 
 * Ишутин Андрей (Ishutin Andrei)
 * telegram: @looks_amazing
 * Github: Andy-823 
 */

#include <iostream>
#include <iomanip>

#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>
#include <eigen3/eigen/core>

/*
 *  строк в матрице не больше чем столбцов
 *  статья, ссылки не нащел: On Speeding-up Parallel Jacobi Iterations for SVD
 *  параллельности пока нет
 *  про матрицы https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
 *  возвращаемый тип такой, поскольку матрицы имеют разный размер
 *  нужны такие матрицы:
 *      Eigen::Matrix<T, N, N> U
 *      Eigen::Matrix<T, N, M> S
 *      Eigen::Matrix<T, M, M> V
 *  видно различие типов, поэтому взять std::tuple
 *
 *  TODO: изменение матриц по ссылке
 */
template<typename T, int N, int M>
std::tuple<Eigen::Matrix<T, N, N>,    // U
           Eigen::Matrix<T, N, M>,    // S
           Eigen::Matrix<T, M, M>>    // V
JTS_SVD_base(Eigen::Matrix<T, N, M> &B, 
             const float tau,
             const float eps,
             const float sweeps_factor = 1)
{
    using matrix_nn = Eigen::Matrix<T, N, N>;
    using matrix_nm = Eigen::Matrix<T, N, M>;
    using matrix_mm = Eigen::Matrix<T, M, M>;
    using pivot     = std::tuple<T, int, int>;

    const int n = B.rows();
    const int m = B.cols();

    // тривиальная проверка на адекватность размеров матриц
    static_assert(N >= M || N == Eigen::Dynamic || M == Eigen::Dynamic,
                  "JTS_SVD_base Compile Error: N must be greater or equal M");
    assert((n >= m, "JTS_SVD_base Runtime Error: N must be greater or equal M"));

    // https://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#ac8da566526419f9742a6c471bbd87e0a
    const float delta = eps * B.squaredNorm();
    // https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
    matrix_mm V = matrix_mm::Identity(m, m);

    const int SWEEPS_MAX = (int64_t)m * (m - 1) / 2 * sweeps_factor + 1;
    for (int sweeps_cur = 0; sweeps_cur < SWEEPS_MAX; ++sweeps_cur)
    {
        // так сделано, чтобы можно было в конец вставлять
        std::vector<pivot> pivots;
        pivots.reserve(m * (m - 1) / 2);
        for (int i = 0; i < m - 1; ++i)
        {
            for (int j = i + 1; j < m; ++j)
            {
                // простое умножение падает
                // https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
                // в будущем это надо будет параллелить
                T b = std::abs(B.col(i).dot(B.col(j)));
                pivots.emplace_back(b, i, j);
            }
        }

        int count_left = pivots.size() * tau;
        count_left += count_left < pivots.size(); // чтобы хотя бы один нашелся
        // https://en.cppreference.com/w/cpp/algorithm/nth_element
        std::nth_element(pivots.begin(),
                         pivots.begin() + count_left,
                         pivots.end(),
                         std::greater<pivot>());
        pivots.resize(count_left);
        std::sort(pivots.begin(), pivots.end(), std::greater<pivot>());
        if (std::get<0>(pivots[0]) < delta) { break; }

        /*
         *  в статье предлагается сначала предлагается поместить тройки в очередь
         *  и только после того как они были забиты, производить подсчеты
         *  поскольку в очереди FIFO то можно по этому принципу сделать и без очереди
         */
        for (const auto &[ _, i, j ] : pivots)
        {                
            T gamma = (B.col(j).squaredNorm() - B.col(i).squaredNorm()) / (2 * B.col(i).dot(B.col(j)));
            T t = 1 / (std::abs(gamma) + std::sqrt(gamma * gamma + 1));
            t *= gamma != 0 ? gamma / std::abs(gamma) : 0; // для работы с комлексными
            T c = 1 / std::sqrt(1 + t * t);
            T s = t * c;

            // умножение на матрицу якоби, меняются только 2 столбца
            // матрицу нельзя менять in-place
            using col_b = Eigen::Matrix<T, 1, N>;            
            col_b col_bi = B.col(i);
            col_b col_bj = B.col(j);
            B.col(i) = col_bi * c - col_bj * s;
            B.col(j) = col_bi * s + col_bj * c;

            using col_v = Eigen::Matrix<T, 1, M>;
            col_v col_vi = V.col(i);
            col_v col_vj = V.col(j);
            V.col(i) = col_vi * c - col_vj * s;
            V.col(j) = col_vi * s + col_vj * c;
        }
    }

    std::vector<std::pair<T, int>> norms(m);
    T norm;
    for (int i = 0; i < m; ++i)
    {
        norm = B.col(i).norm();
        norms[i] = std::make_pair(norm, i);
        B.col(i) /= norm; // нормировка столбца заранее
    }
    std::sort(norms.begin(), norms.end(), std::greater<std::pair<T, int>>());
    
    Eigen::PermutationMatrix<M> P(m);
    matrix_nm S = matrix_nm::Zero(n, m);
    matrix_nn U = matrix_nn::Zero(n, n);
    for (int i = 0; i < m; ++i)
    {
        S(i, i) = norms[i].first;
        P.indices()[i] = norms[i].second;
        // P.indices()[norms[i].second] = i;
    }
    U.block(0, 0, n, m) = B * P;
    V = V * P;

    return std::make_tuple(U, S, V.transpose());
}

template<typename T, int N, int M>
std::tuple<Eigen::Matrix<T, N, N>,    // U
           Eigen::Matrix<T, N, M>,    // S
           Eigen::Matrix<T, M, M>>    // V
JTS_SVD(const Eigen::Matrix<T, N, M> &A,
        const float tau,
        const float eps,
        const float max_sweep_factor = 1)
{
    if (A.rows() < A.cols())
    {
        // как будто тут можно получше сделать
        Eigen::Matrix<T, M, N> B = A.transpose();
        auto [ U, S, V ] = JTS_SVD_base(B, tau, eps, max_sweep_factor);
        return std::make_tuple(V.transpose(), S.transpose(), U.transpose());
    }
    else
    {
        Eigen::Matrix<T, N, M> B = A;
        return JTS_SVD_base(B, tau, eps, max_sweep_factor);
    }
}


int main()
{
    using namespace std;
    using namespace Eigen;

    // PermutationMatrix<Dynamic> P(3);
    // // P.setIdentity();
    // for (int i = 0; i < 3; ++i)
    // {
    //     P.indices()[i] = 3 - i - 1;
    // }
    // auto P2 = P * P.transpose();
    // for (int i = 0; i < 3; i++)
    //     cout << P2.indices()[i] << "\n";

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(10, 9);
    m << 1,  2,  3,  4,  5,  6,  7,  8,  9,
         10, 11, 12, 13, 14, 15, 16, 17, 18,
         19, 20, 21, 22, 23, 24, 25, 26, 27,
         28, 29, 30, 31, 32, 33, 34, 35, 36,
         37, 38, 39, 40, 41, 42, 43, 44, 45,
         46, 47, 48, 49, 50, 51, 52, 53, 54,
         55, 56, 57, 58, 59, 60, 61, 62, 63,
         64, 65, 66, 67, 68, 68, 70, 71, 72,
         73, 74, 75, 76, 77, 78, 79, 80, 81,
         82, 83, 84, 85, 86, 87, 88, 89, 90;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m2 = m.transpose();
    auto [U, S, V] = JTS_SVD(m2, 0.1, 1e-9);

    cout << U << "\n\n" << S << "\n\n" << V << "\n\n";
    cout << U * S * V << "\n\n";
    
    for (int i = 0; i < U.cols(); i++)
    {
        for (int j = 0; j < U.cols(); j++)
        {
            cout << setprecision(4) << fixed << U.col(i).dot(U.col(j)) << "\t";
        }
        cout << "\n";
    }
    cout << "\n";
    cout << U * U.transpose();

    char c;
    cin >> c;
    return 0;
}