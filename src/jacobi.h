#ifndef JACOBI_H
#define JACOBI_H

/* 
 * Ишутин Андрей (Ishutin Andrei)
 * telegram: @looks_amazing
 * Github: Andy-823 
 */

#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>

// не у всех эйген одинаково подключается
// но этот костыль не везде заработает
// https://stackoverflow.com/questions/142877/can-the-c-preprocessor-be-used-to-tell-if-a-file-exists
// #if __has_include(<eigen/core>)
// # include <eigen/core>
// #else
#include <Eigen/Core>
// #endif

/*
 * Пусть пока все будет в одном хедере
 * это класс-обертка
 */
class JTS_SVD
{
  protected:
    template<typename T, int N, int M>
    std::tuple<Eigen::Matrix<T, N, N>,    // U
               Eigen::Matrix<T, N, M>,    // S
               Eigen::Matrix<T, M, M>>    // V
    _get_SVD(Eigen::Matrix<T, N, M> &B);

  public:
    JTS_SVD() {}
    JTS_SVD(float tau, float eps, float sweeps_factor = 1)
        : tau(tau), eps(eps), sweeps_factor(sweeps_factor) {}

    template<typename T, int N, int M>
    std::tuple<Eigen::Matrix<T, N, N>,    // U
               Eigen::Matrix<T, N, M>,    // S
               Eigen::Matrix<T, M, M>>    // V
    operator()(const Eigen::Matrix<T, N, M> &A);

    float tau = 0.1;
    float eps = 1e-9;
    float sweeps_factor = 1;
};

/*
 *  строк в матрице не больше чем столбцов
 *  статья, хорошей ссылки не нащел: On Speeding-up Parallel Jacobi Iterations for SVD
 *  параллельности пока нет
 *  про матрицы https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
 *  возвращаемый тип такой, поскольку матрицы имеют разный размер
 *  нужны такие матрицы:
 *      Eigen::Matrix<T, N, N> U
 *      Eigen::Matrix<T, N, M> S
 *      Eigen::Matrix<T, M, M> V
 *  видно различие типов, поэтому взят std::tuple
 *
 *  TODO: сделать что-то получше тюпла
 */
template<typename T, int N, int M>
std::tuple<Eigen::Matrix<T, N, N>,    // U
           Eigen::Matrix<T, N, M>,    // S
           Eigen::Matrix<T, M, M>>    // V
JTS_SVD::_get_SVD(Eigen::Matrix<T, N, M> &B)
{
    using matrix_nn = Eigen::Matrix<T, N, N>;
    using matrix_nm = Eigen::Matrix<T, N, M>;
    using matrix_mm = Eigen::Matrix<T, M, M>;
    using pivot     = std::tuple<T, int, int>;

    const int n = B.rows();
    const int m = B.cols();

    // тривиальная проверка на адекватность размеров матриц
    // на всякий случай
    static_assert(N >= M || N == Eigen::Dynamic || M == Eigen::Dynamic,
                  "JTS_SVD Compile Error: N must be greater or equalthan M");
    assert((n >= m, "JTS_SVD Runtime Error: N must be greater or equal than M"));

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
            T gamma = (B.col(j).squaredNorm() - B.col(i).squaredNorm()) 
                      / (2 * B.col(i).dot(B.col(j)));
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
    /*
     * U должна состоять из ортогональных столбцов
     * но здесь, на самом деле, это не обязательно получится
     * пока было решено оставить так, как есть
     */
    U.block(0, 0, n, m) = B * P;
    V = V * P;

    return std::make_tuple(U, S, V.transpose());
}

template<typename T, int N, int M>
std::tuple<Eigen::Matrix<T, N, N>,    // U
           Eigen::Matrix<T, N, M>,    // S
           Eigen::Matrix<T, M, M>>    // V
JTS_SVD::operator()(const Eigen::Matrix<T, N, M> &A)
{
    if (A.rows() < A.cols())
    {
        // как будто тут можно получше сделать
        // но самый примитивный вариант не скомпилировался
        Eigen::Matrix<T, M, N> B = A.transpose();
        auto [ U, S, V ] = _get_SVD(B);
        return std::make_tuple(V.transpose(), S.transpose(), U.transpose());
    }
    else
    {
        Eigen::Matrix<T, N, M> B = A;
        return _get_SVD(B);
    }
}

#endif // JACOBI_H
