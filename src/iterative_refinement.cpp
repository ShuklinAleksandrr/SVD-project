#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <vector>
#include <Eigen/LU> 
#include <Eigen/SVD>





/*Название статьи:  Uchino Yuki, Terao Takeshi, Ozaki Katsuhisa. 2022.08.05 – Acceleration of 
Iterative Refinement for Singular Value Decomposition.
(https://www.researchgate.net/publication/362642883_Acceleration_of_Iterative_Refinement_for_Singular_Value_Decomposition)
* Обязательные условия: 1)в матрицах колво строк больше или равному количеству столбцов
*2)каждое сингулярное значение больше последующего (d1>d2>...>dn)
* 3)каждое приближенное сингулярное значение отлично друг от друга (di != dj для i!=j)
* 
* Выполнил Едренников Д.А. КМБО-01-20
*/



template<typename T>
 std::tuple<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
MSVD_SVD(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, const  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &Ui, const  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &Vi)
{
    const int m = A.rows();//Получаем размеры изначальной матрицы
    const int n = A.cols();


    if ( A.rows()< A.cols())
    {
        std::cout << "Attention! Number of the rows must be greater or equal  than  number of the columns";
        return std::make_tuple(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Constant(m, m, (T)0.0), Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Constant(n, n, (T)0.0),
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Constant(m, n, (T)0.0));
    }

   using matrix_dd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


    matrix_dd Ad = A.cast<double>();//Согласно нашему алгоритму все почти все вычисления должны 
    matrix_dd Ud = Ui.cast<double>();//производится с двойной точностью, поэтому переделываем тип
    matrix_dd Vd = Vi.cast<double>();//всех матриц в double
    //Далее все тип всех элементов также будет double. Но в конце, все матрицу будут приведены к
    //изначальному типу.

    matrix_dd Udmn = Ud.block(0, 0, m, n);


    matrix_dd P = Ad * Vd;
    matrix_dd Q = Ad.transpose() * Udmn;

    matrix_dd ViT = Vd.transpose();
    matrix_dd UiT = Ud.transpose();
    std::vector<double> r (n, 0.0);
    std::vector<double> s (n, 0.0);
    std::vector<double> t (n, 0.0);
    matrix_dd Sigma_n(n,n);
    Sigma_n.setZero();

    for (int i = 0; i < n; i++)
    {
        r[i] = 1.0 - UiT.row(i) * Ud.col(i);
        s[i] = 1.0 - ViT.row(i) * Vd.col(i);
        t[i] = UiT.row(i) *P.col(i);
        Sigma_n(i,i) = t[i] / (1 - ((r[i] + s[i]) * 0.5));

    };

    matrix_dd Sigma(m,n);//Матрица сингулярных значений
    Sigma.setZero();
    Sigma.block(0,0,n, n) = Sigma_n.transpose();

    matrix_dd Cgamma = P - Udmn * Sigma_n;
    matrix_dd Cdelta = Q - Vd * Sigma_n;

    matrix_dd Calpha = Udmn.transpose() * Cgamma;
    matrix_dd Cbetta = ViT * Cdelta;

    matrix_dd D = Sigma_n * Calpha + Cbetta * Sigma_n;
    matrix_dd E = Calpha * Sigma_n + Sigma_n * Cbetta;

    matrix_dd G(n,n);//Матрица, необходимая для более точного вычисления правых сингулярных векторов
    matrix_dd F11(n, n);

    double temp1;//Временные переменные, необходимые для более быстрого вычисления элементов матриц
    double temp2;//G и F11

    for (int i = 0; i < n; i++)//Вычисление элементов матриц G и F11
    {
        temp1 = Sigma_n(i, i) * Sigma_n(i, i);

        for (int j = 0; j < n; j++)
        {
            temp2 = Sigma_n(j, j) * Sigma_n(j, j);

            G(i, j) = D(i, j) / (temp2 - temp1);//Значение диагональных элементов не должно равняться 
            F11(i, j) = E(i, j) / (temp2 - temp1);//бесконечности. Оно будет пересчитано
        }
    }

    for (int i = 0; i < n; i++)//Пересчет диагональных элементов для матриц G и F11
    {
        G(i, i) = s[i] * 0.5;
        F11(i, i) = r[i] * 0.5;
    }


    matrix_dd F12(n, m-n);
    F12 = Sigma_n.inverse() * P.transpose() * Ud.block(0, n, m, m - n);
    matrix_dd F21 = Ud.block(0,n-1,m, m - n).transpose() * Cgamma * Sigma_n.inverse();
    matrix_dd F22 = 0.5 * (matrix_dd::Constant(m - n, m - n, 1.0) -
        Ud.block(0,n-1,m, m - n).transpose() * Ud.block(0, n, m, m - n));


    
    matrix_dd F(m,m); //Матрица, необходимая для более точного вычисления левых сингулярных векторов
    F.block(0, 0, n, n) = F11;
    F.block(0, n , n, m-n) = F12;
    F.block(n, 0, m-n, n) = F21;
    F.block(n , n, m - n, m - n) = F22;


    matrix_dd U = Ud+ Ud*F;//Вычисление уточнённых значений левых сингулярных векторов
    matrix_dd V = Vd + Vd * G;//Вычисление уточнённых значений правых сингулярных векторов


   return std::make_tuple(U.cast<T>(), V.cast<T>(), Sigma.cast<T>()); // Приведение матриц к изначальному типу
}


template<typename T>
inline std::tuple<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
Accurate_BDCSVD(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A)
{
    Eigen::BDCSVD< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Неточный расчёт левых и правых сингулярных векторов с помощью функции из библиотеки Eigen
    return  MSVD_SVD(A, svd.matrixU(), svd.matrixV());// Уточнение результата нашей функцией 
}

int main()
{
    using namespace std;
    using namespace Eigen;
    Eigen::Matrix<float, Dynamic, Dynamic> A(10,9);
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25, 26, 27,
        28, 29, 30, 31, 32, 33, 34, 35, 36,
        37, 38, 39, 40, 41, 42, 43, 44, 45,
        46, 47, 48, 49, 50, 51, 52, 53, 54,
        55, 56, 57, 58, 59, 60, 61, 62, 63,
        64, 65, 66, 67, 68, 68, 70, 71, 72,
        73, 74, 75, 76, 77, 78, 79, 80, 81,
        3, 9, (float)4.98942, (float)0.324235,  443534, 345, (float)56.543853, (float)450.435234, (float)43.34353221;

   BDCSVD<MatrixXf> svd(A, ComputeFullU | ComputeFullV);
   std::tuple<Eigen::Matrix<float, 10, 10>, Eigen::Matrix<float, 9, 9>, Eigen::Matrix<float, 10, 9>> a =
       Accurate_BDCSVD(A);
   std::cout << get<0>(a)*get<2>(a)*get<1>(a).transpose() <<"\n" << "\n";
    Array<float,1, Dynamic> sigm(9);
    sigm  = svd.singularValues();
    Eigen::Matrix<float, Dynamic, Dynamic > I(10, 9);
    I.setZero();
    I.block(0,0,9,9) = sigm.matrix().asDiagonal();
    Eigen::Matrix<float, Dynamic, Dynamic > U = svd.matrixU();
    cout << U*I* svd.matrixV().transpose() << "\n" << "\n";

    return 0;
}

