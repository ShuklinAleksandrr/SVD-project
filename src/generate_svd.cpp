//
//  main.cpp
//  SVD
//
//  Created by Victoria Koreshkova on 05.11.2023.
//

#include <iostream>
#include <Eigen/DENSE>
#include "ctime"
#include "tuple"
using namespace std;
using namespace Eigen;

//Получение матриц U,SIGMA,V^T по массиву сингулярных значений
//m - число строк
//n - число столбцов
//lmbd - массив сингулярных значений
template <typename T,int M, int N>
tuple<Eigen::Matrix<T, M, M>, Eigen::Matrix<T, M, N>, Eigen::Matrix<T, N, N>> generate_svd()
{
    tuple<Eigen::Matrix<T, M, M>, Eigen::Matrix<T, M, N>, Eigen::Matrix<T, N, N>> tuple_matrix;
    
    int K = (M <= N) ? M : N;
    
    //U,V-ортогональные матрицы, SIGMA - диагональная матрица
    
    /*Создается две случайные матрицы нужных размеров - T1 и T2,элементы - случайные числа от 0.1 до 10. QR разложение раскладывает матрицу на произведение двух: ортогональной Q и верхнедиагональной R
    С помощью этого разложения случайная матрица превращается в ортогональную и далее эта
    ортогональная матрица Q используется как U или V */
    
    Matrix<T,Dynamic,Dynamic> T_1(M,M), T_2(N,N), Q_1(M,M), Q_2(N,N), R_1(M,M), R_2(N,N), U(M,M), V_T(N,N), SIGMA(M,N);
    
    T lmbd[K]; // lmbd - массив сингулярных значений размера K
    
    //Сингулярные значения нумеруются в порядке убывания
    //Тут на всякий случай сортируем массив сингулярных чисел, чтобы элменты шли в порядке убывания l1 >= l2 >=.....>= lk >= 0
    
    Map<Vector<T,((M <= N) ? M : N)>> l(lmbd,K);
    
    //Пока такая сортировка
    std::sort(l.begin(), l.end(), std::greater_equal<float>());

    //Заполнение матрицы сигма
    //На диагонали стоят сингулярные значения, остальные элементы - 0
    if (M<N)
    {
        SIGMA << Matrix<T,Dynamic,Dynamic>(l.asDiagonal()),
        MatrixXf::Zero(K,N-K);
    }
    
    else
        
        if (M>N)
        {
            SIGMA << Matrix<T,Dynamic,Dynamic>(l.asDiagonal()),
            MatrixXf::Zero(M-K,N);
        }
    
       else //Квадратная матрица
       {
           SIGMA << Matrix<T,Dynamic,Dynamic>(l.asDiagonal());
       }
    
    T rnd_1[M*M];
    
    for (int i {0}; i<M*M; i++)
        rnd_1[i] = (rand() % 101 + 1)/ 10.0;
    
    T rnd_2[N*N];
    
    for (int i {0}; i<N*N; i++)
       rnd_2[i] = (rand() % 101 + 1)/ 10.0;


    Map<Matrix<T,Dynamic,Dynamic>> t_1 (rnd_1,M,M);
    Map<Matrix<T,Dynamic,Dynamic>>  t_2 (rnd_2,N,N);
    
    T_1 = t_1;
    T_2 = t_2;
    
    //QR_разложение матриц T1 и T2
    Q_1 = (FullPivHouseholderQR<MatrixXf>(T_1)).matrixQ();
    Q_2 = (FullPivHouseholderQR<MatrixXf>(T_2)).matrixQ();
    U = Q_1;
    V_T = Q_2;
    
    tuple_matrix = make_tuple(U,SIGMA,V_T);
    
    return tuple_matrix;

}

int main(int argc, const char * argv[]) {
    
    //int m, n, k;
    
    // k = (m <= n) ? m : n; //k=min(m,n), равенство для квадратной матрицы
    
    //Вызов функции пока закомментирован, так как нужны конкретные значения для m и n
    generate_svd<float,4,3>();
    
    return 0;
}
