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
tuple<Eigen::Matrix<T, M, M>, Eigen::Matrix<T, M, N>, Eigen::Matrix<T, N, N>> matrices()
{
    tuple<Eigen::Matrix<T, M, M>, Eigen::Matrix<T, M, N>, Eigen::Matrix<T, N, N>> tuple_matrix;
    
    //U,V-ортогональные матрицы, SIGMA - диагональная матрица
    
    /*Создается две случайные матрицы нужных размеров - T1 и T2,элементы - случайные числа от 0.1 до 10. QR разложение раскладывает матрицу на произведение двух: ортогональной Q и верхнедиагональной R
    С помощью этого разложения случайная матрица превращается в ортогональную и далее эта
    ортогональная матрица Q используется как U или V */
    
    Matrix<T,Dynamic,Dynamic> T_1(M,M), T_2(N,N), Q_1(M,M), Q_2(N,N), R_1(M,M), R_2(N,N), U(M,M), V_T(N,N), SIGMA(M,N);
    
    int K = (M <= N) ? M : N;
    T lmbd[K]; // lmbd - массив сингулярных значений размера K
    
    //Сингулярные значения нумеруются в порядке убывания
    //Тут на всякий случай сортируем массив сингулярных чисел, чтобы элменты шли в порядке убывания l1 >= l2 >=.....>= lk >= 0
    
    Map<Vector<T,Dynamic>> l(lmbd,K);
    
    //Пока такая сортировка
    std::sort(l.begin(), l.end(), std::greater_equal<float>());
    
    
    //Заполнение матрицы сигма
    //На диагонали стоят сингулярные значения, остальные элементы - 0
    
    SIGMA.setZero();
    
    for(int i = 0; i < K; i++)
        SIGMA(i,i) = l[i];
        
    
    //Заполнение матриц T_1 и T_2 случайными элементами от 0.1 до 10
    float HI = 10;
    float LO = 0.1;
    float range= HI-LO;
    
    T_1 = Matrix<T,Dynamic,Dynamic>::Random(M,M);
    T_1 = (T_1 + Matrix<T,Dynamic,Dynamic>::Constant(M,M,1))*range/2;
    T_1 = (T_1 + Matrix<T,Dynamic,Dynamic>::Constant(M,M,LO));
    
    T_2 = Matrix<T,Dynamic,Dynamic>::Random(N,N);
    T_2 = (T_2 + Matrix<T,Dynamic,Dynamic>::Constant(N,N,1))*range/2;
    T_2 = (T_2 + Matrix<T,Dynamic,Dynamic>::Constant(N,N,LO));
    
    
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
    
    //k = (m <= n) ? m : n; //k=min(m,n), равенство для квадратной матрицы
    
    //Вызов функции пока закомментирован, так как нужны конкретные значения для m и n
    //matrices<float,m,n>();
    
    return 0;
}


