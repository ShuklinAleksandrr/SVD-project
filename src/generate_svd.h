//
//  main.cpp
//  SVD
//
//  Created by Victoria Koreshkova on 05.11.2023.
//
#include <iostream>
#include <Eigen/Dense>
#include "ctime"
#include "tuple"
#include <random>
#include <cassert>
//Получение матриц U,SIGMA,V^T по массиву сингулярных значений
//m - число строк
//n - число столбцов
//lmbd - массив сингулярных значений
template <typename T,int M = Eigen::Dynamic,  int N = Eigen::Dynamic>
class SVDGenerator
{
    private:

    using MatrixUType = Eigen::Matrix<T, M, M>;
    using MatrixVType = Eigen::Matrix<T, N, N>;
    using MatrixSType = Eigen::Matrix<T, M, N>;
    using DynamicMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
    using SingValVector =  std::vector<T>;

    bool generatedFLG = false;
    bool sigma_init_flg = false;
    MatrixUType U;
    MatrixVType V;
    MatrixSType S;
    SingValVector sigmas;  

    public:

    SVDGenerator(int rows, int cols)
    {
        assert(rows > 0);
        assert(cols > 0);

        U = MatrixUType::Zero(rows, rows);
        V = MatrixVType::Zero(cols, cols);
        S = MatrixSType::Zero(rows, cols);
    }

    SVDGenerator(int rows, int cols, const SingValVector& sigmas1)
    {
        assert(rows > 0);
        assert(cols > 0);

        U = MatrixUType::Zero(rows, rows);
        V = MatrixVType::Zero(cols, cols);
        S = MatrixSType::Zero(rows, cols);

        set_sing_vals(sigmas1);

    }

    void set_sing_vals(const SingValVector& sigmas1)
    {
        int K = U.rows() < V.rows() ? U.rows() : V.rows();

        assert(K == sigmas1.size());

        sigmas.clear();
        for(int i = 0; i < K; ++i)
        {
            sigmas.push_back(sigmas1[i]);
        }

        std::sort(sigmas.begin(), sigmas.end(), std::greater_equal<T>());
    }


    MatrixUType MatrixU()
    {
        if (!generatedFLG)
            generate();
        return U;
    }
    
    MatrixVType MatrixV()
    {
        if (!generatedFLG)
            generate();
        return V;
    }

    MatrixSType MatrixS()
    {
        if (!generatedFLG)
            generate();
        return S;
    }   

    void generate()
    {
        generatedFLG = true;
        
        //U,V-ортогональные матрицы, SIGMA - диагональная матрица
        
        /*Создается две случайные матрицы нужных размеров - T1 и T2,элементы - случайные числа от 0.1 до 10. QR разложение раскладывает матрицу на произведение двух: ортогональной Q и верхнедиагональной R
        С помощью этого разложения случайная матрица превращается в ортогональную и далее эта
        ортогональная матрица Q используется как U или V */

        int m = U.rows();
        int n = V.rows();
        
        DynamicMatrix T_1(m,m), T_2(n,n), Q_1(m,m), Q_2(n,n), R_1(m,m), R_2(n,n);
        //Сингулярные значения нумеруются в порядке убывания
        //Тут на всякий случай сортируем массив сингулярных чисел, чтобы элменты шли в порядке убывания l1 >= l2 >=.....>= lk >= 0

        S.setZero();
        
        for(int i = 0; i < sigmas.size(); i++)
            S(i,i) = sigmas[i];
            
        
        //Заполнение матриц T_1 и T_2 случайными элементами от 0.1 до 10
        T HI = static_cast<T>(10);
        T LO = static_cast<T>(0.1);
        T range= HI-LO;
        
        T_1 = MatrixUType::Random(U.rows(),U.rows());
        T_1 = (T_1 + MatrixUType::Constant(U.rows(),U.rows(),static_cast<T>(1)))*range/2;
        T_1 = (T_1 + MatrixUType::Constant(U.rows(),U.rows(),LO));
        
        T_2 = MatrixVType::Random(V.rows(),V.rows());
        T_2 = (T_2 + MatrixVType::Constant(V.rows(),V.rows(),static_cast<T>(1)))*range/2;
        T_2 = (T_2 + MatrixVType::Constant(V.rows(),V.rows(),LO));
        
        
        //QR_разложение матриц T1 и T2
        Q_1 = (Eigen::FullPivHouseholderQR<Eigen::Matrix<T, N, M>>(T_1)).matrixQ();
        Q_2 = (Eigen::FullPivHouseholderQR<Eigen::Matrix<T, N, M>>(T_2)).matrixQ();
        U = Q_1;
        V = Q_2.transpose();

    }
};
