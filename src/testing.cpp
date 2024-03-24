//Александр Нам, КМБО-04-20

#include "generate_svd.h"
#include "jacobi.h"

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cassert>

/*
 Функция возвращает матрицу-таблицу формата:
 | размерность | sigma_max/sigma_min | диап. синг. чисел | 
 | sum_n(norm[I - U.transpose*U])/n | sum_n(norm[I - U*U.transpose])/n |
 | sum_n(norm[I - V.transpose*V])/n | sum_n(norm[I - V*V.transpose])/n | 
 | max(abs(sigma_true_i - sigma_calc_i)) |;
 Размер выборки фиксированного размера матриц определяется задаваемым параметром 'n'.
*/

/*
Данный код писался, опираясь на generate_svd.h и jacobi.h, но не исключительно ПОД них, поэтому:
Функция написана с теми условиями, что:
    1. Нешаблонный функтор SVD-разложения, но с шаблонными свойствами, принимает 
       матрицу Eigen::Matrix<T, N, M> и возвращает
       std::tuple<Eigen::Matrix<T, N, N>, Eigen::Matrix<T, N, M>, Eigen::Matrix<T, M, M>>
            1.1. Важно, что возвращаемые функцией матрицы - это U, S, V_t. 
    2. Генерация матриц происходит методом класса, который передаётся в аргумент параметра шаблона 
            2.1. Класс должен содержать конструктор, принимающий в аргументы int размеры матрицы и std::vector синг. чисел
            2.2. Класс должен содержать метод generate(), генерирующий матрицы U, S, V
            2.3. Класс должен содержать геттеры MatrixU(), MatrixS(), MatrixV()
    3. В функцию передаётся std::vector соотношений максимумального и минимального сингулярного числа
    4. В функцию передаётся std::vector<std::pair<int,int>> размеров матриц для исследования
    5. В функцию передаётся int n размер выборки фиксированного размера матриц для подсчёта средних
    6. Функция работает достаточно долго, особенно для матриц больших размеров, поэтому выводится прогресс в процентах
    7. Результат исследования не печатается в консоль, а сохраняется в файл с названием fileName
*/

template<typename T, template <typename> class Gen_Class, class SVD_Class> 
void CMP_SVD_Method(std::string fileName, SVD_Class&& svd_functor,const std::vector<T>& SigmaMaxMinRatiosVec, const std::vector<std::pair<int,int>>& MatSizesVec, const int n){
/*
Алгоритм функции:
Для каждого размера из вектора размеров матриц:
    Для каждого соотношения синг. чисел из вектора отношений синг. чисел:
        Для каждого диапазона{[0,1] и [1,100]}:
            1. Генерация сингулярных чисел sigma_i из диапазона с заданным соотношением sigma_max/sigma_min.
            2. Сделать n раз:
                2.1. Генерация матрицы A с заранее заданными синг. числами.
                2.2. Восстановить исходную матрицу A путем перемножения, т.е. A=USV*.
                2.3. Использовать метод на матрице A для получения U_calc, S_calc,  V_calc_transpose.
                2.4. Вычислить 4 квадратичных нормы(мера унитарности) матриц U_calc и V_calc, 
                     разделить на n и суммировать в аккумуляторы.
                2.5. Вычислить max(|sigma_true_i-sigma_jacobi_i|).
            3. Занести в матрицу-таблицу соответствующую строку.
*/

    auto printTable = [](std::ostream& out, const std::vector<std::vector<std::string>>& data){
        if (data.empty()) return;
        // Определение ширины столбцов
        std::vector<size_t> widths;
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                if (widths.size() <= i) widths.push_back(0);
                widths[i] = std::max(widths[i], row[i].size());
            }
        }
        // Вывод данных в поток
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                out << std::left << std::setw(widths[i] + 3) << row[i];
            }
            out << "\n";
        }
    };

    auto num2str = [](T value){
        std::ostringstream oss;
        oss << value;
        return oss.str();
    };

    //нужно для того, чтобы более менее равномерно выводить прогресс работы
    T generalSum=0;
    for (const auto& MatSize : MatSizesVec){
        generalSum += (MatSize.first*MatSize.second);
    } 
    
    using MatrixDynamic = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const std::vector<std::pair<T,T>> Intervals = {{0,1}, {1,100}};
    std::random_device rd; //случайное число - значение сида
    std::mt19937 gen(rd()); //генерация последовательности псевдослучайных цифр

    //название столбцов таблицы
    std::vector<std::vector<std::string>> table = {{"Dimension", "Sigma-max/min-ratio", "SV interval" , 
                                                   "AVG ||I-U_t*U||", "AVG ||I-U*U_t||", "AVG ||I-V_t*V||",
                                                   "AVG ||I-V*V_t||", "AVG sigma dev."}};

    T ProgressCoeff = n*Intervals.size()*SigmaMaxMinRatiosVec.size()*generalSum/100.0; //наибольшее значение прогресса
    T progress=0; //инициализация аккумулятора прогресса

    std::vector<T> sv_true_vec; //инициализация вектора задаеваемых сингулярных значений

    //инициализация матриц
    MatrixDynamic U_true, S_true, V_true, U_calc, S_calc, V_calc, V_calc_transpose;

    for (const auto& MatSize : MatSizesVec){
        const int N = MatSize.first;
        const int M = MatSize.second;
        int minNM = std::min(N, M);

        U_true.resize(N,N); U_calc.resize(N,N);
        S_true.resize(N,M); S_calc.resize(N,M);
        V_true.resize(M,M); V_calc.resize(M,M); V_calc_transpose.resize(M,M);     

        for (const auto& SigmaMaxMinRatio : SigmaMaxMinRatiosVec){
            for (const auto& interval : Intervals){
                //interval.first и interval.second - l и r границы интервала
                assert((interval.first<interval.second, "\nError: left boundary >= right boundary\n"));

                //пункт 1.
                //определение промежутка [sigma_min, sigma_max]:
                /*
                    проверка на то, что отношение sigmaMax/sigmaMin <= b/a, где a и b - границы интервала,
                    очевидно, что если неравенство не выполняется, то сингулярных чисел с заданным соотношением
                    и в заданном интервале не существует
                */
                assert((interval.first*SigmaMaxMinRatio <= interval.second, 
                       "\nError: no sigma values exist with such ratio in such interval\n"));

                /*
                  случайным образом генерируем sigma_min в допустимом интервале,
                  заметим, что прошлый assert гарантирует, что a < b/SigmaRatio
                */
                std::uniform_real_distribution<T> distrSigmaMin(interval.first, interval.second/SigmaMaxMinRatio);
                T sigma_min=distrSigmaMin(gen); 
                T sigma_max=SigmaMaxMinRatio*sigma_min;

                //заполнение вектора истинных сингулярных чисел из уже определенного диапазона:
                std::uniform_real_distribution<T> distr(sigma_min, sigma_max);
                assert((minNM >= 2, "\nError: no columns and rows allowed\n"));
                sv_true_vec.clear(); sv_true_vec.reserve(minNM);
                sv_true_vec.emplace_back(sigma_min);
                for (int i=0; i<minNM-2; ++i){
                    sv_true_vec.emplace_back(distr(gen));
                }
                sv_true_vec.emplace_back(sigma_max);

                //инициализация аккумуляторов выводимых в таблицу значений характеристик
                T avg_dev_UUt=0, avg_dev_UtU=0, avg_dev_VVt=0, avg_dev_VtV=0, avg_dev_sigma=0;

                //пункты 2.1. - 2.5.
                for (int i=1; i<=n; ++i){
                    Gen_Class svd_gen(N, M, sv_true_vec);  
                    svd_gen.generate();
                    U_true = svd_gen.MatrixU(); S_true = svd_gen.MatrixS(); V_true = svd_gen.MatrixV();
                    std::tie(U_calc, S_calc, V_calc_transpose) = svd_functor((U_true*S_true*V_true.transpose()).eval());
                    V_calc = V_calc_transpose.transpose();
                    avg_dev_UUt += (MatrixDynamic::Identity(N,N) - U_calc*U_calc.transpose()).squaredNorm()/n;
                    avg_dev_UtU += (MatrixDynamic::Identity(N,N) - U_calc.transpose()*U_calc).squaredNorm()/n;
                    avg_dev_VVt += (MatrixDynamic::Identity(M,M) - V_calc*V_calc.transpose()).squaredNorm()/n;
                    avg_dev_VtV += (MatrixDynamic::Identity(M,M) - V_calc.transpose()*V_calc).squaredNorm()/n;
                    avg_dev_sigma += (S_true.diagonal() - S_calc.diagonal()).cwiseAbs().maxCoeff()/n;

                    //вывод прогресса
                    progress += M*N/ProgressCoeff; std::cout << "\n"+num2str(progress)+"% was done.\n";
                }

                //заполнение одной строки таблицы
                table.emplace_back(std::vector<std::string>{num2str(N)+"x"+num2str(M), num2str(SigmaMaxMinRatio), 
                                   "["+num2str(interval.first)+", "+num2str(interval.second)+"]",
                                   num2str(avg_dev_UUt), num2str(avg_dev_UtU), 
                                   num2str(avg_dev_VVt), num2str(avg_dev_VtV), 
                                   num2str(avg_dev_sigma)});
            }
        }
    }

    //сохранение таблицы в файл, в консоль такая громоздкая таблица просто не выводится в нормальном виде
    std::ofstream file(fileName);
    if (file) {
        printTable(file, table);
        file.close();
    } else {
        std::cerr << "Error while creating/opening file!\n";
    }
};

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    CMP_SVD_Method<double, SVDGenerator>("jacobi_test_table.txt",
                                 JTS_SVD(0.1, 1e-9),
                                 {1.01, 1.2, 2, 5, 10, 50},
                                 {{3,3}, {5,5}, {10,10}, {20,20}, {50,50}, {100,100}}, 
                                 20); 
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationGlobal = end - start;
    std::cout << "\nFull execution time = " << durationGlobal.count() << " seconds.\n";

    char c; std::cin >> c;
    return 0;
}
