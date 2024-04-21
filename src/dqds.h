// guardspirit@protonmail.com
// КМБО-03-22

#include <lapacke.h>
#include <Eigen/Core>
#include <Eigen/SVD>

// Класс, реализующий DQDS с использованием LAPACKE.
// (LAPACKE — обёртка над LAPACK)
// Наследует Eigen::SVDBase.
template<typename _MatrixType> class DQDS_SVD
 : public Eigen::SVDBase<DQDS_SVD<_MatrixType> >
{
private:
    typedef Eigen::SVDBase<DQDS_SVD<_MatrixType>> Base;
public:
    // Конструктор, выполняющий вычисление U, сингулярных значений и V
    DQDS_SVD (const _MatrixType &matrix, unsigned int computationOptions) {
        compute(matrix);
        // TODO: computationOptions на данный момент игнорируются.
    }

    DQDS_SVD<_MatrixType>& compute(const _MatrixType& matrix) {
        auto bid = Eigen::internal::UpperBidiagonalization(matrix); // DQDS принимает на вход бидиагональную матрицу
        typename Eigen::internal::UpperBidiagonalization<_MatrixType>::BidiagonalType bid_matrix = bid.bidiagonal();
        // Копируем итоговую бидиагональную матрицу
        auto diagonal = bid_matrix.diagonal(0); // Главная диагональ. diagonal.data() содержит массив только из элементов диагонали
        auto offdiagonal = bid_matrix.diagonal(1); // Наддиагональ
        this->m_matrixU = bid.householderU(); // U. diagonal.data() содержит U в формате Column Major
        this->m_matrixV = bid.householderV(); // V. diagonal.data() содержит V в формате Row Major => при интерпретировании данных
        // как Column Major diagonal.data() содержит V^T.
        int n = bid.bidiagonal().cols();
        int info = LAPACKE_dbdsqr(
                LAPACK_COL_MAJOR,
                'U', // Upper bidiagonal
                n, // N: Order of B
                this->m_matrixV.cols(), // NCVT: Columns in VT
                this->m_matrixU.rows(), // NRU: Rows in U
                0, // NCC: Columns in C
                diagonal.data(), // [in, out]
                offdiagonal.data(), // [in, out]
                this->m_matrixV.data(),
                this->m_matrixV.rows(), // LDVT: Leading dimension of VT
                // т. е. как много элементов необходимо пропустить LAPACK'у, чтобы добраться до начала следующей колонки
                this->m_matrixU.data(), // [in, out]
                this->m_matrixU.rows(), // LDU
                nullptr, // C
                1
        );
        if (info != 0)
            throw std::runtime_error("LAPACK error: " + std::to_string(info));
        this->m_singularValues = diagonal;
        this->m_isInitialized = true;
        this->m_computeFullU = true; // TODO
        this->m_computeFullV = true;
        return *this; 
    }
};

// 
template<typename _MatrixType> 
struct Eigen::internal::traits<DQDS_SVD<_MatrixType>>
        : Eigen::internal::traits<_MatrixType>
{
    typedef _MatrixType MatrixType;
};
