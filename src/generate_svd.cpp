#include "generate_svd.h"

int main()
{

    SVDGenerator<double> SVD(4, 4, generate_standard_distribution<double>);
    SVD.generate(2);
    
    std::cout << SVD.MatrixU() * SVD.MatrixU().transpose();

    return 0;
}
