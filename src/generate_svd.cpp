#include "generate_svd.h"

int main()
{

    std::vector<double> src;

    src.push_back(1.1);
    src.push_back(1.09);
    src.push_back(1.08);
    src.push_back(1);

    SVDGenerator<double> svd_gen(5, 4, src);

    std::cout << svd_gen.MatrixU();

    return 0;
}
