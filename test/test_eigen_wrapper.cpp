#include "UnitTest++.h"
#include "blas_wrapper.h"

using namespace tbem;

TEST(LU) {
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto lu = LU_decompose(matrix);
    auto soln = LU_solve(lu, {1,1});
    std::vector<double> correct{
        -0.25, 1.5
    };
    CHECK_ARRAY_CLOSE(soln, correct, 2, 1e-14);
}

TEST(LUCondition) {
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto lu = LU_decompose(matrix);
    double cond = condition_number(matrix);
    CHECK_CLOSE(cond, 2.7630857945186595, 1e-12);
}

int main() {
    return UnitTest::RunAllTests();
}
