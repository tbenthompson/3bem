#include "UnitTest++.h"
#include "petsc_facade.h"
#include "util.h"
#include <iostream>

using namespace tbem;

TEST(SolveSystem) {
    int n = 10;
    int divide = 3.0;
    std::vector<double> rhs(n, 1.0);

    auto res = solve_system(rhs, 1e-3, 
        [&](std::vector<double>& x, std::vector<double>& y) {
            for(unsigned int i = 0; i < x.size(); i++) {
                y[i] = x[i] * divide;
            }
        });

    std::vector<double> exact(n, 1.0 / divide);
    CHECK_ARRAY_CLOSE(res, exact, n, 1e-12);
}

TEST(Solve2X2) {
    std::vector<double> rhs = random_list(2);
    std::vector<double> row1 = random_list(2);
    std::vector<double> row2 = random_list(2);
    std::array<std::vector<double>,2> matrix = {row1,row2};

    auto res = solve_system(rhs, 1e-12, 
        [&](std::vector<double>& x, std::vector<double>& y) {
            for(unsigned int i = 0; i < 2; i++) {
                y[i] = 0;
                for(unsigned int j = 0; j < 2; j++) {
                    y[i] += matrix[i][j] * x[j];
                }
            }
        });

    double denom = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]);
    std::vector<double> exact = {
        (rhs[0] * matrix[1][1] - rhs[1] * matrix[0][1]) / denom,
        (rhs[1] * matrix[0][0] - rhs[0] * matrix[1][0]) / denom,
    };
    CHECK_ARRAY_CLOSE(res, exact, 2, 1e-5);
}

int main(int, char const* args[])
{
    return UnitTest::RunAllTests();
}
