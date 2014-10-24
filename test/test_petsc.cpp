#include "UnitTest++.h"
#include "petsc_interface.h"
#include "util.h"
#include <iostream>

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

// TEST(Solve2X2) {
//     std::vector<double> rhs = {random 
// }

int main(int, char const* args[])
{
    return UnitTest::RunAllTests();
}
