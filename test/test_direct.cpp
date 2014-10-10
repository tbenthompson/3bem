#include "UnitTest++.h"
#include <iostream>
#include <math.h>
#include "test_shared.h"
#include "direct.h"

const double PI = 4.0 * std::atan(1.0);

TEST(Direct) {
    CHECK_CLOSE(one_kernel(0, 0, 0, 0, 0, 0), 1.0, 1e-14);

    int n = 8 * 250;
    auto src = random_pts(n);
    auto obs = random_pts(n);
    std::vector<double> values(n, 1.0);

    Kernel k = laplace_single;
    TIC
    auto result = direct_n_body(src, obs, k, values);
    TOC("Direct N Body");
    // std::cout << "number of interaction: " << n * n << std::endl;
    


    // std::vector<double> correct(n);
    // for (int i = 0; i < n; ++i) correct[i] = n;
    // CHECK_ARRAY_CLOSE(correct, result, n - 1, 1e-8);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
