#include "UnitTest++.h"
#include <iostream>
#include <cmath>
#include "util.h"
#include "direct.h"

using namespace tbem;

const double PI = 4.0 * std::atan(1.0);

TEST(Direct) {
    UNITTEST_TIME_CONSTRAINT(50);
    CHECK_CLOSE(one_kernel(0, 0, 0, 0, 0, 0), 1.0, 1e-14);

    int n = 8 * 25 * 2;
    auto src = random_pts(n);
    auto obs = random_pts(n);
    std::vector<double> values(n, 1.0);

    Kernel k = laplace_single;
    auto result = direct_n_body(src, obs, k, values);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
