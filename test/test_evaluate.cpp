#include "UnitTest++.h"
#include <iostream>
#include <math.h>
#include "test_shared.h"
#include "direct.h"

const double PI = 4.0 * std::atan(1.0);

inline double one_kernel(double ox, double oy, double oz,
                  double sx, double sy, double sz) {
    return 1.0;
}

inline double laplace_single(double ox, double oy, double oz,
                             double sx, double sy, double sz) {
    const double r = sqrt((ox - sx) * (ox - sx) + 
                          (oy - sy) * (oy - sy) +
                          (oz - sz) * (oz - sz));
    return 1.0 / (4.0 * PI * r);
}

TEST(Direct) {
    CHECK_CLOSE(one_kernel(0, 0, 0, 0, 0, 0), 1.0, 1e-14);

    int n = (int)5e2;
    std::array<std::vector<double>,3> src =
        {random_list(n), random_list(n), random_list(n)};
    std::array<std::vector<double>,3> obs =
        {random_list(n - 1), random_list(n - 1), random_list(n - 1)};
    std::vector<double> values(n);
    for (int i = 0; i < n; ++i) values[i] = 1.0;

    TIC
    auto result = direct_n_body(src, obs, laplace_single, values);
    TOC("Direct N Body");


    // std::vector<double> correct(n);
    // for (int i = 0; i < n; ++i) correct[i] = n;
    // CHECK_ARRAY_CLOSE(correct, result, n - 1, 1e-8);
}
