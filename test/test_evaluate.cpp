#include "UnitTest++.h"
#include <iostream>
#include <math.h>
#include "test_shared.h"
#include "direct.h"
#include "cheb.h"

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
    if (r < 0.4) {
        return 0.0;
    }
    return 1.0 / (4.0 * PI * r);
}

TEST(ChebPolys) {
    std::vector<double> theta({0.0, 0.25, 0.5, 0.75});
    for (int i = 0; i < 4; i++) {
        std::vector<double> exact(11);
        for (int j = 0; j < 11; j++) {
            exact[j] = cos(theta[i] * j);
        }
        auto res = cheb_polys(exact[1], 10);
        CHECK_ARRAY_CLOSE(res, exact, 11, 1e-12);
    }
}

TEST(ChebPtsFirstKind) {
    auto pts = cheb_pts_first_kind(4);
    CHECK_CLOSE(pts[0], 0.92388, 1e-4);
}
TEST(ChebPtsSecondKind) {
    auto pts = cheb_pts_second_kind(4);
    CHECK_CLOSE(pts[0], 1.0, 1e-4);
    CHECK_CLOSE(pts[1], 0.5, 1e-4);
}

TEST(ChebInterp) {
    int n_pts_interp = 32;
    int n_test_locs = 300;
    auto fnc = [](double x){return std::sin(16 * x);};
    std::vector<double> est(n_test_locs);
    std::vector<double> exact(n_test_locs);
    auto cheb_nodes = cheb_pts_first_kind(n_pts_interp);
    for(int i = 0; i < n_test_locs; i++) {
        double loc = -1.0 + (2.0 / (n_test_locs - 1)) * i;
        for (int m = 0; m < n_pts_interp; m++) {
            est[i] += fnc(cheb_nodes[m]) * s_n(cheb_nodes[m], loc, n_pts_interp);
        }
        exact[i] = fnc(loc);
    }
    CHECK_ARRAY_CLOSE(est, exact, n_test_locs, 1e-5);
}

TEST(Direct) {
    CHECK_CLOSE(one_kernel(0, 0, 0, 0, 0, 0), 1.0, 1e-14);

    int n = (int)5e3;
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
