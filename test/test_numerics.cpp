#include "UnitTest++.h"
#include "numerics.h"
#include <iostream>
#include "test_shared.h"

TEST(LinearMapping) {
    int n = 10;
    auto vals = random_list(n);
    for (auto v: vals) {
        const double ref = real_to_ref(v, vals[0], vals[1]);
        const double real = ref_to_real(ref, vals[0], vals[1]);
        // std::cout << v << " " << ref << " " << real << " " << 
        //              vals[0] << " " << vals[1] << std::endl;
        CHECK_CLOSE(real, v, 1e-14);
    }
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

TEST(TanhSinh) {
    auto qr = double_exp(8, 0.3); 
    double result = integrate(qr, [](double x) {return std::pow(x, 10) * std::log(x + 1);});
    CHECK_CLOSE(-0.215466, result, 1e-4);

    result = integrate(qr, [](double x) {return std::pow(x, 4) * std::log(x + 1);});
    CHECK_CLOSE(-0.336074, result, 1e-5);

    result = integrate(qr, [](double x) {return std::log(x + 1);});
    CHECK_CLOSE(-0.613706, result, 1e-6);
}

TEST(GaussQuadrature) {
    auto qr = gauss(4);
    CHECK(qr.size() == 4);

    //Check that the weights are in the right range
    for (auto xw: qr) {
        CHECK(std::fabs(xw.first) <= 1);
        // std::cout << "x: " << xw.first << "   w: " << xw.second << std::endl;
    }

    double result = integrate(qr, [](double x) {return 3 * x * x;});
    CHECK_CLOSE(2.0, result, 1e-12);

    auto qr_high_order = gauss(50);
    CHECK(qr_high_order.size() == 50);
    result = integrate(qr_high_order, [](double x) {return 101 * std::pow(x, 100);});
    CHECK_CLOSE(2.0, result, 1e-11);

    auto qr_odd = gauss(5);
    CHECK(qr_odd.size() == 5);
}
