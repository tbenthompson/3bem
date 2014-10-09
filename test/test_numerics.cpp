#include "UnitTest++.h"
#include "numerics.h"
#include <iostream>
#include "test_shared.h"

TEST(Naturals) {
    auto nats5 = naturals(3, 8);
    double correct[] = {3, 4, 5, 6, 7};
    CHECK_ARRAY_EQUAL(nats5, correct, 5);
}

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

TEST(SnFast) {
    int n = 50;
    auto x = random_list(n); 
    auto y = random_list(n); 
    for(int i = 0; i < n; i++) {
        for (int j = 1;j < 9; j++) {
            CHECK_CLOSE(s_n_fast(x[i], y[i], j),
                        s_n(x[i], y[i], j), 1e-12);
        }
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

    //Check that the weights are in the right range and sum properly
    double wt_sum = 0.0;
    for (auto xw: qr) {
        CHECK(std::fabs(xw.first) <= 1);
        wt_sum += xw.second;
        // std::cout << "x: " << xw.first << "   w: " << xw.second << std::endl;
    }
    CHECK_CLOSE(wt_sum, 2.0, 1e-12);

    double result = integrate(qr, [](double x) {return 3 * x * x;});
    CHECK_CLOSE(2.0, result, 1e-12);

    auto qr_high_order = gauss(50);
    CHECK(qr_high_order.size() == 50);
    result = integrate(qr_high_order, [](double x) {return 101 * std::pow(x, 100);});
    CHECK_CLOSE(2.0, result, 1e-11);

    auto qr_odd = gauss(5);
    CHECK(qr_odd.size() == 5);
}

TEST(DiligentiMapping) {
    auto quad = diligenti_mapping(18, -0.3, 7);
    double result = integrate(quad, 
        [](double x) { return log(std::fabs(x + 0.3));});
    CHECK_CLOSE(result, -1.908598917, 1e-6);

    quad = diligenti_mapping(18, 0.0, 7);
    result = integrate(quad, 
        [](double x) { return log(sqrt(pow(x, 4) + pow(x, 3) + pow(x,2)));});
    CHECK_CLOSE(result, -1.81569, 1e-5);

    result = integrate(quad, 
        [](double x) { return pow(x, 4) * log(sqrt(pow(x, 4) + pow(x, 3) + pow(x,2)));});
    CHECK_CLOSE(result, -0.000611395, 1e-6);
}
