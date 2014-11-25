#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "numerics.h"
#include <iostream>
#include "util.h"

TEST(Naturals) {
    auto nats5 = naturals(3, 8);
    double correct[] = {3, 4, 5, 6, 7};
    CHECK_ARRAY_EQUAL(nats5, correct, 5);
}

TEST(LinearMapping) {
    auto gen100 = ac::fix(100, ac::generator<double>());
    auto arb = ac::make_arbitrary(gen100, gen100, gen100);
    ac::check<double, double, double>(
        [](double real, double v0, double v1) {
            double ref = real_to_ref(real, v0, v1);
            std::array<Vec2<double>,2> locs = {{{v0,0}, {v1,0}}};
            double real2 = ref_to_real({ref}, locs)[0];
            return std::fabs(real - real2) < 1e-12;
        }, 100, arb);
}

TEST(ChebPolys) {
    auto gen6pi = ac::fix(20, ac::generator<double>());
    auto genint10 = ac::fix(10, ac::generator<unsigned int>());
    auto arb = ac::make_arbitrary(gen6pi, genint10);
    ac::check<double, unsigned int>(
        [](double theta, unsigned int n_polys) {
            auto res = cheb_polys(theta, n_polys);
            bool prop_true = true;
            for (unsigned int j = 0; j < n_polys + 1; j++) {
                prop_true |= std::fabs(cos(theta * j) - res[j]) < 1e-12;
            }
            return prop_true; 
        }, 100, arb);
}

TEST(SnFast) {
    auto gen11 = ac::fix(1, ac::generator<double>());
    auto genint10 = ac::fix(10, ac::generator<unsigned int>());
    auto arb = ac::make_arbitrary(gen11, gen11, genint10);
    ac::check<double, double, unsigned int>(
        [](double x, double y, unsigned int pts){
            double val1 = s_n(x, y, pts + 1);
            double val2 = s_n_fast(x, y, pts + 1);
            return std::fabs(val1 - val2) < 1e-12;
        }, 100, arb);
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

TEST(AreaTri) {
    double result = tri_area({{{0,0,0},{1,0,0},{0,1,0}}});
    CHECK_CLOSE(result, 0.5, 1e-12);
    result = tri_area({{{1,1,1},{3,1,1},{3,3,1}}});
    CHECK_CLOSE(result, 2.0, 1e-12);
}


TEST(LinearInterp) {
    CHECK_CLOSE(linear_interp<3>({0,0},{1,0,0}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp<3>({1,0},{0,1,0}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp<3>({0,1},{0,0,1}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp<3>({0.5,0.5},{0,0,1}), 0.5, 1e-12);
    CHECK_CLOSE(linear_interp<3>({0.0,0.5},{0,0,1}), 0.5, 1e-12);
}

TEST(LinearInterpOnes) {
    auto gen1 = ac::fix(1, ac::generator<double>());
    auto arb = ac::make_arbitrary(gen1, gen1);
    ac::check<double, double>(
        [](double x_hat, double y_hat) {
            double result = linear_interp<3>({x_hat, y_hat}, {1,1,1});
            double exact = 1.0;
            return std::fabs(result - exact) < 1e-12;
        }, 30, arb);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
