#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "numerics.h"
#include "numbers.h"
#include <iostream>
#include "util.h"

using namespace tbem;


TEST(LinearMapping) {
    auto gen100 = ac::fix(100, ac::generator<double>());
    auto arb = ac::make_arbitrary(gen100, gen100, gen100);
    ac::check<double, double, double>(
        [](double real, double v0, double v1) {
            double ref = real_to_ref(real, v0, v1);
            std::array<Vec2<double>,2> locs = {{{v0,0}, {v1,0}}};
            double real2 = ref_to_real({ref}, locs)[0];
            return std::fabs(real - real2) < 1e-11;
        }, 100, arb);
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

TEST(RefToRealGradient) {
    Vec<Vec<double,3>,3> f{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    auto result = ref_to_real_gradient({0, 0}, f);
    Vec<Vec<double,2>,3> correct{{{1, 0}, {0, 1}, {0, 0}}};
    CHECK_EQUAL(result, correct);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
