#include "ad.h"
#include "UnitTest++.h"
#include "util.h"

using namespace ad;

TEST(CreateVar) {
    auto x = D_d<5>::D_from_id(3.0);
    CHECK_EQUAL(x.v, 3.0);
    CHECK_EQUAL(x.d.v, 1.0);
}

TEST(CreateConstant) {
    auto x = D_d<1>::D_from_const(2.0);
    CHECK_EQUAL(x.v, 2.0);
    CHECK_EQUAL(x.d.v, 0.0);
}

TEST(AddSub) {
    auto t = D_d<2>::D_from_id(3.0);
    auto t2 = t + (t + 3) - 2 - t + -1;
    CHECK_EQUAL(t2.v, 3.0);
    CHECK_EQUAL(t2.d.v, 1.0);
}

int factorial(int n) {
    return (n == 0) ? 1 : (n * factorial(n - 1));
}

TEST(PolyDerivative) {
    auto t = D_d<11>::D_from_id(-1.5);
    auto t2 = t * t;
    auto t4 = t2 * t2;
    auto t8 = t4 * t4;
    auto res = t8 * t2 * t + t4 * t + 1;

    //taylor series
    double exact[12] = {
        -93.0913, 659.628, -2148.13, 4251.27, -5645.86, 
        5263.47, -3508.31, 1670.63, -556.875, 123.75, -16.5, 1
    };
    std::array<double,12> vals;
    res.get_all(vals);
    for (int i = 0; i < 12;i++) {
        double deriv_exact = exact[i] * factorial(i);
        double error = std::fabs((vals[i] - deriv_exact) / deriv_exact);
        CHECK_CLOSE(error, 0, 1e-5);
    }
}

TEST(ExpDerivative) {
    auto t = D_d<3>::D_from_id(-1.5);
    auto exp_t = exp(t);
    std::cout << exp_t << std::endl;
}

TEST(RecursiveDataStructure) {
    auto x2 = D_d<5>::D_from_id(3.0);
    std::cout << x2 << std::endl;
}

int main() {
    UnitTest::RunAllTests();
}
