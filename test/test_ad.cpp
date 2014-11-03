#include "taylor.h"
#include "UnitTest++.h"

using namespace taylor_ad;

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

    std::cout << res << std::endl;
}

TEST(RecursiveDataStructure) {
    auto x2 = D_d<5>::D_from_id(3.0);
    std::cout << x2 << std::endl;
}

int main() {
    UnitTest::RunAllTests();
}
