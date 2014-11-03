#include "UnitTest++.h"
#include "ad2.h"
#include "util.h"
#include <iostream>

TEST(Create) {
    auto t = var<double,5>(-1.5);
    auto b = t + 5 + t - 3 - t + (-2);
    CHECK_ARRAY_CLOSE(t.v, b.v, 4, 1e-15);
}

TEST(Cube) {
    //TODO: Move over some tests.
    TIC
    auto t = var<double,18>(-1.5);
    auto t2 = t * t;
    auto t4 = t2 * t2;
    auto t8 = t4 * t4;
    auto res = t8 * t2 * t + t4 * t + 1;
    TOC("Mult");
    std::cout << res << std::endl;
}

TEST(Recip) {
    auto t = var<double,5>(1.0);
    auto t_r = recip(t);
    auto t2 = recip(t_r);
    CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
}

TEST(Log) {
    auto t = var<double,5>(1.0);
    auto log_t = log(t);
    auto t2 = exp(log_t);
    CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
}

TEST(Sqrt) {
    auto t = var<double,5>(1.0);
    auto sqrt_t = sqrt(t);
    auto t2 = sqrt_t * sqrt_t;
    CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
}

TEST(SinInv) {
    auto t = var<double,5>(0.5);
    auto asin_t = asin(t);
    auto t2 = sin(asin_t);
    CHECK_ARRAY_CLOSE(t.v,t2.v,5,1e-15);
}

int main() {
    return UnitTest::RunAllTests();
}
