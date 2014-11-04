#include "UnitTest++.h"
#include "taylor.h"
#include "util.h"
#include <iostream>

TEST(CreateConst) {
    auto t = Td<4>::constant(1.0);
    double exact[5] = {1.0, 0,0,0,0};
    CHECK_ARRAY_CLOSE(t.c, exact, 5, 1e-15);
}

TEST(CreateVar) {
    auto t = Td<4>::var(-2.0);
    double exact[5] = {-2.0, 1.0, 0,0,0};
    CHECK_ARRAY_CLOSE(t.c, exact, 5, 1e-15);
}

TEST(AddSub) {
    auto t = Td<4>::var(-4.0);
    auto res = 5 + t - 3 + t - 2 - t + t - t;
    double exact[5] = {-4.0, 1.0, 0,0,0};
    CHECK_ARRAY_CLOSE(res.c, exact, 5, 1e-15);
}

TEST(AddSub2) {
    auto t = Td<4>::var(-4.0);
    auto res = t + t;
    double exact[5] = {-8.0, 2.0, 0,0,0};
    CHECK_ARRAY_CLOSE(res.c, exact, 5, 1e-15);
}

TEST(AddSub3) {
    auto t = Td<4>::var(-4.0);
    auto res = 3 - t;
    double exact[5] = {7.0, -1.0, 0,0,0};
    CHECK_ARRAY_CLOSE(res.c, exact, 5, 1e-15);
}

TEST(Mul) {
    auto t = Td<7>::var(-1.5);
    auto t2 = t * t; 
    CHECK_EQUAL(t2.n_coeffs, 3);
    auto t4 = t2 * t2; 
    CHECK_EQUAL(t4.n_coeffs, 5);
    auto res = t4 * t2 * t + t2 * t + 8;
    double exact[8] = {
        -1595 / 128.0, 5535/64.0, -5247/32.0, 2851/16.0,
        -945/8.0, 189/4.0, -21/2.0, 1.0
    };
    CHECK_EQUAL(res.n_coeffs, 8);
    CHECK_ARRAY_CLOSE(res.c, exact, 8, 1e-15);
}

TEST(Div) {
    auto t = Td<7>::var(-1.5);
    auto t2 = t * t;
    auto divided = t2 / t;
    CHECK_EQUAL(divided.n_coeffs, 8);
    CHECK_ARRAY_CLOSE(t.c, divided.c, 8, 1e-15);
}

TEST(MulDivScalars) {
    auto t = Td<3>::var(3.0);
    auto res = 3 / t;
    double exact[4] = {1, -1 / 3.0, 1 / 9.0, -1.0 / 27.0};
    CHECK_ARRAY_CLOSE(res.c, exact, 4, 1e-15);
}

TEST(Sqrt) {
    auto t = Td<3>::var(3.0);
    auto res = sqrt(t);
    auto rt3 = std::sqrt(3);
    double exact[4] = {rt3, (1 / (2 * rt3)), -1 / (24 * rt3), 1 / (144 * rt3)};
    CHECK_ARRAY_CLOSE(res.c, exact, 4, 1e-15);
}

TEST(Log) {
    auto t = Td<3>::var(2);
    auto res = log(t);
    double exact[4] = {std::log(2), 1.0 / 2, -1.0 / 8, 1.0 / 24.0};
    CHECK_ARRAY_CLOSE(res.c, exact, 4, 1e-15);
}

TEST(Exp) {
    auto t = Td<3>::var(2);
    auto res = exp(log(t));
    double exact[4] = {2.0, 1.0, 0,0};
    CHECK_ARRAY_CLOSE(res.c, exact, 4, 1e-15);
    res = exp(t);
    auto exp2 = std::exp(2);
    double exact2[4] = {exp2, exp2, exp2 / 2, exp2 / 6};
    CHECK_ARRAY_CLOSE(res.c, exact2, 4, 1e-15);
}

TEST(Pow) {
    auto t = Td<3>::var(2);
    auto res = pow(t, 7/5.);
    auto A = std::pow(2, 2/5.);
    double exact[4] = {A * 2, A * (7 / 5.), A * (7 / 50.), -A * (7 / 500.)};
    CHECK_ARRAY_CLOSE(res.c, exact, 4, 1e-15);
}

TEST(Eval) {
    auto t = Td<3>::var(1);
    auto t2 = pow(t, 2);
    std::cout << t2 << std::endl;
    std::cout << t2.eval(3) << std::endl; 
}

TEST(Profiling) {
    const int d = 5;
    const int n = (int)1e6;

    TIC
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        auto t = Td<d>::var(0.8);
        auto t2 = t * t;
        auto t4 = t2 * t2;
        auto invt4 = pow(1 / sqrt(t4), 1.034299);
        t = invt4 + 10 + t2;
    }
    TOC("Lots of mults")
    // std::cout << t << std::endl;
}

int main() {
    UnitTest::RunAllTests();
}
