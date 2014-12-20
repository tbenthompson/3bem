#include <iostream>

#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "util.h"
#include "numerics.h"
#include "quadrature.h"
#include "adaptive_quad.h"

using namespace tbem;

TEST(GaussQuadrature) {
    auto qr = gauss(4);
    CHECK(qr.size() == 4);

    //Check that the weights are in the right range and sum properly
    double wt_sum = 0.0;
    for (auto xw: qr) {
        CHECK(std::fabs(xw.x_hat[0]) <= 1);
        wt_sum += xw.w;
        // std::cout << "x: " << xw.first << "   w: " << xw.second << std::endl;
    }
    CHECK_CLOSE(wt_sum, 2.0, 1e-12);

    double result = integrate<double,1>(qr, [](std::array<double,1> x) {
            return 3 * x[0] * x[0];
        });
    CHECK_CLOSE(2.0, result, 1e-12);

    auto qr_high_order = gauss(50);
    CHECK(qr_high_order.size() == 50);
    result = integrate<double,1>(qr_high_order, [](std::array<double,1> x) {
            return 101 * std::pow(x[0], 100);
        });
    CHECK_CLOSE(2.0, result, 1e-11);

    auto qr_odd = gauss(5);
    CHECK(qr_odd.size() == 5);
}

TEST(GaussExactness) {
    auto genint = ac::fix(250, ac::generator<unsigned int>());
    auto arb = ac::make_arbitrary(genint);
    ac::check<unsigned int>(
        [](unsigned int n) {
            if (n == 0) {
                return true;
            }
            int g = 2 * n - 1;
            auto q = gauss(n);
            double result = integrate<double,1>(q, [&](std::array<double,1> x) {
                    return (g + 1) * pow(x[0], g);
                });
            double exact = 2.0 * ((g + 1) % 2);
            return std::fabs(exact - result) < 1e-13;
        }, 30, arb);
}

TEST(QuadRule2dConstructor) {
    QuadRule<2> q(10);
    CHECK_EQUAL(q.size(), 10);
}

TEST(TensorProduct) {
    auto g1d = gauss(2);
    auto g2d = tensor_product(g1d, g1d);
    double x_hat[4] = {-0.57735, -0.57735, 0.57735, 0.57735};
    double y_hat[4] = {-0.57735, 0.57735, -0.57735, 0.57735};
    double weights[4] = {1,1,1,1};
    for (unsigned int i = 0; i < g2d.size(); i++) {
        CHECK_CLOSE(g2d[i].x_hat[0], x_hat[i], 1e-4);
        CHECK_CLOSE(g2d[i].x_hat[1], y_hat[i], 1e-4);
        CHECK_CLOSE(g2d[i].w, weights[i], 1e-4);
    }
}

TEST(TensorProductIntegrate) {
    auto g2d = tensor_gauss(2);
    double result = integrate<double,2>(g2d, [](std::array<double,2> x) {
        return pow(x[0] - 0.5,3) + pow(x[1], 2);
    });
    CHECK_CLOSE(result, -(7.0/6), 1e-6);
}

TEST(TensorProductIntegrate2) {
    auto g2d = tensor_gauss(35);
    double result = integrate<double,2>(g2d, [](std::array<double,2> x) {
        return std::exp(x[0] / (x[1] + 1.1));
    });
    CHECK_CLOSE(result, 38.6995, 1e-4);
}

TEST(IntegrateAreaOfUnitTriangle) {
    auto q2d = tensor_gauss(2);
    auto q2d_tri = square_to_tri(q2d);
    double result = integrate<double,2>(q2d_tri, [](std::array<double,2> x) {return 1.0;});
    CHECK_CLOSE(result, 0.5, 1e-10);
}

TEST(IntegrateTriPoly) {
    auto g2d = tri_gauss(13);
    double result = integrate<double,2>(g2d, [](std::array<double,2> x) {
        return pow(x[0],23) + pow(x[1], 19);
    });
    CHECK_CLOSE(result, 17. / 4200, 1e-15);
}

void test_tri_integrate(QuadRule<2> q2d_tri) {
    double result = integrate<double,2>(q2d_tri, [](std::array<double,2> x) {
        return std::exp(x[0] / (x[1] - 1.1));
    });
    CHECK_CLOSE(result, 0.337429, 1e-6);
    result = integrate<double,2>(q2d_tri, [](std::array<double,2> x) {
            return std::exp(x[0] / (x[1] + 1.1));
        });
    CHECK_CLOSE(result, 0.656602, 1e-6);
}

TEST(TriangleIntegrate) {
    auto q2d_tri = tri_gauss(8);
    test_tri_integrate(q2d_tri);
}

TEST(AdaptiveQuad) {
    double result = adaptive_integrate<double>([] (double x) {
            return std::exp(x);
        }, 0, 2, 1e-3);
    CHECK_CLOSE(result, std::exp(2) - 1, 1e-6);
}

TEST(AdaptiveQuadCloseToSingular) {
    double a = 0.5;
    double result = adaptive_integrate<double>([a] (double x) {
            return 1.0 / std::sqrt(x*x + a);
        }, -0.5, 0.5, 1e-3);
    double exact = -std::log(-1 + std::sqrt(1 + 4 * a)) +
                    std::log(1 + std::sqrt(1 + 4 * a));
    CHECK_CLOSE(result, exact, 1e-6);
}

TEST(AdaptiveQuadVec) {
    Vec3<double> result = adaptive_integrate<Vec3<double>>([] (double x) {
            return Vec3<double>{
                std::exp(x), std::exp(2*x), std::exp(3*x)
            };
        }, 0, 2, 1e-6);
    Vec3<double> exact = {
        std::exp(2) - 1, (std::exp(4) - 1) / 2.0, (std::exp(6) - 1) / 3.0
    };
    CHECK_ARRAY_CLOSE(result, exact, 3, 1e-4);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
