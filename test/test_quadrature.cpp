#include <iostream>

#include "UnitTest++.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

#include "util.h"
#include "numerics.h"
#include "quadrature.h"

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

TEST(GaussExactness) {
    auto genint = ac::fix(500, ac::generator<unsigned int>());
    auto arb = ac::make_arbitrary(genint);
    ac::check<unsigned int>(
        [](unsigned int n) {
            int g = 2 * n - 1;
            auto q = gauss(n);
            double result = integrate(q, [=](double x) {return (g + 1) * pow(x, g);});
            double exact = 2.0 * ((g + 1) % 2);
            return std::fabs(exact - result) < 1e-13;
        }, 30, arb);
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

TEST(QuadratureRule2DConstructor) {
    QuadratureRule2D q(10);
    CHECK_EQUAL(q.x_hat.size(), 10);
    CHECK_EQUAL(q.y_hat.size(), 10);
    CHECK_EQUAL(q.weights.size(), 10);
}

TEST(TensorProduct) {
    auto g1d = gauss(2);
    auto g2d = tensor_product(g1d, g1d);
    double x_hat[4] = {-0.57735, -0.57735, 0.57735, 0.57735};
    double y_hat[4] = {-0.57735, 0.57735, -0.57735, 0.57735};
    double weights[4] = {1,1,1,1};
    CHECK_ARRAY_CLOSE(g2d.x_hat, x_hat, 4, 1e-4);
    CHECK_ARRAY_CLOSE(g2d.y_hat, y_hat, 4, 1e-4);
    CHECK_ARRAY_CLOSE(g2d.weights, weights, 4, 1e-4);
}

TEST(TensorProductIntegrate) {
    auto g2d = tensor_gauss(2);
    double result = integrate(g2d, [](double x,double y) {
        return pow(x - 0.5,3) + pow(y, 2);
    });
    CHECK_CLOSE(result, -(7.0/6), 1e-6);
}

TEST(TensorProductIntegrate2) {
    auto g2d = tensor_gauss(35);
    double result = integrate(g2d, [](double x,double y) {
        return std::exp(x / (y + 1.1));
    });
    CHECK_CLOSE(result, 38.6995, 1e-4);
}

TEST(IntegrateAreaOfUnitTriangle) {
    auto q2d = tensor_gauss(2);
    auto q2d_tri = square_to_tri(q2d);
    double result = integrate(q2d_tri, [](double x, double y) {return 1.0;});
    CHECK_CLOSE(result, 0.5, 1e-10);
}

TEST(IntegrateTriPoly) {
    auto g2d = tri_gauss(13);
    double result = integrate(g2d, [](double x,double y) {
        return pow(x,23) + pow(y, 19);
    });
    CHECK_CLOSE(result, 17. / 4200, 1e-15);
}

void test_tri_integrate(QuadratureRule2D q2d_tri) {
    double result = integrate(q2d_tri, [](double x,double y) {
        return std::exp(x / (y - 1.1));
    });
    CHECK_CLOSE(result, 0.337429, 1e-6);
    result = integrate(q2d_tri, [](double x,double y) {
        return std::exp(x / (y + 1.1));
    });
    CHECK_CLOSE(result, 0.656602, 1e-6);
}

TEST(TriangleIntegrate) {
    auto q2d_tri = tri_gauss(8);
    test_tri_integrate(q2d_tri);
}

TEST(TriangleIntegrateDoubleExp) {
    auto q2d_tri = tri_double_exp(7, 0.3);
    test_tri_integrate(q2d_tri);
}

TEST(AreaTri) {
    double result = tri_area({{{0,0,0},{1,0,0},{0,1,0}}});
    CHECK_CLOSE(result, 0.5, 1e-12);
    result = tri_area({{{1,1,1},{3,1,1},{3,3,1}}});
    CHECK_CLOSE(result, 2.0, 1e-12);
}

TEST(LinearInterp) {
    CHECK_CLOSE(linear_interp(0,0,{1,0,0}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp(1,0,{0,1,0}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp(0,1,{0,0,1}), 1.0, 1e-12);
    CHECK_CLOSE(linear_interp(0.5,0.5,{0,0,1}), 0.5, 1e-12);
    CHECK_CLOSE(linear_interp(0.0,0.5,{0,0,1}), 0.5, 1e-12);
}

TEST(LinearInterpOnes) {
    auto gen1 = ac::fix(1, ac::generator<double>());
    auto arb = ac::make_arbitrary(gen1, gen1);
    ac::check<double, double>(
        [](double x_hat, double y_hat) {
            double result = linear_interp(x_hat, y_hat, {1,1,1});
            double exact = 1.0;
            return std::fabs(result - exact) < 1e-12;
        }, 30, arb);
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
