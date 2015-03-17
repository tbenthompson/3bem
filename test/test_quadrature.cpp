#include <iostream>
#include <iomanip>

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

TEST(AdaptiveQuadTensor) {
    auto result = adaptive_integrate<Vec2<Vec2<double>>>([] (double x) {
            return Vec2<Vec2<double>>{{
                {{x, 1}}, {{0, 2 * x}} 
            }};
        }, 0, 2, 1e-6);
    Vec2<Vec2<double>> correct{{
        {{2, 2}}, {{0, 4}}
    }};
    CHECK_ARRAY_CLOSE((double*)(&result[0]), (double*)(&correct[0]), 4, 1e-12);
}

double test_sinh_transform(bool iterated, double scale_factor, double a, double b,
    std::function<double(Vec<double,1>)> f) 
{
    double exact_error = 1e-15;
    auto n = 10;
    auto q = sinh_transform(n, a, b, iterated);
    decltype(q) q_scaled;
    for (const auto& pt: q) {
        q_scaled.push_back({pt.x_hat / scale_factor, pt.w / scale_factor});
    }

    auto res = integrate<double,1>(q_scaled, f);

    auto correct = adaptive_integrate<double>([&](double x) {
            return f({x});
        }, -1.0 / scale_factor, 1.0 / scale_factor, exact_error);

    double error = fabs(res - correct) / fabs(correct);
    return error;
}

TEST(SinhTransform) {
    double a = 0.0;
    double b = 0.001;
    auto fnc = [&](std::array<double,1> x) {
        return std::log(std::pow(x[0] - a, 2) + b * b);
    };
    CHECK(test_sinh_transform(false, 1, a, b, fnc) < 1e-4);
}

void ElliotJohnston2007Test(double a, bool iterated, std::vector<double> correct_error,
    std::function<std::function<double(Vec<double,1>)>(double)> f_builder) {
    size_t idx = 0;
    for (int log_b = -1; log_b >= -6; log_b--) {
        double b = std::pow(10, log_b);
        auto fnc = f_builder(b);
        double error = test_sinh_transform(iterated, 1, a, b, fnc);
        double error_error = fabs(error - correct_error[idx]) / fabs(correct_error[idx]);
        CHECK_CLOSE(error_error, 0.0, 1e-3);
        idx++;
    }
}

void ElliotJohnston2007I2(bool iterated, std::vector<double> correct_error) {
    double a = 0.25;
    ElliotJohnston2007Test(a, iterated, correct_error, 
        [=](double b) {
            return [=] (Vec<double,1> xs) {
                double x = xs[0];
                return (1 - x * x) / std::sqrt(std::pow(x - a, 2) + b * b);
            };
        }
    );
}

void ElliotJohnston2007I3(bool iterated, std::vector<double> correct_error) {
    double a = 0.25;
    ElliotJohnston2007Test(a, iterated, correct_error, 
        [=](double b) {
            return [=] (Vec<double,1> xs) {
                double x = xs[0];
                return (1 - x * x) / (std::pow(x - a, 2) + b * b);
            };
        }
    );
}

TEST(SinhTransformElliotJohnston2007I2NotIterated) {
    std::vector<double> correct_error = {
        1.8346e-11, 3.8039e-8, 1.8163e-6, 1.8068e-5, 8.0392e-5, 2.2492e-4
    };
    ElliotJohnston2007I2(false, correct_error);
}

TEST(SinhTransformElliotJohnston2007I2Iterated) {
    std::vector<double> correct_error = {
        2.2359e-6, 2.9066e-4, 2.2667e-3, 6.3398e-3, 1.0995e-2, 1.4802e-2
    };
    ElliotJohnston2007I2(true, correct_error);
}

TEST(SinhTransformElliotJohnston2007I3NotIterated) {
    std::vector<double> correct_error = {
        3.3753e-6, 4.6976e-3, 4.0194e-2, 1.2006e-1, 2.2969e-1, 3.4904e-1
    };
    ElliotJohnston2007I3(false, correct_error);
}

TEST(SinhTransformElliotJohnston2007I3Iterated) {
    std::vector<double> correct_error = {
        4.1270e-9, 8.6692e-7, 8.7829e-6, 3.0635e-5, 7.2012e-5, 1.4502e-4
    };
    ElliotJohnston2007I3(true, correct_error);
}

TEST(SinhTransformScaled) {
    double a = 0.0;
    double b = 0.001;
    auto fnc = [&](std::array<double,1> x) {
        return std::log(std::pow(x[0] - a, 2) + b * b);
    };
    CHECK(test_sinh_transform(false, 100, a, b, fnc) < 1e-4);
}

TEST(SinhSigmoidal2D) {
    size_t nt = 19;
    size_t nr = 9;
    std::vector<double> bs;
    for (int i = 1; i <= 6; i++) {
        bs.push_back(std::pow(10, -i) / std::sqrt(2));
    }
    std::vector<double> x0s = {0.0};
    std::vector<double> y0s = {0.0};
    std::vector<std::vector<double>> correct_errors{
        {1.02e-12, 8.97e-13, 1.14e-12, 1.33e-10, 3.04e-9, 3.44e-8}
    };
    double lambda = 1.5;

    for (size_t b_idx = 0; b_idx < bs.size(); b_idx++) {
        double b = bs[b_idx];
        for (size_t x0_idx = 0; x0_idx < x0s.size(); x0_idx++) {
            double x0 = x0s[x0_idx];
            double y0 = y0s[x0_idx];

            auto q = sinh_sigmoidal_transform(nt, nr, x0, y0, b, true);
            auto fnc = [&](std::array<double,2> x) {
                return 1.0 / std::pow(
                    std::pow(x[0] - x0, 2) + std::pow(x[1] - y0, 2) + b * b, lambda
                );
            };


            auto res = integrate<double,2>(q, [&] (Vec<double,2> x_hat) {
                    return fnc(x_hat);
                });

            double exact_error = 1e-9;
            auto correct = adaptive_integrate<double>([&](double x) {
                    return adaptive_integrate<double>([=](double y) {
                        return fnc({x, y});
                    }, 0.0, 1 - x, exact_error);
                }, 0.0, 1.0, exact_error);

            auto error = std::fabs(res - correct) / std::fabs(correct);
            std::cout << std::setprecision(16) << res << " " << correct << std::endl;
            double correct_error = correct_errors[x0_idx][b_idx];
            std::cout << error << " " << correct_error << std::endl;
            double error_error = fabs(error - correct_error) / fabs(correct_error);
            CHECK_CLOSE(error_error, 0.0, 1e-2);
        }
    }
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
