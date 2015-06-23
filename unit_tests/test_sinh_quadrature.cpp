#include "catch.hpp"
#include "sinh_quadrature.h"
#include "quadrature.h"
#include "vec_ops.h"
#include "adaptive_quad.h"
#include <functional>

using namespace tbem;

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, size_t dim>
T integrate(const QuadRule<dim>& qr, 
    const std::function<T(std::array<double,dim>)>& fnc) 
{
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (size_t i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

double test_sinh_transform(bool iterated, double scale_factor, double a, double b,
    std::function<double(Vec<double,1>)> f) 
{
    double exact_error = 1e-15;
    auto n = 10;
    auto q = sinh_transform(gauss(n), a, b, iterated);
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

TEST_CASE("SinhTransform", "[quadrature]") {
    double a = 0.0;
    double b = 0.001;
    auto fnc = [&](std::array<double,1> x) {
        return std::log(std::pow(x[0] - a, 2) + b * b);
    };
    REQUIRE(test_sinh_transform(false, 1, a, b, fnc) < 1e-4);
}

TEST_CASE("SinhTransformScaled", "[quadrature]") {
    double a = 0.0;
    double b = 0.001;
    auto fnc = [&](std::array<double,1> x) {
        return std::log(std::pow(x[0] - a, 2) + b * b);
    };
    REQUIRE(test_sinh_transform(false, 100, a, b, fnc) < 1e-4);
}

void test_sinh_sigmoidal(double lambda, size_t nt, size_t nr, double x0, double y0) {
    std::vector<double> bs;
    for (int i = 1; i <= 6; i++) {
        bs.push_back(std::pow(10, -i));
    }

    size_t sinh_eval = 0;
    size_t adapt_eval = 0;
    for (size_t b_idx = 0; b_idx < bs.size(); b_idx++) {
        double b = bs[b_idx];
        auto q = sinh_sigmoidal_transform(gauss(nt), gauss(nr), x0, y0, b, false);

        auto res = integrate<double,2>(q, [&] (Vec<double,2> x_hat) {
                sinh_eval++;
                double dx = x_hat[0] - x0;
                double dy = x_hat[1] - y0;
                return 1.0 / std::pow(dx * dx + dy * dy + b * b, lambda);
            });

        double exact_error = 1e-4;
        auto correct = adaptive_integrate<double>([&](double x) {
                return adaptive_integrate<double>([&](double y) {
                    adapt_eval++;
                    double dx = x - x0;
                    double dy = y - y0;
                    return std::pow(dx * dx + dy * dy + b * b, -lambda);
                }, 0.0, 1 - x, exact_error);
            }, 0.0, 1.0, exact_error);

        auto error = std::fabs(res - correct) / std::fabs(correct);
        REQUIRE_CLOSE(error, 0.0, 1e-2);
    }
}

TEST_CASE("SinhSigmoidal2D", "[quadrature]") {
    test_sinh_sigmoidal(1.5, 17, 11, 0.0, 0.0);
    test_sinh_sigmoidal(0.5, 11, 5, 0.2, 0.4);
    test_sinh_sigmoidal(1.0, 20, 10, 0.1, 0.1);
}
