#include "util.h"
#include "quadrature.h"
#include "numerics.h"
#include "vec_ops.h"

using namespace tbem;

double single_tri(const std::function<double(Vec<double,2>)>& f,
    const QuadRule<1>& gauss_theta, const QuadRule<1>& gauss_r,
    double b, const Vec<Vec<double,3>,3>& tri) 
{
    double tri_jacobian = 2.0 * tri_area(tri);
    double out = 0.0;
    std::vector<double> r_unscaleds(gauss_r.size());
    std::vector<double> r_jacobians(gauss_r.size());
    double mu_0 = 0.5 * std::asinh(1.0 / b);
    //THIS SECTION: 32ms
    for (size_t j = 0; j < gauss_r.size(); j++) {
        double s = gauss_r[j].x_hat[0];
        double arg = mu_0 * (s + 1);
        double expneg = std::exp(-arg);
        double exppos = 1.0 / expneg;
        double sinh_val = 0.5 * (exppos - expneg);
        double cosh_val = 0.5 * (exppos + expneg);
        double r = b * sinh_val;
        r_unscaleds[j] = r;
        r_jacobians[j] = r * b * mu_0 * cosh_val * gauss_r[j].w;
    }

    //THIS SECTION: 120ms
    const double quarterpi = M_PI / 4.0;
    const double sinquarterpi = std::sin(quarterpi);
    for (size_t i = 0; i < gauss_theta.size(); i++) {
        double sigma = (gauss_theta[i].x_hat[0] + 1.0) / 2.0;
        double sig_transform = std::pow(sigma, 2) /
            (std::pow(sigma, 2) + std::pow(1 - sigma, 2));
        double theta = sig_transform * (M_PI / 2.0);
        double theta_jacobian = (M_PI / 2.0) * sigma * (1 - sigma) /
            std::pow(std::pow(sigma, 2) + std::pow(1 - sigma, 2), 2);
        double cos_theta = std::cos(theta);
        double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

        double R_theta = sinquarterpi / std::sin(theta + quarterpi);
        double outer_jacobian = R_theta * tri_jacobian * theta_jacobian * gauss_theta[i].w;
        double inner_sum = 0.0;
        for (size_t j = 0; j < gauss_r.size(); j++) {
            double r = R_theta * r_unscaleds[j];
            double x_h = r * cos_theta;
            double y_h = r * sin_theta;
            double x = (1.0 - x_h - y_h) * tri[0][0] + x_h * tri[1][0] + y_h * tri[2][0];
            double y = (1.0 - x_h - y_h) * tri[0][1] + x_h * tri[1][1] + y_h * tri[2][1];
            double w = r_jacobians[j];
            inner_sum += (x * y + 1 / std::pow(y, 3)) * w;
        }
        out += outer_jacobian * inner_sum;
    }

    return out;
};

double sinh_sigmoidal_integrate(const std::function<double(Vec<double,2>)>& f,
    const QuadRule<1>& gauss_theta, const QuadRule<1>& gauss_r, double x0,
    double y0, double b) 
{
    Vec<double,3> pt0{0.0, 0.0, 0.0};
    Vec<double,3> pt1{1.0, 0.0, 0.0};
    Vec<double,3> pt2{0.0, 1.0, 0.0};
    Vec<double,3> singular_pt{x0, y0, 0};

    double result = single_tri(f, gauss_theta, gauss_r, b, {singular_pt, pt0, pt2});
    result += single_tri(f, gauss_theta, gauss_r, b, {singular_pt, pt2, pt1});
    result += single_tri(f, gauss_theta, gauss_r, b, {singular_pt, pt0, pt1});

    return result;
}

int main() {
    size_t N = 180000;

    auto rules = gauss_set(10, 40);
    TIC
    double val = 0.0;
#pragma omp parallel for
    for(size_t i = 0; i < N; i++) {
        size_t n_q = i % 30 + 10;
        // auto q = sinh_sigmoidal_transform(rules.at(n_q), rules.at(n_q),
        //     0.5, 0.1, 0.1, false);
        val += sinh_sigmoidal_integrate([](Vec<double,2> x){return 1.0;},
            rules.at(n_q), rules.at(n_q), 0.5, 0.1, 0.1);
    }
    std::cout << val << std::endl;
    TOC("Constructing " + std::to_string(N) + " sinh-sigmoidal quadrature rules");
}
