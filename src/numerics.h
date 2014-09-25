#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <vector>
#include <cmath>
#include <functional>

inline double real_to_ref(double x, double a, double b) {
    return 2.0 * ((x - a) / (b - a)) - 1.0;
}

inline double ref_to_real(double x_hat, double a, double b) {
    return a + (b - a) * ((x_hat + 1.0) / 2.0);
}

inline double hypot(double x, double y) {
    return sqrt(x * x + y * y);
}
inline double hypot(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

typedef std::vector<std::pair<double, double>> QuadratureRule;

std::vector<int> naturals(int min, int max);
std::vector<int> naturals(int max);

std::vector<double> cheb_polys(double x_hat, int n_max);
double s_n(double x_hat, double y_hat, int n);
std::vector<double> cheb_pts_first_kind(unsigned int n_pts);
std::pair<double, double> legendre_and_n_minus_1(unsigned int n, double x);
QuadratureRule double_exp(int n, double h);
QuadratureRule gauss(unsigned int n);
QuadratureRule diligenti_mapping(unsigned int n, double x0, int q);
double integrate(QuadratureRule& qr, std::function<double (double)> fnc);
#endif
