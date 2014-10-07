#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <vector>
#include <cmath>
#include <functional>

// Map from [-1, 1] to [0, 1].
inline double from_11_to_01(double x) {
    return 0.5 * x + 0.5;
}

// Map from [0, 1] to [-1, 1].
inline double from_01_to_11(double x) {
    return 2 * x - 1;
}

// Map from [a, b] to [-1, 1]
inline double real_to_ref(double x, double a, double b) {
    return 2.0 * ((x - a) / (b - a)) - 1.0;
}

// Map from [-1, 1] to [a, b]
inline double ref_to_real(double x_hat, double a, double b) {
    return a + (b - a) * ((x_hat + 1.0) / 2.0);
}

inline double hypot2(double x, double y) {return x * x + y * y;}
inline double hypot2(double x, double y, double z) {return x * x + y * y + z * z;}
inline double hypot2(std::array<double,3> v0) {return hypot2(v0[0], v0[1], v0[2]);}
inline double hypot(double x, double y) {return sqrt(hypot2(x, y));}
inline double hypot(double x, double y, double z) {return sqrt(hypot2(x, y, z));}

template <int dim>
inline double dist2(const std::array<double,dim> v0, 
                    const std::array<double,dim> v1) {
    double result = 0.0;
    for (int d = 0; d < dim; d++) {
        const double delta = (v0[d] - v1[d]);
        result += delta * delta;
    }
    return result;
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
