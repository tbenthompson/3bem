#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <vector>
#include <cmath>
#include <functional>
#include "vec.h"

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

inline double linear_interp(double x_hat, double y_hat,
                            const std::array<double,3> corner_vals) {
    return (1 - x_hat - y_hat) * corner_vals[0] + 
           x_hat * corner_vals[1] + 
           y_hat * corner_vals[2];
}


inline std::array<double,3> 
tri_unscaled_normal(const std::array<std::array<double,3>,3> corners) {
    return cross(corners[2] - corners[0], corners[2] - corners[1]);
}

inline std::array<double,3> 
tri_normal(const std::array<std::array<double,3>,3> corners) {
    auto unscaled = tri_unscaled_normal(corners);
    return normalized(unscaled);
}

inline double tri_area(const std::array<double,3> unscaled_normal) {
    return 0.5 * hypot(unscaled_normal);
}
inline double tri_area(const std::array<std::array<double,3>,3> corners) {
    return tri_area(tri_unscaled_normal(corners));
}
 

std::vector<double> cheb_polys(double x_hat, int n_max);
double s_n(double x_hat, double y_hat, unsigned int n);
double s_n_fast(double x_hat, double y_hat, unsigned int n);
std::vector<double> cheb_pts_first_kind(unsigned int n_pts);
std::pair<double, double> legendre_and_n_minus_1(unsigned int n, double x);

#endif
