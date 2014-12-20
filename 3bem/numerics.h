#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <vector>
#include <cmath>
#include <functional>
#include "vec.h"

namespace tbem {

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

// 2D linear basis on [-1,1]
inline Vec2<double> linear_basis(const std::array<double,1>& x_hat) {
    return {0.5 - 0.5 * x_hat[0], 0.5 + 0.5 * x_hat[0]};
}

// 3D linear basis on (0,0)-(1,0)-(0,1)
inline Vec3<double> linear_basis(const std::array<double,2>& x_hat) {
    return {1 - x_hat[0] - x_hat[1], x_hat[0], x_hat[1]};
}

template <int dim>
inline double linear_interp(const std::array<double,dim-1>& x_hat,
                            const Vec<double,dim>& corner_vals) {
    return dot(linear_basis(x_hat), corner_vals);
}

inline Vec2<double> ref_to_real(const std::array<double,1>& x_hat,
                                const std::array<Vec2<double>,2>& locs) {
    auto basis = linear_basis(x_hat);
    return {
        dot(basis, {locs[0][0], locs[1][0]}),
        dot(basis, {locs[0][1], locs[1][1]}),
    };
}

inline Vec3<double> ref_to_real(const std::array<double,2>& x_hat,
                                const std::array<Vec3<double>,3>& locs) {
    auto basis = linear_basis(x_hat);
    return {
        dot(basis, {locs[0][0], locs[1][0], locs[2][0]}),
        dot(basis, {locs[0][1], locs[1][1], locs[2][1]}),
        dot(basis, {locs[0][2], locs[1][2], locs[2][2]})
    };
}
 

std::vector<double> cheb_polys(double x_hat, int n_max);
double s_n(double x_hat, double y_hat, unsigned int n);
double s_n_fast(double x_hat, double y_hat, unsigned int n);
std::vector<double> cheb_pts_first_kind(unsigned int n_pts);
std::pair<double, double> legendre_and_n_minus_1(unsigned int n, double x);

} //END NAMESPACE tbem

#endif
