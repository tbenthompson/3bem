#ifndef TBEMNUMERICS_H
#define TBEMNUMERICS_H

#include <vector>
#include <cmath>
#include <functional>
#include "vec_ops.h"
#include "geometry.h"

namespace tbem {

template <size_t dim>
inline double inv_ref_facet_area();

template <>
inline double inv_ref_facet_area<2>() 
{
    return 0.5;
}

template <>
inline double inv_ref_facet_area<3>() 
{
    return 2.0;
}

template <size_t dim>
double facet_jacobian(const Vec<Vec<double,dim>,dim>& facet);
template <>
inline double facet_jacobian<2>(const Vec<Vec<double,2>,2>& facet) {
    return inv_ref_facet_area<2>() * hypot(unscaled_normal(facet));
}
template <>
inline double facet_jacobian<3>(const Vec<Vec<double,3>,3>& facet) {
    return inv_ref_facet_area<3>() * tri_area(facet);
}

// Map from [-1, 1] to [0, 1].
inline double from_11_to_01(double x) 
{
    return 0.5 * x + 0.5;
}

// Map from [0, 1] to [-1, 1].
inline double from_01_to_11(double x) 
{
    return 2 * x - 1;
}

// Map from [a, b] to [-1, 1]
inline double real_to_ref(double x, double a, double b) 
{
    return 2.0 * ((x - a) / (b - a)) - 1.0;
}

// Map from [-1, 1] to [a, b]
inline double ref_to_real(double x_hat, double a, double b) 
{
    return a + (b - a) * ((x_hat + 1.0) / 2.0);
}

// 2D linear basis on [-1,1]
inline Vec2<double> linear_basis(const Vec<double,1>& x_hat) 
{
    return {{0.5 - 0.5 * x_hat[0], 0.5 + 0.5 * x_hat[0]}};
}

// 3D linear basis on (0,0)-(1,0)-(0,1)
inline Vec3<double> linear_basis(const Vec<double,2>& x_hat) 
{
    return {{1 - x_hat[0] - x_hat[1], x_hat[0], x_hat[1]}};
}

template <size_t dim>
inline double linear_interp(const Vec<double,dim-1>& x_hat,
    const Vec<double,dim>& corner_vals) 
{
    return dot_product(linear_basis(x_hat), corner_vals);
}

inline Vec2<double> ref_to_real(const Vec<double,1>& x_hat,
    const std::array<Vec2<double>,2>& locs) 
{
    auto basis = linear_basis(x_hat);
    return {{
        dot_product(basis, Vec2<double>{{locs[0][0], locs[1][0]}}),
        dot_product(basis, Vec2<double>{{locs[0][1], locs[1][1]}}),
    }};
}

inline Vec3<double> ref_to_real(const Vec<double,2>& x_hat,
    const Vec<Vec<double,3>,3>& locs) 
{
    auto basis = linear_basis(x_hat);
    return {{
        dot_product(basis, Vec3<double>{{locs[0][0], locs[1][0], locs[2][0]}}),
        dot_product(basis, Vec3<double>{{locs[0][1], locs[1][1], locs[2][1]}}),
        dot_product(basis, Vec3<double>{{locs[0][2], locs[1][2], locs[2][2]}})
    }};
}
 
} //END NAMESPACE tbem

#endif
