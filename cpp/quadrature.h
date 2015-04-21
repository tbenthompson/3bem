#ifndef __QEWKJQWEMNZ_QUADRATURE_H
#define __QEWKJQWEMNZ_QUADRATURE_H
#include <vector>
#include <utility>
#include <functional>
#include <map>

namespace tbem {

template <size_t dim>
struct QuadPt {
    const std::array<double,dim> x_hat;
    const double w;
};

template <size_t dim>
using QuadRule = std::vector<QuadPt<dim>>;

QuadRule<1> gauss(size_t n);
QuadRule<2> tensor_gauss(size_t n_pts);
QuadRule<2> tri_gauss(size_t n_pts);
template <size_t dim>
QuadRule<dim-1> gauss_facet(size_t n_q);
template <>
inline QuadRule<2> gauss_facet<3>(size_t n_q) {
    return tri_gauss(n_q);
}
template <>
inline QuadRule<1> gauss_facet<2>(size_t n_q) {
    return gauss(n_q);
}

QuadRule<1> sinh_transform(const QuadRule<1>& gauss_rule, 
                           double a, double b, bool iterated_sinh);
QuadRule<2> sinh_sigmoidal_transform(const QuadRule<1>& gauss_theta,
    const QuadRule<1>& gauss_r, double x0, double y0, double b, bool iterated_sinh);

template <size_t dim>
struct QuadStrategy {
    QuadStrategy(size_t obs_order);
    QuadStrategy(size_t obs_order, size_t n_singular_steps,
                 double far_threshold, double near_tol);

    const QuadRule<dim-1> obs_quad;
    const QuadRule<dim-1> src_far_quad;
    
    const double far_threshold;
    const size_t n_singular_steps;
    const std::vector<double> singular_steps;
    const double near_tol;
};

} // END namespace tbem
#endif
