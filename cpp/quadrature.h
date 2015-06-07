#ifndef TBEMQEWKJQWEMNZ_QUADRATURE_H
#define TBEMQEWKJQWEMNZ_QUADRATURE_H
#include <vector>
#include <utility>
#include <functional>
#include <map>
#include <array>

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

} // END namespace tbem
#endif
