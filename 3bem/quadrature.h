#ifndef __QEWKJQWEMNZ_QUADRATURE_H
#define __QEWKJQWEMNZ_QUADRATURE_H
#include <vector>
#include <utility>
#include <functional>

namespace tbem {

template <size_t dim>
struct QuadPt {
    const std::array<double,dim> x_hat;
    const double w;
};

template <size_t dim>
using QuadRule = std::vector<QuadPt<dim>>;

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, size_t dim>
T integrate(const std::vector<QuadPt<dim>>& qr, 
            const std::function<T(std::array<double,dim>)>& fnc) {
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (unsigned int i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

/* One dimensional quadrature methods */
QuadRule<1> gauss(unsigned int n);

/* Two dimensional quadrature methods */
QuadRule<2> tensor_product(QuadRule<1> xq, QuadRule<1> yq);
QuadRule<2> tensor_gauss(int n_pts);
QuadRule<2> tri_gauss(int n_pts);
QuadRule<2> square_to_tri(QuadRule<2> square_quad);

template <size_t dim>
struct QuadStrategy {
    QuadStrategy(int obs_order);
    QuadStrategy(int obs_order, int src_far_order, int n_singular_steps,
                 double far_threshold, double near_tol);

    const QuadRule<dim-1> obs_quad;
    const QuadRule<dim-1> src_far_quad;
    
    const double far_threshold;
    const int n_singular_steps;
    const std::vector<double> singular_steps;
    const double near_tol;
};

} // END namespace tbem
#endif
