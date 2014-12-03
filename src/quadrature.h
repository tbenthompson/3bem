#ifndef __QEWKJQWEMNZ_QUADRATURE_H
#define __QEWKJQWEMNZ_QUADRATURE_H
#include <vector>
#include <utility>
#include <functional>

template <int dim>
struct QuadPt {
    std::array<double,dim> x_hat;
    double w;
};

template <int dim>
using QuadRule = std::vector<QuadPt<dim>>;

template <typename T, int dim>
T integrate(const QuadRule<dim>& qr,
            const std::function<T(std::array<double,dim>)>& fnc);

/* One dimensional quadrature methods */
QuadRule<1> gauss(unsigned int n);

/* Two dimensional quadrature methods */
QuadRule<2> tensor_product(QuadRule<1> xq, QuadRule<1> yq);
QuadRule<2> tensor_gauss(int n_pts);
QuadRule<2> tri_gauss(int n_pts);
QuadRule<2> square_to_tri(QuadRule<2> square_quad);

template <int dim>
struct QuadStrategy {
    QuadStrategy(int obs_order);
    QuadStrategy(int obs_order, int src_far_order, int n_singular_steps,
                 double far_threshold, double singular_tol);

    const QuadRule<dim-1> obs_quad;
    const QuadRule<dim-1> src_far_quad;
    
    const double far_threshold;
    const int n_singular_steps;
    const std::vector<double> singular_steps;
    const double singular_tol;
};
#endif
