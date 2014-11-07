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

template <typename T, int dim>
T integrate(const std::vector<QuadPt<dim>>& qr, 
            std::function<T(std::array<double,dim>)> fnc);

typedef std::vector<QuadPt<1>> QuadRule1d;
typedef std::vector<QuadPt<2>> QuadRule2d;

/* One dimensional quadrature methods */
QuadRule1d double_exp(int n, double h);
QuadRule1d gauss(unsigned int n);
QuadRule1d diligenti_mapping(unsigned int n, double x0, int q);

/* Two dimensional quadrature methods */
QuadRule2d tensor_product(QuadRule1d xq, QuadRule1d yq);
QuadRule2d tensor_gauss(int n_pts);
QuadRule2d tensor_double_exp(int n_pts, double h);
QuadRule2d tri_gauss(int n_pts);
QuadRule2d tri_double_exp(int n_pts, double h);
QuadRule2d tri_double_exp(int n_pts);
QuadRule2d square_to_tri(QuadRule2d square_quad);


struct QuadStrategy {
    QuadStrategy(int obs_order);
    QuadStrategy(int obs_order, int src_far_order, int src_near_order,
                 int n_singular_steps, double far_threshold);

    std::vector<const QuadRule2d*> get_near_quad(bool singular) const;

    const QuadRule2d obs_quad;

    const QuadRule2d src_far_quad;
    const QuadRule2d src_near_quad;
    
    const double far_threshold;
    const int n_singular_steps;
    const std::vector<double> singular_steps;
    const std::vector<int> singular_orders;
    const std::vector<QuadRule2d> singular_quads;
};
#endif
