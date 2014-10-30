#ifndef __QEWKJQWEMNZ_QUADRATURE_H
#define __QEWKJQWEMNZ_QUADRATURE_H
#include <vector>
#include <utility>
#include <functional>

typedef std::vector<std::pair<double, double>> QuadratureRule;

QuadratureRule double_exp(int n, double h);
QuadratureRule gauss(unsigned int n);
QuadratureRule diligenti_mapping(unsigned int n, double x0, int q);
double integrate(QuadratureRule& qr, std::function<double (double)> fnc);

struct QuadratureRule2D {
    QuadratureRule2D(int n): x_hat(n), y_hat(n), weights(n) {}
    std::vector<double> x_hat;
    std::vector<double> y_hat;
    std::vector<double> weights;
};


QuadratureRule2D tensor_product(QuadratureRule xq, QuadratureRule yq);
QuadratureRule2D tensor_gauss(int n_pts);
QuadratureRule2D tensor_double_exp(int n_pts, double h);
QuadratureRule2D tri_gauss(int n_pts);
QuadratureRule2D tri_double_exp(int n_pts, double h);
QuadratureRule2D square_to_tri(QuadratureRule2D square_quad);

double integrate(QuadratureRule2D& qr, std::function<double (double,double)> fnc);
#endif
