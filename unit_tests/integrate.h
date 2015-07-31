#ifndef TBEMTEST_INTEGRATE_H
#define TBEMTEST_INTEGRATE_H
#include "quad_rule.h"
#include <functional>

namespace tbem {

/* A helper function for integrating a given function using a quadrature rule.
 * Via templating, can be used with 1D, 2D, double, Vec3<double> quadrature.
 */
template <typename T, size_t dim>
T integrate(const QuadRule<dim>& qr, 
    const std::function<T(std::array<double,dim>)>& fnc) 
{
    T integral_val = qr[0].w * fnc(qr[0].x_hat);;
    for (size_t i = 1; i < qr.size(); i++) {
        integral_val += qr[i].w * fnc(qr[i].x_hat);
    }
    return integral_val;
}

} //end namespace tbem

#endif
