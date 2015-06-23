#ifndef TBEMSINH_QUADRATURE_H
#define TBEMSINH_QUADRATURE_H

#include "quad_rule.h"

namespace tbem {

QuadRule<1> sinh_transform(const QuadRule<1>& gauss_rule, 
                           double a, double b, bool iterated_sinh);
QuadRule<2> sinh_sigmoidal_transform(const QuadRule<1>& gauss_theta,
    const QuadRule<1>& gauss_r, double x0, double y0, double b, bool iterated_sinh);

} // end namespace tbem

#endif
