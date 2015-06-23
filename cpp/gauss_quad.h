#ifndef TBEMGAUSS_QUAD_H
#define TBEMGAUSS_QUAD_H
#include "quad_rule.h"

namespace tbem {

QuadRule<1> gauss(size_t n);
QuadRule<2> tensor_gauss(size_t n_pts);
QuadRule<2> tri_gauss(size_t n_pts);
template <size_t dim> QuadRule<dim-1> gauss_facet(size_t n_q);

} // END namespace tbem
#endif
