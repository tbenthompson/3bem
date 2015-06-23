#ifndef TBEMQEWKJQWEMNZ_QUADRATURE_H
#define TBEMQEWKJQWEMNZ_QUADRATURE_H
#include <utility>
#include <map>
#include "quad_rule.h"

namespace tbem {

QuadRule<1> gauss(size_t n);
QuadRule<2> tensor_gauss(size_t n_pts);
QuadRule<2> tri_gauss(size_t n_pts);


template <size_t dim>
QuadRule<dim-1> gauss_facet(size_t n_q);

//TODO: Move to cpp
template <>
inline QuadRule<2> gauss_facet<3>(size_t n_q) 
{
    return tri_gauss(n_q);
}

template <>
inline QuadRule<1> gauss_facet<2>(size_t n_q) 
{
    return gauss(n_q);
}

} // END namespace tbem
#endif
