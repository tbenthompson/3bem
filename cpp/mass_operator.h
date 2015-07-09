#ifndef TBEM123786172627689192_MASS_OPERATOR_H
#define TBEM123786172627689192_MASS_OPERATOR_H

#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "sparse_operator.h"
#include "gauss_quad.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
SparseOperator
mass_operator(const Mesh<dim>& obs_mesh, size_t n_q)
{
    auto obs_quad = gauss_facet<dim>(n_q);
    auto interp = make_interpolation_operator(C, obs_mesh, obs_quad);
    auto galerkin = make_galerkin_operator(R, obs_mesh, obs_quad);
    return galerkin.right_multiply(interp);
}

} // END NAMESPACE tbem
#endif
