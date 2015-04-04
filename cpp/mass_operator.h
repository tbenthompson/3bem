#ifndef __123786172627689192_MASS_OPERATOR_H
#define __123786172627689192_MASS_OPERATOR_H
#include "dense_builder.h"
#include "sparse_operator.h"
#include "block_operator.h"

namespace tbem {

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <size_t dim, size_t R, size_t C>
BlockSparseOperator 
mass_operator(const BoundaryIntegral<dim,R,C>& p, const QuadStrategy<dim>& qs);

} // END NAMESPACE tbem
#endif
