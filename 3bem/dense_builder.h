#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "mesh.h"
#include "integral_term.h"
#include "identity_kernels.h"
#include "dense_operator.h"
#include "block_operator.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct BoundaryIntegral {
    const Mesh<dim>& obs_mesh;
    const Mesh<dim>& src_mesh;
    const Kernel<dim,R,C>& K;
};

template <size_t dim, size_t R, size_t C> 
BoundaryIntegral<dim,R,C> make_boundary_integral(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const Kernel<dim,R,C>& k) 
{
    return {obs_mesh, src_mesh, k};
}


template <size_t dim> struct ObsPt;

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
//TODO: Don't pass a BoundaryIntegral, because this function only needs the kernel
//and the src_mesh.
template <size_t dim, size_t R, size_t C>
BlockDenseOperator mesh_to_point_operator(const BoundaryIntegral<dim,R,C>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs);

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <size_t dim, size_t R, size_t C>
BlockDenseOperator mesh_to_mesh_operator(const BoundaryIntegral<dim,R,C>& p,
                                     const QuadStrategy<dim>& qs);

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <size_t dim, size_t R, size_t C>
BlockDenseOperator 
mass_operator(const BoundaryIntegral<dim,R,C>& p, const QuadStrategy<dim>& qs);

} // END NAMESPACE tbem
#endif
