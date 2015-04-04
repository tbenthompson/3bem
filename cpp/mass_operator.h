#ifndef __123786172627689192_MASS_OPERATOR_H
#define __123786172627689192_MASS_OPERATOR_H

#include "sparse_operator.h"
#include "block_operator.h"
#include "mesh.h"
#include "kernel.h"
#include "quadrature.h"

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
mass_operator(const Mesh<dim>& obs_mesh, size_t n_q)
{
    auto quad = gauss_facet<dim>(n_q);
    auto n_obs_dofs = obs_mesh.n_dofs();
    auto n_blocks = R * C;
    std::vector<std::vector<SparseMatrixEntry>> entries(n_blocks);

    auto kernel_val = ones<Vec<Vec<double,C>,R>>::make();
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) 
    {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < quad.size(); obs_q++) 
        {
            auto qpt = quad[obs_q];
            auto basis = linear_basis(qpt.x_hat);
            auto weight = obs_face.jacobian * qpt.w;

            for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) 
            {
                auto obs_dof = dim * obs_idx + obs_basis_idx;
                for (size_t src_basis_idx = 0; src_basis_idx < dim; src_basis_idx++) 
                {
                    auto src_dof = dim * obs_idx + src_basis_idx;
                    auto basis_product = basis[obs_basis_idx] * basis[src_basis_idx];
                    auto entry_value = kernel_val * basis_product * weight;
                    for (size_t d1 = 0; d1 < R; d1++) {
                        for (size_t d2 = 0; d2 < C; d2++) {
                            auto val = entry_value[d1][d2];
                            auto blk = d1 * C + d2;
                            entries[blk].push_back({obs_dof, src_dof, val});
                        }
                    }
                }
            }
        }
    }

    std::vector<SparseOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        ops.push_back(SparseOperator(n_obs_dofs, n_obs_dofs, entries[i]));
    }
    return BlockSparseOperator(R, C, std::move(ops));
}

} // END NAMESPACE tbem
#endif
