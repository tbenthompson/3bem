#ifndef __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H
#define __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H

#include "matrix_free_operator.h"
#include "dense_operator_builder.h"
#include "quadrature.h"
#include "facet_info.h"

namespace tbem {

template <size_t n_rows, size_t n_cols> 
void reshape_to_add(std::vector<std::vector<MatrixEntry>>& entries, size_t row, size_t col,
    const Vec<Vec<double,n_cols>,n_rows>& data) 
{
    for (size_t d1 = 0; d1 < n_rows; d1++) {
        for (size_t d2 = 0; d2 < n_cols; d2++) {
            entries[d1 * n_cols + d2].push_back({row, col, data[d1][d2]});
        }
    }
}

struct MatrixFreeBuilder {

};

template <size_t dim, typename KT>
BlockMatrixFreeOperator 
matrix_free_mass_operator(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs) {
    auto n_obs_dofs = p.obs_mesh.n_dofs();
    auto n_blocks = KT::n_rows * KT::n_cols;
    std::vector<std::vector<MatrixEntry>> entries(n_blocks);

    auto kernel_val = p.K.call_with_no_params();
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) 
    {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) 
        {
            auto qpt = qs.obs_quad[obs_q];
            auto basis = linear_basis(qpt.x_hat);
            auto weight = obs_face.jacobian * qpt.w;

            for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) 
            {
                int obs_dof = dim * obs_idx + obs_basis_idx;
                for (size_t src_basis_idx = 0; src_basis_idx < dim; src_basis_idx++) 
                {
                    int src_dof = dim * obs_idx + src_basis_idx;
                    auto basis_product = basis[obs_basis_idx] * basis[src_basis_idx];
                    auto entry_value = kernel_val * basis_product * weight;
                    reshape_to_add(entries, obs_dof, src_dof, entry_value);
                }
            }
        }
    }

    std::vector<MatrixFreeOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        MatrixFreeOperator::NearfieldPtr nearfield(
            new PETScSparseMatWrapper(n_obs_dofs, n_obs_dofs, entries[i]));
        ops.push_back(MatrixFreeOperator(std::move(nearfield)));
    }
    return BlockOperator<MatrixFreeOperator>(KT::n_rows, KT::n_cols, std::move(ops));
}

}//end namespace tbem

#endif
