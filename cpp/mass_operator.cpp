#include "mass_operator.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
BlockSparseOperator
mass_operator(const BoundaryIntegral<dim,R,C>& p, const QuadStrategy<dim>& qs)
{
    auto n_obs_dofs = p.obs_mesh.n_dofs();
    auto n_blocks = R * C;
    std::vector<std::vector<SparseMatrixEntry>> entries(n_blocks);

    auto Z = zeros<Vec<double,dim>>::make();
    auto kernel_val = p.K(0.0, Z, Z, Z);
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

template
BlockSparseOperator
mass_operator(const BoundaryIntegral<2,1,1>& p, const QuadStrategy<2>& qs);
template
BlockSparseOperator
mass_operator(const BoundaryIntegral<2,2,2>& p, const QuadStrategy<2>& qs);
template
BlockSparseOperator
mass_operator(const BoundaryIntegral<3,1,1>& p, const QuadStrategy<3>& qs);
template
BlockSparseOperator
mass_operator(const BoundaryIntegral<3,3,3>& p, const QuadStrategy<3>& qs);

} //namespace tbem
