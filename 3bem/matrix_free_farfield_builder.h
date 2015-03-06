#ifndef __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H
#define __78987898789787_MATRIX_FREE_OPERATOR_BUILDER_H

#include "petsc_facade.h"
#include "dense_builder.h"
#include "quadrature.h"
#include "facet_info.h"

namespace tbem {

template <size_t n_rows, size_t n_cols> 
void reshape_to_add(std::vector<std::vector<MatrixEntry>>& entries, 
    size_t row, size_t col, const Vec<Vec<double,n_cols>,n_rows>& data) 
{
    for (size_t d1 = 0; d1 < n_rows; d1++) {
        for (size_t d2 = 0; d2 < n_cols; d2++) {
            entries[d1 * n_cols + d2].push_back({row, col, data[d1][d2]});
        }
    }
}

template <size_t dim, typename KT>
BlockSparseOperator build_nearfield(const BoundaryIntegral<dim,KT>& p,
    const QuadStrategy<dim>& qs) 
{
    size_t n_obs_dofs = p.obs_mesh.n_dofs();
    size_t n_src_dofs = p.src_mesh.n_dofs();
    auto n_blocks = KT::n_rows * KT::n_cols;

    auto src_facet_info = get_facet_info(p.src_mesh);
    std::vector<std::vector<MatrixEntry>> entries(p.obs_mesh.facets.size());
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

            std::vector<std::pair<size_t,typename KT::OperatorType>> row;
            for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
                FarNearLogic far_near_logic{qs.far_threshold, 3.0};
                if (far_near_logic.decide(pt.loc, src_facet_info[i]).type 
                        == FarNearType::Farfield) {
                    continue; 
                }
                auto term = make_integral_term(qs, p.K, pt, src_facet_info[i]);
                auto integrals = compute_term<dim>(term);
                for (int b = 0; b < dim; b++) {
                    row.push_back(std::make_pair(dim * i + b, integrals[b]));
                }
            }

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

            for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
                int obs_dof = dim * obs_idx + obs_basis_idx;
                for (const auto& e: row) {
                    auto val_to_add = basis[obs_basis_idx] * e.second *
                        qs.obs_quad[obs_q].w * obs_face.jacobian;
#pragma omp critical
                    reshape_to_add(entries, obs_dof, e.first, val_to_add);
                }
            }
        }
    }

    std::vector<SparseOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        ops.push_back(SparseOperator(n_obs_dofs, n_src_dofs, entries[i]));
    }
    return BlockSparseOperator(KT::n_rows, KT::n_cols, std::move(ops));
}

template <size_t dim, typename KT>
struct MatrixFreeFarfieldOperator {
    const BoundaryIntegral<dim,KT> p;
    const QuadStrategy<dim> qs;
    const BlockSparseOperator nearfield;

    MatrixFreeFarfieldOperator(const BoundaryIntegral<dim,KT>& p, const QuadStrategy<dim>& qs):
        p(p), qs(qs), nearfield(build_nearfield(p, qs))
    {}

    BlockVectorX apply(const BlockVectorX& x) {
        size_t n_obs_dofs = p.obs_mesh.n_dofs();
        BlockVectorX farfield(KT::n_rows, VectorX(n_obs_dofs, 0.0));

        auto src_facet_info = get_facet_info(p.src_mesh);
#pragma omp parallel for
        for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
            auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
            for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
                auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

                typename KT::OutType row = zeros<typename KT::OutType>::make();
                for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
                    FarNearLogic far_near_logic{qs.far_threshold, 3.0};
                    if (far_near_logic.decide(pt.loc, src_facet_info[i]).type 
                            != FarNearType::Farfield) {
                        continue; 
                    }
                    auto term = make_integral_term(qs, p.K, pt, src_facet_info[i]);
                    auto integrals = compute_term<dim>(term);
                    for (int b = 0; b < dim; b++) {
                        for (size_t d1 = 0; d1 < KT::n_rows; d1++) {
                            for (size_t d2 = 0; d2 < KT::n_cols; d2++) {
                                row[d1] += integrals[b][d1][d2] * x[d2][dim * i + b];
                            }
                        }
                    }
                }

                const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

                for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
                    int obs_dof = dim * obs_idx + obs_basis_idx;
                    auto val_to_add = basis[obs_basis_idx] * row *
                        qs.obs_quad[obs_q].w * obs_face.jacobian;
                    for (size_t d = 0; d < KT::n_rows; d++) {
                        farfield[d][obs_dof] += val_to_add[d]; 
                    }
                }
            }
        }

        return nearfield.apply(x) + farfield;
    }
};

template <size_t dim, typename KT>
MatrixFreeFarfieldOperator<dim,KT>
make_matrix_free(const BoundaryIntegral<dim,KT>& p, const QuadStrategy<dim>& qs) {
    return MatrixFreeFarfieldOperator<dim,KT>(p,qs);
}

template <size_t dim, typename KT>
BlockSparseOperator 
matrix_free_mass_operator(const BoundaryIntegral<dim,KT>& p, const QuadStrategy<dim>& qs) {
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

    std::vector<SparseOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        ops.push_back(SparseOperator(n_obs_dofs, n_obs_dofs, entries[i]));
    }
    return BlockSparseOperator(KT::n_rows, KT::n_cols, std::move(ops));
}

}//end namespace tbem

#endif
