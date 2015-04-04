#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "dense_operator.h"
#include "obs_pt.h"
#include "integral_term.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
std::vector<Vec<Vec<double,C>,R>> 
mesh_to_point_vector(const BoundaryIntegral<dim,R,C>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs, 
    const std::vector<FacetInfo<dim>>& facet_info) 
{
    std::vector<Vec<Vec<double,C>,R>> result(p.src_mesh.n_dofs());
    FarNearLogic<dim> far_near_logic{qs.far_threshold, 1.0};
    AdaptiveIntegrationMethod<dim,R,C> mthd(qs, p.K);
    for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        IntegralTerm<dim,R,C> term{obs, facet_info[i]};
        auto nearest_pt = far_near_logic.decide(obs.loc, facet_info[i]);
        auto integrals = mthd.compute_term(term, nearest_pt);
        for (int b = 0; b < dim; b++) {
            result[dim * i + b] = integrals[b];
        }
    }
    return result;
}

inline BlockDenseOperator build_operator_shape(size_t n_comp_rows, size_t n_comp_cols,
    size_t n_rows, size_t n_cols) 
{
    std::vector<DenseOperator> ops;
    for (size_t i = 0; i < n_comp_rows * n_comp_cols; i++) {
        ops.push_back(DenseOperator(n_rows, n_cols, 0.0));
    }
    return {n_comp_rows, n_comp_cols, ops}; 
}

template <size_t n_rows, size_t n_cols> 
void reshape_to_add(BlockDenseOperator& block_op, size_t idx,
    const Vec<Vec<double,n_cols>,n_rows>& data) 
{
    for (size_t d1 = 0; d1 < n_rows; d1++) {
        for (size_t d2 = 0; d2 < n_cols; d2++) {
            block_op.ops[d1 * n_cols + d2][idx] += data[d1][d2];
        }
    }
}

template <size_t dim, size_t R, size_t C>
BlockDenseOperator mesh_to_points_operator(const std::vector<ObsPt<dim>>& obs, 
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& int_mthd)
{
    size_t n_out_dofs = dim * p.src_mesh.facets.size();
    auto block_op = build_operator_shape(R, C, obs_pts.size(), n_out_dofs);
#pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < obs_pts.size(); pt_idx++) {
        auto pt = obs_pts[pt_idx];
        auto result = mesh_to_point_vector(p, qs, pt, get_facet_info(p.src_mesh));
        auto start_idx = pt_idx * n_out_dofs;
        for (size_t i = 0; i < result.size(); i++) {
            reshape_to_add(block_op, start_idx + i, result[i]);
        }
    }
    return block_op;
}

template <size_t dim, size_t R, size_t C>
BlockDenseOperator mesh_to_mesh_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& int_mthd)
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    auto block_op = build_operator_shape(
        R, C, n_obs_dofs, n_src_dofs
    );
    auto src_facet_info = get_facet_info(src_mesh);
    const auto& obs_quad = int_mthd.get_obs_quad();
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);

        std::vector<Vec<Vec<Vec<double,C>,R>,dim>> row(n_src_dofs, 
                zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make());
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(obs_quad[obs_q].x_hat, obs_face);

            const auto basis = linear_basis(obs_quad[obs_q].x_hat);
            auto add_to_row = mesh_to_point_vector(p, qs, pt, src_facet_info);
            for (size_t dof = 0; dof < n_src_dofs; dof++) {
                row[dof] += outer_product(basis,
                    add_to_row[dof] * obs_quad[obs_q].w * obs_face.jacobian);
            }
        }


        for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
            int obs_dof = dim * obs_idx + obs_basis_idx;
            for (size_t src_dof = 0; src_dof < n_src_dofs; src_dof++) {
                auto val_to_add = row[src_dof][obs_basis_idx];
                auto idx = obs_dof * n_src_dofs + src_dof;
                reshape_to_add(block_op, idx, val_to_add);
            }
        }
    }
    return block_op;
}



} // END NAMESPACE tbem
#endif
