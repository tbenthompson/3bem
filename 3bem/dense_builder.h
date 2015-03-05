#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "mesh.h"
#include "integral_term.h"
#include "identity_kernels.h"
#include "dense_operator.h"
#include "block_operator.h"
#include "obs_pt.h"
#include "facet_info.h"

namespace tbem {

//TODO: To REALLY REALLY clean this up, ever, I'm going to need to get rid of the 
//templating on KT. This may simply involve pushing the templating into a different
//layer of compilation process. Or alternatively, I could template on the tensor 
//rows and columns?

//TODO: Perform the inner integral for each observation point and then apply the
//galerkin operator?

//TODO: PointSetI - Interface for a bunch of observation points?

//TODO: Is factory method the appropriate pattern here?
//Should each operator factory return an abstract operator interface
//without any knowledge of what the object's operator type actually is? 
//or abstract factory?

struct BasisPlaceholder {};
template <size_t dim, typename KT>
struct BoundaryIntegral {
    const Mesh<dim>& obs_mesh;
    const Mesh<dim>& src_mesh;
    const KT& K;
    const BasisPlaceholder src_basis;
};

template <size_t dim, typename KT> 
BoundaryIntegral<dim,KT> make_boundary_integral(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const KT& k) 
{
    return {obs_mesh, src_mesh, k};
}

template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> mesh_to_point_vector(const BoundaryIntegral<dim,KT>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs, 
    const std::vector<FacetInfo<dim>>& facet_info) 
{
    std::vector<typename KT::OperatorType> result(p.src_mesh.n_dofs());
    for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        auto term = make_integral_term(qs, p.K, obs, facet_info[i]);
        auto integrals = compute_term<dim>(term);
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

template <size_t dim>
std::vector<FacetInfo<dim>> get_facet_info(const Mesh<dim>& m) {
    std::vector<FacetInfo<dim>> out;
    for (const auto& f: m.facets) {
        out.push_back(FacetInfo<dim>::build(f));
    }
    return out;
}

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
//TODO: Don't pass a BoundaryIntegral, because this function only needs the kernel
//and the src_mesh.
template <size_t dim, typename KT>
BlockDenseOperator mesh_to_point_operator(const BoundaryIntegral<dim,KT>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs) 
{
    size_t n_out_dofs = dim * p.src_mesh.facets.size();
    auto result = mesh_to_point_vector(p, qs, obs, get_facet_info(p.src_mesh));
    auto block_op = build_operator_shape(KT::n_rows, KT::n_cols, 1, n_out_dofs);
    for (size_t i = 0; i < result.size(); i++) {
        reshape_to_add(block_op, i, result[i]);
    }
    return block_op;
}

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <size_t dim, typename KT>
BlockDenseOperator mesh_to_mesh_operator(const BoundaryIntegral<dim,KT>& p,
                                     const QuadStrategy<dim>& qs) 
{
    size_t n_obs_dofs = p.obs_mesh.n_dofs();
    size_t n_src_dofs = p.src_mesh.n_dofs();
    auto block_op = build_operator_shape(
        KT::n_cols, KT::n_rows, n_obs_dofs, n_src_dofs
    );
    auto src_facet_info = get_facet_info(p.src_mesh);
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);

        std::vector<Vec<typename KT::OperatorType,dim>> row(n_src_dofs, 
                zeros<Vec<typename KT::OperatorType,dim>>::make());
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);
            auto add_to_row = mesh_to_point_vector(p, qs, pt, src_facet_info);
            for (size_t dof = 0; dof < n_src_dofs; dof++) {
                row[dof] += outer_product(basis,
                    add_to_row[dof] * qs.obs_quad[obs_q].w * obs_face.jacobian);
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

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <size_t dim, typename KT>
BlockDenseOperator mass_operator(const BoundaryIntegral<dim,KT>& p, const QuadStrategy<dim>& qs)
{
    auto n_obs_dofs = p.obs_mesh.n_dofs();
    auto block_op = build_operator_shape(KT::n_rows, KT::n_cols, n_obs_dofs, n_obs_dofs);

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
                    auto idx = obs_dof * n_obs_dofs + src_dof;
                    reshape_to_add(block_op, idx, entry_value);
                }
            }
        }
    }

    return block_op;
}



} // END NAMESPACE tbem
#endif
