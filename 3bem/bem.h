#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "vec.h"
#include "numerics.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "quadrature.h"
#include "integral_term.h"
#include "identity_kernels.h"
#include "dense_operator.h"

namespace tbem {

template <size_t dim, typename KT>
struct Problem {
    const Mesh<dim>& obs_mesh;
    const Mesh<dim>& src_mesh;
    const KT& K;
};

template <size_t dim, typename KT> Problem<dim,KT> 
make_problem(const Mesh<dim>& obs_mesh, const Mesh<dim>& src_mesh, const KT& k) 
{
    return {obs_mesh, src_mesh, k};
}

template <size_t dim>
struct FacetInfo {
    //The responsibility is on the user to maintain the lifetime of the facet.
    const Facet<dim> face;

    const double area_scale;
    const double length_scale;
    const double jacobian;
    const Vec<double,dim> normal;

    static FacetInfo<dim> build(const Facet<dim>& facet);
};


template <>
inline FacetInfo<3> FacetInfo<3>::build(const Facet<3>& facet){
    const double inv_ref_facet_area = 2.0;
    auto unscaled_n = unscaled_normal(facet);
    auto area = tri_area(unscaled_n);
    auto length_scale = std::sqrt(area);
    auto jacobian = area * inv_ref_facet_area;
    auto normal = unscaled_n / jacobian;
    return FacetInfo<3>{facet, area, length_scale, jacobian, normal};
}

template <>
inline FacetInfo<2> FacetInfo<2>::build(const Facet<2>& facet){
    const double inv_ref_facet_area = 0.5;
    auto unscaled_n = unscaled_normal(facet);
    auto area_scale = hypot2(unscaled_n);
    auto length = std::sqrt(area_scale);
    auto jacobian = length * inv_ref_facet_area;
    auto normal = unscaled_n / length;
    return FacetInfo<2>{facet, area_scale, length, jacobian, normal};
}


template <size_t dim>
struct ObsPt {
    static ObsPt<dim> from_face(const Vec<double,dim-1>& ref_loc,
                                const FacetInfo<dim>& obs_face);

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};


template <size_t dim>
ObsPt<dim> ObsPt<dim>::from_face(const Vec<double,dim-1>& ref_loc,
                                 const FacetInfo<dim>& obs_face) {
    const int basis_order = 1;
    return {
        obs_face.length_scale / basis_order,
        ref_to_real(ref_loc, obs_face.face),
        obs_face.normal,
        obs_face.normal 
    };
}

/* It may be that the exact distance to an element is not required,
 * but a low precision approximate distance is useful.
 * This function simply approximates the distance by the distance from
 * the given point to each vertex of the element.
 * A better approximation to the distance to a face might include
 * the centroid (see appx_face_dist2)
 */
template <size_t dim>
double appx_face_dist2(const Vec<double,dim>& pt,
                       const std::array<Vec<double,dim>,dim>& vs) {
    double res = dist2(pt, vs[0]);
    for (int d = 1; d < dim; d++) {
        res = std::min(res, dist2(pt, vs[d])); 
    }
    return res;
}

template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> mesh_to_point_vector(const Problem<dim,KT>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs) 
{
    size_t n_out_dofs = dim * p.src_mesh.facets.size();
    std::vector<typename KT::OperatorType> result(n_out_dofs);
    for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        auto src_face = FacetInfo<dim>::build(p.src_mesh.facets[i]);
        const double dist2 = appx_face_dist2<dim>(obs.loc, src_face.face);
        auto term = make_integral_term(qs, p.K, obs, src_face, dist2);
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

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
template <size_t dim, typename KT>
BlockDenseOperator mesh_to_point_operator(const Problem<dim,KT>& p,
    const QuadStrategy<dim>& qs, const ObsPt<dim>& obs) 
{
    size_t n_out_dofs = dim * p.src_mesh.facets.size();
    auto result = mesh_to_point_vector(p, qs, obs);
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
BlockDenseOperator mesh_to_mesh_operator(const Problem<dim,KT>& p,
                                     const QuadStrategy<dim>& qs) 
{
    size_t n_obs_dofs = p.obs_mesh.n_dofs();
    size_t n_src_dofs = p.src_mesh.n_dofs();
    auto block_op = build_operator_shape(KT::n_cols, KT::n_rows, n_obs_dofs, n_src_dofs);
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

            const auto row = mesh_to_point_vector(p, qs, pt);
            assert(row.size() == n_src_dofs);

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

            for (int obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
                int obs_dof = dim * obs_idx + obs_basis_idx;
                for (size_t src_dof = 0; src_dof < n_src_dofs; src_dof++) {
                    auto val_to_add = basis[obs_basis_idx] * row[src_dof] *
                        qs.obs_quad[obs_q].w * obs_face.jacobian;
                    auto idx = obs_dof * n_src_dofs + src_dof;
                    reshape_to_add(block_op, idx, val_to_add);
                }
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
BlockDenseOperator mass_operator(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs)
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

template <size_t dim>
double get_len_scale(Mesh<dim>& mesh, int which_face, int q);

template <>
inline double get_len_scale<3>(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face])) / q;
}

template <>
inline double get_len_scale<2>(Mesh<2>& mesh, int which_face, int q) {
    return dist(mesh.facets[which_face][1],
                mesh.facets[which_face][0]) / q;
}



} // END NAMESPACE tbem
#endif
