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

namespace tbem {
template <size_t dim, typename KT>
struct Problem {
    const Mesh<dim> src_mesh;
    const Mesh<dim> obs_mesh;
    const KT K;
};

template <size_t dim, typename KT>
Problem<dim,KT> make_problem(const Mesh<dim>& src_mesh,
                             const Mesh<dim>& obs_mesh, const KT& k) 
{
    return {src_mesh, obs_mesh, k};
}

template <size_t dim>
struct FacetInfo {
    //The responsibility is on the user to maintain the lifetime of the facet.
    const Facet<dim>& face;

    const double area_scale;
    const double length_scale;
    const double jacobian;
    const Vec<double,dim> normal;

    static const double inv_ref_facet_area;
    static FacetInfo<dim> build(const Facet<dim>& facet);
};

template <> const double FacetInfo<3>::inv_ref_facet_area = 2.0;
template <> const double FacetInfo<2>::inv_ref_facet_area = 0.5;

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
        ref_to_real(ref_loc, obs_face.face.vertices),
        obs_face.normal,
        obs_face.normal 
    };
}

template <>
FacetInfo<3> FacetInfo<3>::build(const Facet<3>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area = tri_area(unscaled_n);
    auto length_scale = std::sqrt(area);
    auto jacobian = area * inv_ref_facet_area;
    auto normal = unscaled_n / jacobian;
    return FacetInfo<3>{facet, area, length_scale, jacobian, normal};
}

template <>
FacetInfo<2> FacetInfo<2>::build(const Facet<2>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area_scale = hypot2(unscaled_n);
    auto length = std::sqrt(area_scale);
    auto jacobian = length * inv_ref_facet_area;
    auto normal = unscaled_n / length;
    return FacetInfo<2>{facet, area_scale, length, jacobian, normal};
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

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> 
integral_equation_vector(const Problem<dim,KT>& p, const QuadStrategy<dim>& qs,
                         const ObsPt<dim>& obs) {
    int n_out_dofs = dim * p.src_mesh.facets.size();
    std::vector<typename KT::OperatorType> result(n_out_dofs);
    for (size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        auto src_face = FacetInfo<dim>::build(p.src_mesh.facets[i]);
        const double dist2 = appx_face_dist2<dim>(obs.loc, src_face.face.vertices);
        auto term = make_integral_term(qs, p.K, obs, src_face, dist2);
        auto integrals = compute_term<dim>(term);
        for (int b = 0; b < dim; b++) {
            result[dim * i + b] = integrals[b];
        }
    }
    return result;
}

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <size_t dim, typename KT>
std::vector<typename KT::OperatorType> interact_matrix(const Problem<dim,KT>& p,
                                    const QuadStrategy<dim>& qs) 
{
    size_t n_obs_dofs = dim * p.obs_mesh.facets.size();
    size_t n_src_dofs = dim * p.src_mesh.facets.size();
    std::vector<typename KT::OperatorType> matrix(n_obs_dofs * n_src_dofs, 
            zeros<typename KT::OperatorType>::make());
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad[obs_q].x_hat, obs_face);

            const auto row = integral_equation_vector(p, qs, pt);

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

            for (int v = 0; v < dim; v++) {
                int b = dim * obs_idx + v;
                for (size_t i = 0; i < n_src_dofs; i++) {
                    matrix[b * n_src_dofs + i] +=
                        basis[v] *
                        row[i] *
                        qs.obs_quad[obs_q].w *
                        obs_face.jacobian;
                }
            }
        }
    }
    return matrix;
}


/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <size_t dim, typename KT>
std::vector<typename KT::OutType> mass_term(const Problem<dim,KT>& p,
    const QuadStrategy<dim>& qs, const std::vector<typename KT::InType> function) 
{
    int n_obs_dofs = dim * p.obs_mesh.facets.size();
    std::vector<typename KT::OutType> integrals(
        n_obs_dofs, 
        zeros<typename KT::OutType>::make()
    );

    for (size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];

            int dof = dim * obs_idx;
            Vec<typename KT::InType,dim> face_vals;
            for (int d = 0; d < dim; d++) {
                face_vals[d] = function[dof + d];
            }

            auto interp_val = dot_product(linear_basis(qpt.x_hat), face_vals);
            auto kernel_val = p.K.call_with_no_params();
            auto out_val = dot_product(interp_val, kernel_val);

            auto basis = linear_basis(qpt.x_hat);
            for(int v = 0; v < dim; v++) {
                integrals[dof + v] += obs_face.jacobian * basis[v] * 
                                      out_val * qpt.w;
            }
        }
    }
    return integrals;
}

template <typename KT>
std::vector<typename KT::OutType>
bem_mat_mult(const std::vector<typename KT::OperatorType>& A, const KT& k,
             int n_obs_dofs, const std::vector<typename KT::InType>& x) {

    assert(n_obs_dofs * x.size() == A.size());
    std::vector<typename KT::OutType> res(n_obs_dofs, 
                                          zeros<typename KT::OutType>::make());
#pragma omp parallel for
    for (int i = 0; i < n_obs_dofs; i++) {
        for (size_t j = 0; j < x.size(); j++) {
            res[i] += dot_product(x[j], A[i * x.size() + j]);
        }
    }
    return res;
}

template <size_t dim>
double get_len_scale(Mesh<dim>& mesh, int which_face, int q);

template <>
double get_len_scale<3>(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face].vertices)) / q;
}

template <>
double get_len_scale<2>(Mesh<2>& mesh, int which_face, int q) {
    return dist(mesh.facets[which_face].vertices[1],
                mesh.facets[which_face].vertices[0]) / q;
}


} // END NAMESPACE tbem
#endif
