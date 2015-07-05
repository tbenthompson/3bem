#ifndef TBEMJLQJWE67155151_NBODY_OPERATOR_H
#define TBEMJLQJWE67155151_NBODY_OPERATOR_H

#include <cstdlib>
#include <vector>
#include "vec.h"
#include "kernel.h"
#include "operator.h"
#include "dense_operator.h"
#include "mesh.h"
#include "quad_rule.h"
#include "numerics.h"

namespace tbem {

template <size_t dim>
struct NBodyObservationPoints {
    std::vector<Vec<double,dim>> locs;
    std::vector<Vec<double,dim>> normals;
};

template <size_t dim>
struct NBodySourcePoints {
    std::vector<Vec<double,dim>> locs;
    std::vector<Vec<double,dim>> normals;
    std::vector<double> weights;
};

template <size_t dim>
struct NBodyData {
    std::vector<Vec<double,dim>> obs_locs;
    std::vector<Vec<double,dim>> obs_normals;
    std::vector<Vec<double,dim>> src_locs;
    std::vector<Vec<double,dim>> src_normals;
    std::vector<double> src_weights;
};

template <size_t dim>
NBodyObservationPoints<dim> nbody_obs_from_bem(const Mesh<dim>& obs_mesh,
    const QuadRule<dim-1>& obs_quad)
{
    auto n_obs_dofs = obs_mesh.n_facets() * obs_quad.size();
    std::vector<Vec<double,dim>> obs_locs(n_obs_dofs);
    std::vector<Vec<double,dim>> obs_normals(n_obs_dofs);
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = obs_mesh.facets[obs_idx];
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto obs_dof = obs_idx * obs_quad.size() + obs_q;
            obs_normals[obs_dof] = facet_normal<dim>(obs_face);
            obs_locs[obs_dof] = ref_to_real(obs_quad[obs_q].x_hat, obs_face);
        }
    }
    return {obs_locs, obs_normals};
}

template <size_t dim>
NBodySourcePoints<dim> nbody_src_from_bem(const Mesh<dim>& src_mesh,
    const QuadRule<dim-1>& src_quad)
{
    auto n_src_dofs = src_mesh.n_facets() * src_quad.size();
    std::vector<Vec<double,dim>> src_locs(n_src_dofs);
    std::vector<Vec<double,dim>> src_normals(n_src_dofs);
    std::vector<double> src_weights(n_src_dofs);
    for (size_t src_idx = 0; src_idx < src_mesh.facets.size(); src_idx++) {
        auto src_face = src_mesh.facets[src_idx];
        for (size_t src_q = 0; src_q < src_quad.size(); src_q++) {
            auto src_dof = src_idx * src_quad.size() + src_q;
            auto jacobian = facet_jacobian<dim>(src_face);
            src_weights[src_dof] = src_quad[src_q].w * jacobian;
            src_normals[src_dof] = facet_normal<dim>(src_face);
            src_locs[src_dof] = ref_to_real(src_quad[src_q].x_hat, src_face);
        }
    }

    return {src_locs, src_normals, src_weights};
}

template <size_t dim>
NBodyData<dim> nbody_data_from_bem(const Mesh<dim>& obs_mesh, const Mesh<dim>& src_mesh,
    const QuadRule<dim-1>& obs_quad, const QuadRule<dim-1>& src_quad)
{
    auto src = nbody_src_from_bem(src_mesh, src_quad);
    auto obs = nbody_obs_from_bem(obs_mesh, obs_quad);
    return NBodyData<dim>{obs.locs, obs.normals, src.locs, src.normals, src.weights};
}

template <size_t dim, size_t R, size_t C>
std::vector<double>
nbody_matrix(const Kernel<dim,R,C>& K, const NBodyData<dim>& data, bool parallel = false) 
{
    auto n_pairs = data.obs_locs.size() * data.src_locs.size();
    auto n_blocks = R * C;
    std::vector<double> op(n_pairs * n_blocks);

#pragma omp parallel for if(parallel)
    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = data.src_weights[j] * K(
                data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]
            );

            for (size_t d1 = 0; d1 < R; d1++) {
                auto row = d1 * data.obs_locs.size() + i;
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto col = d2 * data.src_locs.size() + j;
                    auto matrix_idx = row * C * data.src_locs.size() + col;
                    op[matrix_idx] = kernel_val[d1][d2];
                }
            }
        }
    }

    return op;
}

template <size_t dim, size_t R, size_t C>
std::vector<double>
nbody_eval(const Kernel<dim,R,C>& K, const NBodyData<dim>& data,
           double const* x) 
{
    std::vector<double> out(R * data.obs_locs.size(), 0.0);
    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = data.src_weights[j] * K(
                data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]
            );

            for (size_t d1 = 0; d1 < R; d1++) {
                auto row = d1 * data.obs_locs.size() + i;
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto col = d2 * data.src_locs.size() + j;
                    out[row] += kernel_val[d1][d2] * x[col];
                }
            }
        }
    }

    return out;
}

template <size_t dim, size_t R, size_t C>
DenseOperator make_direct_nbody_operator(const NBodyData<dim>& data,
    const Kernel<dim,R,C>& K)
{
    return DenseOperator(
        R * data.obs_locs.size(),
        C * data.src_locs.size(),
        nbody_matrix(K, data, true)
    );
}

} // end namespace tbem

#endif
