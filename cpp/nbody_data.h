#ifndef __QQQQQQQQQQQQQQQZZZZZZZZZZZ_NBODY_DATA_H
#define __QQQQQQQQQQQQQQQZZZZZZZZZZZ_NBODY_DATA_H

#include "quadrature.h"
#include "mesh.h"
#include "facet_info.h"

namespace tbem {

template <size_t dim>
struct NBodyData {
    std::vector<Vec<double,dim>> obs_locs;
    std::vector<Vec<double,dim>> obs_normals;
    std::vector<Vec<double,dim>> src_locs;
    std::vector<Vec<double,dim>> src_normals;
    std::vector<double> src_weights;
};

template <size_t dim>
NBodyData<dim> nbody_data_from_bem(const Mesh<dim>& obs_mesh, const Mesh<dim>& src_mesh,
    const QuadRule<dim-1>& obs_quad, const QuadRule<dim-1>& src_quad)
{
    auto n_obs_dofs = obs_mesh.n_facets() * obs_quad.size();
    std::vector<Vec<double,dim>> obs_locs(n_obs_dofs);
    std::vector<Vec<double,dim>> obs_normals(n_obs_dofs);
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto obs_dof = obs_idx * obs_quad.size() + obs_q;
            obs_normals[obs_dof] = obs_face.normal;
            obs_locs[obs_dof] = ref_to_real(obs_quad[obs_q].x_hat, obs_face.face);
        }
    }

    auto n_src_dofs = src_mesh.n_facets() * src_quad.size();
    std::vector<Vec<double,dim>> src_locs(n_src_dofs);
    std::vector<Vec<double,dim>> src_normals(n_src_dofs);
    std::vector<double> src_weights(n_src_dofs);
    for (size_t src_idx = 0; src_idx < src_mesh.facets.size(); src_idx++) {
        auto src_face = FacetInfo<dim>::build(src_mesh.facets[src_idx]);
        for (size_t src_q = 0; src_q < src_quad.size(); src_q++) {
            auto src_dof = src_idx * src_quad.size() + src_q;
            src_weights[src_dof] = src_quad[src_q].w * src_face.jacobian;
            src_normals[src_dof] = src_face.normal;
            src_locs[src_dof] = ref_to_real(src_quad[src_q].x_hat, src_face.face);
        }
    }

    return NBodyData<dim>{obs_locs, obs_normals, src_locs, src_normals, src_weights};
}

} //end namespace tbem

#endif
