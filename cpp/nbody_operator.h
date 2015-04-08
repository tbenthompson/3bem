#ifndef __JLQJWE67155151_NBODY_OPERATOR_H
#define __JLQJWE67155151_NBODY_OPERATOR_H

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct BlockNBodyOperator {
    const OperatorShape shape;
    const std::vector<Vec<double,dim>> obs_locs;
    const std::vector<Vec<double,dim>> obs_normals;

    const std::vector<Vec<double,dim>> src_locs;
    const std::vector<Vec<double,dim>> src_normals;
    const std::vector<double> src_weights;

    const Kernel<dim,R,C>& K;

    BlockVectorX apply(const BlockVectorX& x) const {
        BlockVectorX out(shape.n_rows, VectorX(obs_locs.size(), 0.0));
#pragma omp parallel for
        for (size_t i = 0; i < obs_locs.size(); i++) {
            for (size_t j = 0; j < src_locs.size(); j++) {
                auto d = src_locs[j] - obs_locs[i];
                auto r2 = dot_product(d, d);
                auto kernel_val = K(r2, d, src_normals[j], obs_normals[i]);
                auto entry = src_weights[j] * kernel_val;
                for (size_t d1 = 0; d1 < R; d1++) {
                    for (size_t d2 = 0; d2 < C; d2++) {
                        out[d1][i] += entry[d1][d2] * x[d2][j];
                    }
                }
            }
        }
        return out;
    }
};

template <size_t dim, size_t R, size_t C>
BlockNBodyOperator<dim,R,C> nbody_farfield(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) 
{
    auto obs_quad = mthd.get_obs_quad();
    auto src_quad = mthd.get_src_quad();
    const auto& K = mthd.get_kernel();

    auto n_obs_dofs = obs_mesh.n_facets() * obs_quad.size();
    std::vector<Vec<double,dim>> obs_locs(n_obs_dofs);
    std::vector<Vec<double,dim>> obs_normals(n_obs_dofs);
#pragma omp parallel for
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
#pragma omp parallel for
    for (size_t src_idx = 0; src_idx < src_mesh.facets.size(); src_idx++) {
        auto src_face = FacetInfo<dim>::build(src_mesh.facets[src_idx]);
        for (size_t src_q = 0; src_q < src_quad.size(); src_q++) {
            auto src_dof = src_idx * src_quad.size() + src_q;
            src_weights[src_dof] = src_quad[src_q].w * src_face.jacobian;
            src_normals[src_dof] = src_face.normal;
            src_locs[src_dof] = ref_to_real(src_quad[src_q].x_hat, src_face.face);
        }
    }

    return BlockNBodyOperator<dim,R,C>{
        {R, C},
        obs_locs, obs_normals, 
        src_locs, src_normals, src_weights,
        K
    };
}

} // end namespace tbem

#endif
