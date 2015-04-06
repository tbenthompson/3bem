#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "dense_operator.h"
#include "block_operator.h"
#include "sparse_operator.h"
#include "galerkin_operator.h"
#include "interpolation_operator.h"
#include "obs_pt.h"
#include "integral_term.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
BlockDenseOperator mesh_to_points_operator(const std::vector<ObsPt<dim>>& obs_pts, 
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd)
{
    auto n_obs_dofs = obs_pts.size();
    auto n_src_dofs = src_mesh.n_dofs();
    std::vector<DenseOperator> ops;
    for (size_t i = 0; i < R * C; i++) {
        ops.push_back(DenseOperator(n_obs_dofs, n_src_dofs, 0.0));
    }
    BlockDenseOperator block_op{R, C, ops}; 
    auto src_facet_info = get_facet_info(src_mesh);
#pragma omp parallel for
    for (size_t pt_idx = 0; pt_idx < obs_pts.size(); pt_idx++) {
        auto pt = obs_pts[pt_idx];
        std::vector<Vec<Vec<double,C>,R>> result(src_mesh.n_dofs());
        FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
        for (size_t i = 0; i < src_mesh.facets.size(); i++) {
            IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
            auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
            auto integrals = mthd.compute_term(term, nearest_pt);
            for (int b = 0; b < dim; b++) {
                result[dim * i + b] = integrals[b];
            }
        }
        auto start_idx = pt_idx * n_src_dofs;
        for (size_t i = 0; i < result.size(); i++) {
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    block_op.ops[d1 * C + d2][start_idx + i] += result[i][d1][d2];
                }
            }
        }
    }
    return block_op;
}

template <size_t dim, size_t R, size_t C>
BlockSparseOperator galerkin_nearfield(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) 
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    auto n_blocks = R * C;

    auto src_facet_info = get_facet_info(src_mesh);
    std::vector<std::vector<SparseMatrixEntry>> entries(obs_mesh.facets.size());
    const auto& obs_quad = mthd.get_obs_quad();
#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(obs_quad[obs_q].x_hat, obs_face);

            std::vector<std::pair<size_t,Vec<Vec<double,C>,R>>> row;
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                if (nearest_pt.type == FarNearType::Farfield) {
                    continue; 
                }
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto integrals = mthd.compute_term(term, nearest_pt) -
                                mthd.compute_farfield(term, nearest_pt);
                for (int b = 0; b < dim; b++) {
                    row.push_back(std::make_pair(dim * i + b, integrals[b]));
                }
            }

            const auto basis = linear_basis(obs_quad[obs_q].x_hat);

            for (size_t obs_basis_idx = 0; obs_basis_idx < dim; obs_basis_idx++) {
                auto obs_dof = dim * obs_idx + obs_basis_idx;
                for (const auto& e: row) {
                    auto val_to_add = basis[obs_basis_idx] * e.second *
                        obs_quad[obs_q].w * obs_face.jacobian;
#pragma omp critical
                    for (size_t d1 = 0; d1 < R; d1++) {
                        for (size_t d2 = 0; d2 < C; d2++) {
                            entries[d1 * C + d2].push_back(
                                {obs_dof, e.first, val_to_add[d1][d2]}
                            );
                        }
                    }
                }
            }
        }
    }

    std::vector<SparseOperator> ops;
    for (size_t i = 0; i < n_blocks; i++) {
        ops.push_back(SparseOperator(n_obs_dofs, n_src_dofs, entries[i]));
    }
    return BlockSparseOperator(R, C, std::move(ops));
}

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

template <size_t dim, size_t R, size_t C>
struct BlockIntegralOperator: public BlockOperatorI {
    const BlockSparseOperator nearfield;
    const BlockGalerkinOperator<dim> galerkin;
    const BlockNBodyOperator<dim,R,C> farfield;
    const BlockInterpolationOperator<dim> interp;

    BlockIntegralOperator(const BlockSparseOperator& nearfield,
        const BlockGalerkinOperator<dim>& galerkin,
        const BlockNBodyOperator<dim,R,C>& farfield,
        const BlockInterpolationOperator<dim>& interp):
        nearfield(nearfield),
        galerkin(galerkin),
        farfield(farfield),
        interp(interp)
    {}

    virtual size_t n_block_rows() const {return nearfield.n_block_rows();}
    virtual size_t n_block_cols() const {return nearfield.n_block_cols();}
    virtual size_t n_total_rows() const {return nearfield.n_total_rows();} 
    virtual size_t n_total_cols() const {return nearfield.n_total_cols();}
    virtual BlockVectorX apply(const BlockVectorX& x) const {
        auto interpolated = interp.apply(x);
        auto nbodied = farfield.apply(interpolated);
        auto galerkin_far = galerkin.apply(nbodied);
        auto near_eval = nearfield.apply(x);
        return near_eval + galerkin_far;
    }
};

template <size_t dim, size_t R, size_t C>
BlockIntegralOperator<dim,R,C> integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) {

    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd);
    auto farfield = nbody_farfield(obs_mesh, src_mesh, mthd);

    auto obs_quad = mthd.get_obs_quad();
    auto src_quad = mthd.get_src_quad();
    BlockGalerkinOperator<dim> galerkin({R,C}, obs_mesh, obs_quad);
    BlockInterpolationOperator<dim> interp({R,C}, src_mesh, src_quad);

    return BlockIntegralOperator<dim,R,C>(nearfield, galerkin, farfield, interp);
}

template <size_t dim, size_t R, size_t C>
BlockDenseOperator dense_integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd)
{
    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    std::vector<DenseOperator> ops;
    for (size_t i = 0; i < R * C; i++) {
        ops.push_back(DenseOperator(n_obs_dofs, n_src_dofs, 0.0));
    }
    BlockDenseOperator block_op{R, C, ops}; 
    auto src_facet_info = get_facet_info(src_mesh);
    const auto& obs_quad = mthd.get_obs_quad();

#pragma omp parallel for
    for (size_t obs_idx = 0; obs_idx < obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FacetInfo<dim>::build(obs_mesh.facets[obs_idx]);

        std::vector<Vec<Vec<Vec<double,C>,R>,dim>> row(n_src_dofs, 
                zeros<Vec<Vec<Vec<double,C>,R>,dim>>::make());
        for (size_t obs_q = 0; obs_q < obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(obs_quad[obs_q].x_hat, obs_face);

            const auto basis = linear_basis(obs_quad[obs_q].x_hat);
            std::vector<Vec<Vec<double,C>,R>> add_to_row(src_mesh.n_dofs());
            FarNearLogic<dim> far_near_logic{mthd.far_threshold(), 1.0};
            for (size_t i = 0; i < src_mesh.facets.size(); i++) {
                IntegralTerm<dim,R,C> term{pt, src_facet_info[i]};
                auto nearest_pt = far_near_logic.decide(pt.loc, src_facet_info[i]);
                auto integrals = mthd.compute_term(term, nearest_pt);
                for (int b = 0; b < dim; b++) {
                    add_to_row[dim * i + b] = integrals[b];
                }
            }
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
                for (size_t d1 = 0; d1 < R; d1++) {
                    for (size_t d2 = 0; d2 < C; d2++) {
                        block_op.ops[d1 * C + d2][idx] += val_to_add[d1][d2];
                    }
                }
            }
        }
    }
    return block_op;
}



} // END NAMESPACE tbem
#endif
