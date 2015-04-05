#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include "dense_operator.h"
#include "block_operator.h"
#include "sparse_operator.h"
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
                auto integrals = mthd.compute_term(term, nearest_pt);
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

// template <size_t dim>
// struct BlockNBodyOperator: public BlockOperatorI {
//     const OperatorShape shape;
//     const std::vector<Vec<double,dim>> obs_locs;
//     const std::vector<Vec<double,dim>> src_locs;
//     const std::function<double(const Vec<double,dim>&, const Vec<double,dim>&)> kernel;
// 
//     virtual size_t n_block_rows() const {return shape.n_rows;}
//     virtual size_t n_block_cols() const {return shape.n_cols;}
//     virtual size_t n_total_rows() const {return shape.n_rows * obs_locs.size();}
//     virtual size_t n_total_cols() const {return shape.n_cols * src_locs.size();}
// 
//     virtual BlockVectorX apply(const BlockVectorX& x) const {
// 
//     }
// };
// 
// template <typename T>
// struct BlockOperator;
// typedef BlockOperator<IntegralOperator> IntegralOperator;

struct BlockIntegralOperator: public BlockOperatorI {
    const BlockSparseOperator nearfield;
    const BlockDenseOperator farfield;

    BlockIntegralOperator(const BlockSparseOperator& nearfield,
        const BlockDenseOperator& farfield):
        nearfield(nearfield),
        farfield(farfield)
    {}

    virtual size_t n_block_rows() const {return farfield.n_block_rows();}
    virtual size_t n_block_cols() const {return farfield.n_block_cols();}
    virtual size_t n_total_rows() const {return farfield.n_total_rows();} 
    virtual size_t n_total_cols() const {return farfield.n_total_cols();}
    virtual BlockVectorX apply(const BlockVectorX& x) const {
        return nearfield.apply(x) + farfield.apply(x);
    }
};

template <size_t dim, size_t R, size_t C>
BlockIntegralOperator integral_operator(const Mesh<dim>& obs_mesh,
    const Mesh<dim>& src_mesh, const IntegrationMethodI<dim,R,C>& mthd) {

    size_t n_obs_dofs = obs_mesh.n_dofs();
    size_t n_src_dofs = src_mesh.n_dofs();
    std::vector<DenseOperator> ops;
    for (size_t i = 0; i < R * C; i++) {
        ops.push_back(DenseOperator(n_obs_dofs, n_src_dofs, 0.0));
    }
    BlockDenseOperator block_op{R, C, ops}; 
    auto src_facet_info = get_facet_info(src_mesh);
    const auto& obs_quad = mthd.get_obs_quad();

    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd);

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
                if (nearest_pt.type != FarNearType::Farfield) {
                    continue;
                }
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
    return BlockIntegralOperator{nearfield, block_op};
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

    auto nearfield = galerkin_nearfield(obs_mesh, src_mesh, mthd);
    for (size_t blk_idx = 0; blk_idx < nearfield.ops.size(); blk_idx++) {
        auto& dst_blk = block_op.ops[blk_idx];
        const auto& src_blk = nearfield.ops[blk_idx];
        for (size_t i = 0; i < src_blk.storage.size(); i++) {
            const auto& entry = src_blk.storage[i];
            dst_blk[entry.loc[0] * dst_blk.n_cols() + entry.loc[1]] += entry.value;
        }
    }

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
                if (nearest_pt.type != FarNearType::Farfield) {
                    continue;
                }
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
