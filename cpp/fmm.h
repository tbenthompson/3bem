#ifndef __KASDJLKJKJ_FMM_H
#define __KASDJLKJKJ_FMM_H

#include <cassert>
#include <memory>
#include "kernel.h"
#include "nbody_data.h"
#include "octree.h"
#include "numbers.h"
#include "vectorx.h"
#include "util.h"

namespace tbem {

extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv,
    int* info);
extern "C" void dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA,
    int* IPIV, double* B, int* LDB, int* INFO);

struct LUDecomposition {
    std::vector<double> LU;
    std::vector<int> pivots;
};
LUDecomposition LU_decompose(const std::vector<double>& matrix);
std::vector<double> LU_solve(LUDecomposition& lu, const std::vector<double>& b);

// TODO: Deal with the normal equations
// struct NormalEqtns {
//     LUDecomposition LU_ATA;
//     std::vector<double> AT;
// };
// std::vector<double> form_regularized_normal_eqtns(const std::vector<double>& matrix);

template <size_t dim>
struct TranslationSurface {
    const std::vector<Vec<double,dim>> pts;
    const std::vector<Vec<double,dim>> normals;

    std::vector<Vec<double,dim>> move(const Box<dim>& box, double r_ref) const;

    std::vector<Vec<double,dim>> 
    upward_check_points(const Box<dim>& box, double d) const
    {
        auto r_ref = 4.0 - std::sqrt(2) - 2 * d;
        return move(box, r_ref);
    }

    std::vector<Vec<double,dim>> 
    upward_equiv_points(const Box<dim>& box, double d) const
    {
        auto r_ref = std::sqrt(2) + d;
        return move(box, r_ref);
    }
};
template <size_t dim>
TranslationSurface<dim> make_surrounding_surface(size_t expansion_order);

template <size_t dim, size_t R, size_t C>
std::vector<double>
nbody_matrix(const Kernel<dim,R,C>& K, const NBodyData<dim>& data) 
{
    auto n_pairs = data.obs_locs.size() * data.src_locs.size();
    auto n_blocks = R * C;
    std::vector<double> op(n_pairs * n_blocks);

    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = K(
                data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]
            );

            auto pair_idx = i * data.src_locs.size() + j;
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto block_idx = d1 * C + d2;
                    auto matrix_idx = block_idx * n_pairs + pair_idx;
                    op[matrix_idx] = kernel_val[d1][d2];
                }
            }
        }
    }

    return op;
}

static size_t m2p = 0;
static size_t p2p = 0;

template <size_t dim, size_t R, size_t C>
struct TreeNBodyOperator {
    typedef OctreeData<dim,std::vector<double>> P2MData;

    const Kernel<dim,R,C>& K;
    const NBodyData<dim> data;
    const TranslationSurface<dim> surface;
    Octree<dim> src_oct;
    const double mac;
    const double d = 0.1;
    const double alpha = 1e-8;

    TreeNBodyOperator(const Kernel<dim,R,C>& K, const NBodyData<dim>& data,
        size_t min_pts_per_cell, size_t expansion_order, double mac):
        K(K),
        data(data),
        surface(make_surrounding_surface<dim>(expansion_order)),
        src_oct(build_octree(data.src_locs, min_pts_per_cell)),
        mac(mac)
    {}

    LUDecomposition build_check_to_equiv_operator(const Octree<dim>& cell) const
    {
        auto check_pts = surface.upward_check_points(cell.data.bounds, d);
        auto equiv_pts = surface.upward_equiv_points(cell.data.bounds, d);

        auto op = nbody_matrix(K, {
            check_pts, surface.normals,
            equiv_pts, surface.normals, {}
        });

        auto lu = LU_decompose(op);        

        return lu;
    }

    std::vector<double> apply_check_to_equiv_operator(const Octree<dim>& cell,
        const std::vector<double>& x, 
        std::map<size_t,LUDecomposition> check_to_equiv_ops) const
    {
        if (check_to_equiv_ops.count(cell.data.level) == 0) {
#pragma omp critical
            {
                check_to_equiv_ops[cell.data.level] = 
                    build_check_to_equiv_operator(cell);
            }
        }
        auto op = check_to_equiv_ops[cell.data.level];
        return LU_solve(op, x);
    }

    std::vector<double> apply_src_to_check_operator(const Octree<dim>& cell,
        const BlockVectorX& x) const
    {
        assert(cell.is_leaf());
        auto check_pts = surface.upward_check_points(cell.data.bounds, d);
        auto n_check = check_pts.size();
        std::vector<double> out(R * n_check, 0.0);
        for (size_t i = 0; i < n_check; i++) {
            for (auto j: cell.data.indices) {
                auto kernel_val = K(
                    check_pts[i], data.src_locs[j],
                    surface.normals[i], data.src_normals[j]
                );

                auto entry = data.src_weights[j] * kernel_val;
                for (size_t d1 = 0; d1 < R; d1++) {
                    for (size_t d2 = 0; d2 < C; d2++) {
                        out[d1 * n_check + i] += entry[d1][d2] * x[d2][j];
                    }
                }
            }
        }

        return out;
    }

    std::vector<double> apply_children_to_check_operator(const Octree<dim>& cell,
        const typename P2MData::ChildrenType& children) const
    {
        assert(!cell.is_leaf());
        auto check_pts = surface.upward_check_points(cell.data.bounds, d);
        auto n_check = check_pts.size();
        std::vector<double> out(R * check_pts.size(), 0.0);
        for (size_t i = 0; i < n_check; i++) {
            for (size_t c = 0; c < Octree<dim>::split; c++) {
                if (children[c] == nullptr) {
                    continue;
                }
                auto equiv_pts = surface.upward_equiv_points(
                    cell.children[c]->data.bounds, d);
                auto equiv_srcs = children[c]->data;
                auto n_equiv = equiv_pts.size();

                for (size_t j = 0; j < n_equiv; j++) {
                    auto kernel_val = K(
                        check_pts[i],
                        equiv_pts[j],
                        surface.normals[i], 
                        surface.normals[j]
                    );

                    auto entry = data.src_weights[j] * kernel_val;
                    for (size_t d1 = 0; d1 < R; d1++) {
                        for (size_t d2 = 0; d2 < C; d2++) {
                            auto obs_idx = d1 * n_check + i;
                            auto src_idx = d2 * n_equiv + j;
                            out[obs_idx] += entry[d1][d2] * equiv_srcs[src_idx];
                        }
                    }
                }
            }
        }

        return out;
    }

    std::unique_ptr<P2MData> P2M(const Octree<dim>& cell, const BlockVectorX& x,
        std::map<size_t,LUDecomposition> check_to_equiv_ops) const
    {
        std::vector<double> check_eval;
        typename P2MData::ChildrenType child_P2M;

        if (cell.is_leaf()) {
            check_eval = apply_src_to_check_operator(cell, x);
        } else {
#pragma omp parallel for if(cell.data.level == 0)
            for (size_t i = 0; i < Octree<dim>::split; i++) {
                if (cell.children[i] == nullptr) {
                    continue;
                }
                child_P2M[i] = P2M(*cell.children[i], x, check_to_equiv_ops);
            }
            check_eval = apply_children_to_check_operator(cell, child_P2M);
        }

        auto equiv_srcs = apply_check_to_equiv_operator(cell, check_eval, check_to_equiv_ops);

        return std::unique_ptr<P2MData>(new P2MData{
            equiv_srcs, std::move(child_P2M)
        });
    }

    Vec<double,R> M2P(const Vec<double,dim>& pt, const Vec<double,dim>& normal,
        const Octree<dim>& cell, const P2MData& p2m) const 
    {
        auto out = zeros<Vec<double,R>>::make();
        auto equiv_pts = surface.upward_equiv_points(cell.data.bounds, d);
        auto n_equiv = equiv_pts.size();
        for (size_t j = 0; j < n_equiv; j++) {
            auto kernel_val = K(
                pt, equiv_pts[j],
                normal, surface.normals[j]
            );
            auto entry = kernel_val;
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    out[d1] += entry[d1][d2] * p2m.data[d2 * n_equiv + j];
                }
            }
        }
        return out;
    }

    Vec<double,R> P2P(const Vec<double,dim>& pt, const Vec<double,dim>& normal,
        const Octree<dim>& cell, const BlockVectorX& x) const 
    {
        auto out = zeros<Vec<double,R>>::make();
        for (auto j: cell.data.indices) {
            auto kernel_val = K(
                pt, data.src_locs[j],
                normal, data.src_normals[j]
            );
            auto entry = data.src_weights[j] * kernel_val;
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    out[d1] += entry[d1][d2] * x[d2][j];
                }
            }
        }
        return out;
    }

    Vec<double,R> treecode(const Vec<double,dim>& pt,
        const Vec<double,dim>& normal, const Octree<dim>& cell,
        const BlockVectorX& x, const P2MData& p2m) const 
    {
        const double distance = dist(cell.data.bounds.center, pt);
        const double radius = hypot(cell.data.bounds.half_width); 
        if (distance > mac * radius) {
            m2p++;
            return M2P(pt, normal, cell, p2m);
        } else if (cell.is_leaf()) {
            p2p++;
            return P2P(pt, normal, cell, x);
        } else {
            auto out = zeros<Vec<double,R>>::make();
            for (size_t c = 0; c < Octree<dim>::split; c++) {
                const auto& child = cell.children[c];
                if (child == nullptr) {
                    continue;
                }
                out += treecode(pt, normal, *child, x, *p2m.children[c]);
            }
            return out;
        }
    }

    BlockVectorX apply(const BlockVectorX& x) const 
    {
        assert(x.size() == C);
        const auto p2m = P2M(src_oct, x, {});

        BlockVectorX out(R, VectorX(data.obs_locs.size()));
#pragma omp parallel for
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            auto res = treecode(
                data.obs_locs[i], data.obs_normals[i], src_oct, x, *p2m
            );
            for (size_t d1 = 0; d1 < R; d1++) {
                out[d1][i] = res[d1];
            }
        }
        return out;
    }
};

} // END namespace tbem
#endif
