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
        auto r_ref = 0.9;
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
           const std::vector<double>& x) 
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
        std::map<size_t,LUDecomposition>& check_to_equiv_ops) const
    {
#pragma omp critical
        if (check_to_equiv_ops.count(cell.data.level) == 0) {
            check_to_equiv_ops[cell.data.level] = build_check_to_equiv_operator(cell);
        }
        auto op = check_to_equiv_ops[cell.data.level];
        return LU_solve(op, x);
    }

    std::vector<double> apply_src_to_check_operator(const Octree<dim>& cell,
        const BlockVectorX& x) const
    {
        assert(cell.is_leaf());
        auto n_src = cell.data.indices.size();

        NBodyData<dim> s2c;
        s2c.src_locs.resize(n_src);
        s2c.src_normals.resize(n_src);
        s2c.src_weights.resize(n_src);
        s2c.obs_locs = surface.upward_check_points(cell.data.bounds, d);
        s2c.obs_normals = surface.normals;

        std::vector<double> src_str(n_src * C);
        for (size_t i = 0; i < n_src; i++) {
            s2c.src_locs[i] = data.src_locs[cell.data.indices[i]];
            s2c.src_normals[i] = data.src_normals[cell.data.indices[i]];
            s2c.src_weights[i] = data.src_weights[cell.data.indices[i]];
            for (size_t d = 0; d < C; d++) {
                src_str[d * n_src + i] = x[d][cell.data.indices[i]];
            }
        }

        return nbody_eval(K, s2c, src_str);
    }

    std::vector<double> apply_children_to_check_operator(const Octree<dim>& cell,
        const typename P2MData::ChildrenType& child_p2m) const
    {
        assert(!cell.is_leaf());

        auto check_pts = surface.upward_check_points(cell.data.bounds, d);
        auto n_check = check_pts.size();
        auto n_equiv = n_check;

        auto n_children = cell.count_children();
        auto n_src = n_children * n_check;

        NBodyData<dim> c2c;
        c2c.obs_locs = check_pts;
        c2c.obs_normals = surface.normals;
        c2c.src_locs.resize(n_src);
        c2c.src_normals.resize(n_src);
        c2c.src_weights.resize(n_src);

        std::vector<double> src_str(n_src * C);
        size_t child_idx = 0;
        for (size_t c = 0; c < Octree<dim>::split; c++) {
            auto& child = cell.children[c];
            if (child == nullptr) {
                continue;
            }
            auto equiv_pts = surface.upward_equiv_points(child->data.bounds, d);
            for (size_t i = 0; i < n_equiv; i++) {
                auto idx = child_idx * n_equiv + i;
                c2c.src_locs[idx] = equiv_pts[i];
                c2c.src_normals[idx] = surface.normals[i];
                c2c.src_weights[idx] = 1.0;
                for (size_t d = 0; d < C; d++) {
                    src_str[d * n_src + idx] = child_p2m[c]->data[d * n_equiv + i];
                }
            }
            child_idx++;
        }

        return nbody_eval(K, c2c, src_str);
    }

    std::unique_ptr<P2MData> P2M(const Octree<dim>& cell, const BlockVectorX& x,
        std::map<size_t,LUDecomposition>& check_to_equiv_ops) const
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

        auto equiv_srcs = apply_check_to_equiv_operator(
            cell, check_eval, check_to_equiv_ops
        );

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
        std::map<size_t,LUDecomposition> check_to_equiv;
        const auto p2m = P2M(src_oct, x, check_to_equiv);

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
