#ifndef __KASDJLKJKJ_FMM_H
#define __KASDJLKJKJ_FMM_H

#include <memory>
#include "kernel.h"
#include "nbody_data.h"
#include "octree.h"
#include "numbers.h"
#include "vectorx.h"

namespace tbem {

extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv,
    int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA,
    int *IPIV, double *B, int *LDB, int *INFO);
struct LUDecomposition {
    std::vector<double> LU;
    std::vector<int> pivots;
};
LUDecomposition LU_decompose(const std::vector<double>& matrix);
std::vector<double> LU_solve(LUDecomposition& lu, const std::vector<double>& b);


template <size_t dim>
struct TranslationSurface {
    const std::vector<Vec<double,dim>> pts;
    const std::vector<Vec<double,dim>> normals;
};
template <size_t dim>
TranslationSurface<dim> make_surrounding_surface(size_t expansion_order);



struct P2MOperator {
    // LUDecomposition equiv_to_check;
    // std::vector<double> src_to_check;
};

template <size_t dim, size_t R, size_t C>
struct P2MBuilder {
    const Kernel<dim,R,C>& K;
    const Octree<dim>& src_oct;
    const NBodyData<dim>& nbody_data;
    const TranslationSurface<dim> translation_surface;
    const double d = 0.1;
    std::vector<P2MOperator> ops;

    P2MBuilder(const Kernel<dim,R,C>& K, const Octree<dim>& src_oct,
        const NBodyData<dim>& nbody_data, size_t expansion_order):
        K(K), src_oct(src_oct), nbody_data(nbody_data),
        translation_surface(make_surrounding_surface<dim>(expansion_order))
    {}

    OctreeData<dim,P2MOperator> build() {
        if (src_oct.is_leaf()) {
            std::vector<Vec<double,dim>> src_pts;
            auto n_srcs = src_oct.data.indices.size();
            auto n_chks = translation_surface.pts.size();
            std::vector<std::vector<double>> src_to_check(R * C, 
                    std::vector<double>(n_srcs * n_chks));
            for (size_t i = 0; i < n_chks; i++) {
                auto chk_pt = translation_surface.pts[i];
                auto chk_normal = translation_surface.normals[i];
                for (auto j: src_oct.data.indices) {
                    auto src_pt = nbody_data.src_locs[j];
                    auto src_normal = nbody_data.src_normals[j];

                    auto d = src_pt - chk_pt;
                    auto r2 = dot_product(d, d);
                    auto kernel_val = K(r2, d, src_normal, chk_normal);
                    auto entry = nbody_data.src_weights[j] * kernel_val;
                    for (size_t d1 = 0; d1 < R; d1++) {
                        for (size_t d2 = 0; d2 < C; d2++) {
                            src_to_check[d1 * C + d2][i * n_srcs + j] = entry[d1][d2];
                        }
                    }
                }
            }
        } else if (src_oct.data.level <= 1) {
            // do nothing
            return;    
        } else {
            // mid tree cell
            std::vector<Vec<double,dim>> src_pts;
        }
    }
};

template <size_t dim, size_t R, size_t C>
struct TreeNBodyOperator {
    const Kernel<dim,R,C>& K;
    const NBodyData<dim> data;
    Octree<dim> src_oct;
    // P2MOperator p2m;
    const double mac;

    void P2M(const BlockVectorX& x) const {}

    BlockVectorX P2P(const BlockVectorX& x) const {
        BlockVectorX out(R, VectorX(data.obs_locs.size(), 0.0));
        for (size_t i = 0; i < data.obs_locs.size(); i++) {
            for (size_t j = 0; j < data.src_locs.size(); j++) {
                auto d = data.src_locs[j] - data.obs_locs[i];
                auto r2 = dot_product(d, d);
                auto kernel_val = K(r2, d, data.src_normals[j], data.obs_normals[i]);
                auto entry = data.src_weights[j] * kernel_val;
                for (size_t d1 = 0; d1 < R; d1++) {
                    for (size_t d2 = 0; d2 < C; d2++) {
                        out[d1][i] += entry[d1][d2] * x[d2][j];
                    }
                }
            }
        }
        return out;
    }

    BlockVectorX apply(const BlockVectorX& x) const {
        // auto equivalent_srces = P2M(x);
        // auto farfield = M2P(equivalent_srces);
        auto nearfield = P2P(x);
        return nearfield;// + farfield;
    }
};

template <size_t dim, size_t R, size_t C>
TreeNBodyOperator<dim,R,C> make_tree_nbody_operator(const Kernel<dim,R,C>& K,
    const NBodyData<dim>& data, size_t min_pts_per_cell, size_t expansion_order,
    double mac)
{
    auto src_oct = build_octree(data.src_locs, min_pts_per_cell);
    auto p2m = P2MBuilder<dim,R,C>(K, src_oct, data, expansion_order).build();
    return {K, data, std::move(src_oct), mac};
}

// void FMMInfo::treecode_process_cell(const OctreeCell& cell, int cell_idx, int pt_idx) {
//     const double dist_squared = dist2({{obs_oct.elements[0][pt_idx],
//                                        obs_oct.elements[1][pt_idx],
//                                        obs_oct.elements[2][pt_idx]}},
//                                        cell.bounds.center);
//     const double radius_squared = hypot2(cell.bounds.half_width); 
//     if (dist_squared > mac2 * radius_squared) {
//         M2P_cell_pt(cell.bounds, cell_idx, pt_idx);    
//         return;
//     } else if (cell.is_leaf) {
//         P2P_cell_pt(cell, pt_idx);
//         return;
//     }
//     treecode_helper(cell, pt_idx);
// }
// 
// void FMMInfo::treecode_helper(const OctreeCell& cell, int pt_idx) {
//     for (int c = 0; c < 8; c++) {
//         const int child_idx = cell.children[c];
//         if (child_idx == -1) {
//             continue;
//         }
//         const auto child = src_oct.cells[child_idx];
//         treecode_process_cell(child, child_idx, pt_idx);
//     }
// }
// 
// void FMMInfo::treecode() {
//     const int root_idx = src_oct.get_root_index();
//     const auto root = src_oct.cells[root_idx];
// #pragma omp parallel for
//     for(unsigned int i = 0; i < obs_oct.n_elements(); i++) {
//         treecode_process_cell(root, root_idx, i);
//     }
// }

} // END namespace tbem
#endif
