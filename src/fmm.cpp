#include "fmm.h"
#include "octree.h"
#include "numerics.h"

std::vector<std::array<double,3>> get_3d_expansion_nodes(int n_exp_pts) {
    auto nodes = cheb_pts_first_kind(n_exp_pts);
    std::vector<std::array<double,3>> nodes_3d;
    for(int d0 = 0; d0 < n_exp_pts; d0++) {
        for(int d1 = 0; d1 < n_exp_pts; d1++) {
            for(int d2 = 0; d2 < n_exp_pts; d2++) {
                nodes_3d.push_back({nodes[d0], nodes[d1], nodes[d2]});
            }
        }
    }
    return nodes_3d;
}

double interp_operator(const OctreeCell& cell,
                       const std::array<double,3>& node,
                       const std::array<double,3>& pt,
                       int n_exp_pts) {
    double effect = 1.0;
    for (int d = 0; d < 3; d++) {
        double x_hat = real_to_ref(pt[d], cell.bounds.min_corner[d],
                                   cell.bounds.max_corner[d]);
        effect *= s_n(node[d], x_hat, n_exp_pts);
    }
    return effect;
}


FMMInfo::FMMInfo(Kernel kernel, const Octree& src, std::vector<double>& values, 
                 const Octree& obs, int n_exp_pts):
    n_exp_pts(n_exp_pts),
    nodes(get_3d_expansion_nodes(n_exp_pts)),
    kernel(kernel),
    src_oct(src),
    multipole_weights(src_oct.cells.size() * nodes.size(), 0.0),
    values(values),
    obs_oct(obs),
    local_weights(obs_oct.cells.size() * nodes.size(), 0.0),
    obs_effect(obs_oct.elements.size())
{}

void FMMInfo::P2M_pts_cell(int m_cell_idx) {
    const auto cell = src_oct.cells[m_cell_idx];
    int cell_start_idx = m_cell_idx * nodes.size();
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        for(unsigned int j = 0; j < nodes.size(); j++) {
            int node_idx = cell_start_idx + j;   
            double P2M_kernel = interp_operator(cell, nodes[j],
                                                src_oct.elements[i], 
                                                n_exp_pts);
            multipole_weights[node_idx] += values[i] * P2M_kernel;
        }
    }
}

void FMMInfo::P2M_helper(int m_cell_idx) {
    const auto cell = src_oct.cells[m_cell_idx];
    bool children = false;
    for (int c = 0; c < 8; c++) {
        int child_idx = cell.children[c];
        if (child_idx == -1) {
            continue;
        }
        children = true;
        // bottom-up tree traversal, recurse before doing work.
        P2M_helper(child_idx);
        auto child = src_oct.cells[child_idx];
        //TODO: Extract this function as P2M_cell_cell.
        for(unsigned int i = 0; i < nodes.size(); i++) {
            int child_node_idx = nodes.size() * child_idx + i;
            std::array<double,3> mapped_pt_src;
            for(int d = 0; d < 3; d++) {
                mapped_pt_src[d] = ref_to_real(nodes[i][d],
                                               child.bounds.min_corner[d],
                                               child.bounds.max_corner[d]);
            }
            for(unsigned int j = 0; j < nodes.size(); j++) {
                int node_idx = nodes.size() * m_cell_idx + j;   
                double P2M_kernel = interp_operator(cell, nodes[j],
                                                    mapped_pt_src, n_exp_pts);
                multipole_weights[node_idx] += 
                    multipole_weights[child_node_idx] * P2M_kernel;
            }
        }
    }
    // If no children, do leaf P2M
    if (!children) {
        P2M_pts_cell(m_cell_idx);
    }
}

void FMMInfo::P2M() {
    P2M_helper(src_oct.get_root_index()); 
}

void FMMInfo::P2P_cell_pt(int m_cell_idx, int pt_idx) {
    const auto cell = src_oct.cells[m_cell_idx];
    const auto obs_pt = obs_oct.elements[pt_idx];
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        const auto src_pt = src_oct.elements[i];
        const double P2P_kernel = kernel(obs_pt, src_pt);
        obs_effect[pt_idx] += values[i] * P2P_kernel;
    }
}

void FMMInfo::M2P_cell_pt(const Box& m_cell_bounds,
                 int m_cell_idx, int pt_idx) {
    int cell_start_idx = m_cell_idx * nodes.size();
    for(unsigned int j = 0; j < nodes.size(); j++) {
        int node_idx = cell_start_idx + j;   
        std::array<double,3> src_node_loc;
        for(int d = 0; d < 3; d++) {
            src_node_loc[d] = ref_to_real(nodes[j][d],
                                          m_cell_bounds.min_corner[d],
                                          m_cell_bounds.max_corner[d]);
        }
        double M2P_kernel = 
            kernel(obs_oct.elements[pt_idx], src_node_loc);
        obs_effect[pt_idx] += 
            multipole_weights[node_idx] * M2P_kernel;
    }
}

void FMMInfo::treecode_eval_helper(int m_cell_idx, int pt_idx) {
    auto x = obs_oct.elements[pt_idx];
    const auto cell = src_oct.cells[m_cell_idx];
    const double dist_squared = dist2<3>(x, cell.bounds.center);
    const double radius_squared = hypot2(cell.bounds.half_width); 
    if (dist_squared < 9 * radius_squared) {
        //too close
        bool children = false;
        for (int c = 0; c < 8; c++) {
            int child_idx = cell.children[c];
            if (child_idx == -1) {
                continue;
            }
            children = true;
            treecode_eval_helper(child_idx, pt_idx);
        }
        if (!children) {
            P2P_cell_pt(m_cell_idx, pt_idx);
        }
    } else {
        M2P_cell_pt(cell.bounds, m_cell_idx, pt_idx);
    }
}

void FMMInfo::treecode_eval() {
    auto pts = obs_oct.elements;
#pragma omp parallel for
    for(unsigned int i = 0; i < pts.size(); i++) {
        treecode_eval_helper(src_oct.get_root_index(), i);
    }
}
// 
// void M2L_cell_cell(FMMInfo& fmm_info, const Box& m_cell_bounds, int m_cell_idx, 
//                    const Box& l_cell_bounds, int l_cell_idx) {
//     int m_cell_start_idx = m_cell_idx * fmm_info.nodes.size();
//     int l_cell_start_idx = l_cell_idx * fmm_info.nodes.size();
//     for(unsigned int i = 0; i < fmm_info.nodes.size(); i++) {
//         int l_node_idx = l_cell_start_idx + i;
//         std::array<double,3> obs_node_loc;
//         for(int d = 0; d < 3; d++) {
//             obs_node_loc[d] = ref_to_real(fmm_info.nodes[i][d],
//                                           l_cell_bounds.min_corner[d],
//                                           l_cell_bounds.max_corner[d]);
//         }
// 
//         for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
//             int m_node_idx = m_cell_start_idx + j;   
//             std::array<double,3> src_node_loc;
//             for(int d = 0; d < 3; d++) {
//                 src_node_loc[d] = ref_to_real(fmm_info.nodes[j][d],
//                                               m_cell_bounds.min_corner[d],
//                                               m_cell_bounds.max_corner[d]);
//             }
// 
//             double M2L_kernel = fmm_info.kernel(obs_node_loc, src_node_loc);
//             fmm_info.local_weights[l_node_idx] += 
//                 fmm_info.multipole_weights[m_node_idx] * M2L_kernel;
//         }
//     }
// }
// 
// void L2P_cell_pts(FMMInfo& fmm_info, int l_cell_idx) {
//     const auto cell = fmm_info.obs_oct.cells[l_cell_idx];
//     int cell_start_idx = l_cell_idx * fmm_info.nodes.size();
//     for(unsigned int i = cell.begin; i < cell.end; i++) {
//         for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
//             int node_idx = cell_start_idx + j;   
//             double L2P_kernel = interp_operator(cell, fmm_info.nodes[j],
//                                                 fmm_info.obs_oct.elements[i], 
//                                                 fmm_info.n_exp_pts);
//             fmm_info.obs_effect[i] += fmm_info.values[node_idx] * L2P_kernel;
//         }
//     }
// }
// 
// void L2P_helper(FMMInfo& fmm_info, int l_cell_idx) {
//     const auto cell = fmm_info.obs_oct.cells[l_cell_idx];
//     bool children = false;
//     for (int c = 0; c < 8; c++) {
//         int child_idx = cell.children[c];
//         if (child_idx == -1) {
//             continue;
//         }
//         children = true;
// 
//         // top-down tree traversal, do work before recursing
//         auto child = fmm_info.obs_oct.cells[child_idx];
//         //TODO: Extract this function as L2P_cell_cell.
//         for(unsigned int i = 0; i < fmm_info.nodes.size(); i++) {
//             int child_node_idx = fmm_info.nodes.size() * child_idx + i;
//             std::array<double,3> mapped_local_pt;
//             for(int d = 0; d < 3; d++) {
//                 mapped_local_pt[d] = ref_to_real(fmm_info.nodes[i][d],
//                                                  child.bounds.min_corner[d],
//                                                  child.bounds.max_corner[d]);
//             }
//             for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
//                 int node_idx = fmm_info.nodes.size() * l_cell_idx + j;   
//                 double L2P_kernel = interp_operator(cell, 
//                                                 fmm_info.nodes[j],
//                                                 mapped_local_pt,
//                                                 fmm_info.n_exp_pts);
//                 fmm_info.local_weights[child_node_idx] += 
//                     fmm_info.local_weights[node_idx] * L2P_kernel;
//             }
//         }
//         
//         L2P_helper(fmm_info, child_idx);
//     }
//     // If no children, do leaf P2M
//     if (!children) {
//         L2P_cell_pts(fmm_info, l_cell_idx);
//     }
// }
// 
// void L2P(FMMInfo& fmm_info) {
//     L2P_helper(fmm_info, fmm_info.obs_oct.get_root_index()); 
// }
// 
// void fmm_helper(FMMInfo& fmm_info, int m_cell_idx, int l_cell_idx) {
//     const auto m_cell = fmm_info.src_oct.cells[m_cell_idx]; 
//     const auto l_cell = fmm_info.obs_oct.cells[l_cell_idx]; 
//     const double dist_squared = dist2<3>(l_cell.bounds.center, 
//                                          m_cell.bounds.center);
//     const double m_radius_squared = hypot2(m_cell.bounds.half_width); 
//     const double l_radius_squared = hypot2(l_cell.bounds.half_width); 
// 
//     // squared multipole acceptance criteria
//     double mac2 = (m_radius_squared + l_radius_squared) / dist_squared;
//     std::cout << mac2 << std::endl;
// 
//     if (mac2 < 0.25 || std::isinf(mac2)) {
//         //too close
//         bool children = false;
//         for (int c = 0; c < 8; c++) {
//             // If m_cell is more refined, refine l_cell.
//             // TODO: These if statements are ugly. Add a is_leaf flag to
//             // OctreeCell& and use a queue based implementation like in
//             // Yokota 2012 (FMM based on dtt...)
//             int child_idx;
//             if (m_radius_squared >= l_radius_squared) {
//                 child_idx = l_cell.children[c];
//             } else {
//                 child_idx = m_cell.children[c];
//             }
// 
//             if (child_idx == -1) {
//                 continue;
//             }
//             children = true;
// 
//             if (m_radius_squared >= l_radius_squared) {
//                 fmm_helper(fmm_info, m_cell_idx, child_idx);
//             } else {
//                 fmm_helper(fmm_info, child_idx, l_cell_idx);
//             }
//         }
//         if (!children) {
//             for(unsigned int i = l_cell.begin; i < l_cell.end; i++) {
//                 P2P_cell_pt(fmm_info, m_cell_idx, i);
//             }
//         }
//     } else {
//         M2L_cell_cell(fmm_info, m_cell.bounds, m_cell_idx, 
//                       l_cell.bounds, l_cell_idx);
//     }
// }
// 
// void fmm(FMMInfo& fmm_info) {
//     fmm_helper(fmm_info, fmm_info.src_oct.get_root_index(),
//                fmm_info.obs_oct.get_root_index());
// }
