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

double P2M_pt_node(const OctreeCell& cell,
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

void P2M_pts_cell(FMMInfo& fmm_info, int m_cell_idx) {
    const auto cell = fmm_info.src_oct.cells[m_cell_idx];
    int cell_start_idx = m_cell_idx * fmm_info.nodes.size();
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
            int node_idx = cell_start_idx + j;   
            double P2M_kernel = P2M_pt_node(cell, fmm_info.nodes[j],
                                            fmm_info.src_oct.elements[i], 
                                            fmm_info.n_exp_pts);
            fmm_info.multipole_weights[node_idx] += fmm_info.values[i] * P2M_kernel;
        }
    }
}

void P2M_helper(FMMInfo& fmm_info, int m_cell_idx) {
    const auto cell = fmm_info.src_oct.cells[m_cell_idx];
    bool children = false;
    for (int c = 0; c < 8; c++) {
        int child_idx = cell.children[c];
        if (child_idx == -1) {
            continue;
        }
        children = true;
        P2M_helper(fmm_info, child_idx);
        auto child = fmm_info.src_oct.cells[child_idx];
        //TODO: Extract this function as P2M_cell_cell.
        for(unsigned int i = 0; i < fmm_info.nodes.size(); i++) {
            int child_node_idx = fmm_info.nodes.size() * child_idx + i;
            std::array<double,3> mapped_pt_src;
            for(int d = 0; d < 3; d++) {
                mapped_pt_src[d] = ref_to_real(fmm_info.nodes[i][d],
                                               child.bounds.min_corner[d],
                                               child.bounds.max_corner[d]);
            }
            for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
                int node_idx = fmm_info.nodes.size() * m_cell_idx + j;   
                double P2M_kernel = P2M_pt_node(cell, 
                                                fmm_info.nodes[j],
                                                mapped_pt_src,
                                                fmm_info.n_exp_pts);
                fmm_info.multipole_weights[node_idx] += 
                    fmm_info.multipole_weights[child_node_idx] * P2M_kernel;
            }
        }
    }
    // If no children, do leaf P2M
    if (!children) {
        P2M_pts_cell(fmm_info, m_cell_idx);
    }
}

void P2M(FMMInfo& fmm_info) {
    P2M_helper(fmm_info, fmm_info.src_oct.get_root_index()); 
}

void P2P_cell_pt(FMMInfo& fmm_info, int m_cell_idx, int pt_idx) {
    const auto cell = fmm_info.src_oct.cells[m_cell_idx];
    const auto obs_pt = fmm_info.obs_oct.elements[pt_idx];
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        const auto src_pt = fmm_info.src_oct.elements[i];
        const double P2P_kernel = fmm_info.kernel(obs_pt, src_pt);
        fmm_info.obs_effect[pt_idx] += fmm_info.values[i] * P2P_kernel;
    }
}

void M2P_cell_pt(FMMInfo& fmm_info, int m_cell_idx, int pt_idx) {
    const auto cell = fmm_info.src_oct.cells[m_cell_idx];
    int cell_start_idx = m_cell_idx * fmm_info.nodes.size();
    for(unsigned int j = 0; j < fmm_info.nodes.size(); j++) {
        int node_idx = cell_start_idx + j;   
        std::array<double,3> src_node_loc;
        for(int d = 0; d < 3; d++) {
            src_node_loc[d] = ref_to_real(fmm_info.nodes[j][d],
                                      cell.bounds.min_corner[d],
                                      cell.bounds.max_corner[d]);
        }
        double M2P_kernel = 
            fmm_info.kernel(fmm_info.obs_oct.elements[pt_idx], src_node_loc);
        fmm_info.obs_effect[pt_idx] += 
            fmm_info.multipole_weights[node_idx] * M2P_kernel;
    }
}

void treecode_eval_helper(FMMInfo& fmm_info, int m_cell_idx, int pt_idx) {
    auto x = fmm_info.obs_oct.elements[pt_idx];
    const auto cell = fmm_info.src_oct.cells[m_cell_idx];
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
            treecode_eval_helper(fmm_info, child_idx, pt_idx);
        }
        if (!children) {
            P2P_cell_pt(fmm_info, m_cell_idx, pt_idx);
        }
    } else {
        M2P_cell_pt(fmm_info, m_cell_idx, pt_idx);
    }
}

void treecode_eval(FMMInfo& fmm_info) {
    auto pts = fmm_info.obs_oct.elements;
    for(unsigned int i = 0; i < pts.size(); i++) {
        treecode_eval_helper(fmm_info, fmm_info.src_oct.get_root_index(), i);
    }
}

void M2L(FMMInfo& fmm_info, int m_cell_idx, int l_cell_idx);

void L2P(FMMInfo& fmm_info, int l_cell_idx, int pt_idx);


