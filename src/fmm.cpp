#include "fmm.h"
#include "octree.h"
#include "numerics.h"
#include <cassert>

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
        // Interpolation operators are computed in reference [-1, 1] space, so
        // the point must be linearly transformed to that range.
        double x_hat = real_to_ref(pt[d], 
                                   cell.bounds.min_corner[d],
                                   cell.bounds.max_corner[d]);
        effect *= s_n_fast(node[d], x_hat, n_exp_pts);
    }
    return effect;
}


FMMInfo::FMMInfo(Kernel kernel, const Octree& src, std::vector<double>& values, 
                 const Octree& obs, int n_exp_pts, double mac2):
    mac2(mac2),
    n_exp_pts(n_exp_pts),
    nodes(get_3d_expansion_nodes(n_exp_pts)),
    kernel(kernel),
    src_oct(src),
    multipole_weights(src_oct.cells.size() * nodes.size(), 0.0),
    values(values),
    obs_oct(obs),
    local_weights(obs_oct.cells.size() * nodes.size(), 0.0),
    obs_effect(obs_oct.elements.size(), 0.0)
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
    // If no children, do leaf P2M
    if (cell.is_leaf) {
        P2M_pts_cell(m_cell_idx);
        return;
    }
    for (int c = 0; c < 8; c++) {
        int child_idx = cell.children[c];
        if (child_idx == -1) {
            continue;
        }
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
}

void FMMInfo::P2M() {
    P2M_helper(src_oct.get_root_index()); 
}

void FMMInfo::P2P_cell_pt(const OctreeCell& m_cell, int pt_idx) {
    const auto obs_pt = obs_oct.elements[pt_idx];
    for(unsigned int i = m_cell.begin; i < m_cell.end; i++) {
        const auto src_pt = src_oct.elements[i];
        const double P2P_kernel = kernel(obs_pt, src_pt);
        obs_effect[pt_idx] += values[i] * P2P_kernel;
    }
}

void FMMInfo::P2P_cell_cell(const OctreeCell& m_cell, const OctreeCell& l_cell) {
    for(unsigned int i = l_cell.begin; i < l_cell.end; i++) {
        P2P_cell_pt(m_cell, i);
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
        obs_effect[pt_idx] += multipole_weights[node_idx] * M2P_kernel;
    }
}

void FMMInfo::treecode_process_cell(const OctreeCell& cell, int cell_idx, int pt_idx) {
    const auto x = obs_oct.elements[pt_idx];
    const double dist_squared = dist2<3>(x, cell.bounds.center);
    const double radius_squared = hypot2(cell.bounds.half_width); 
    if (dist_squared > mac2 * radius_squared) {
        M2P_cell_pt(cell.bounds, cell_idx, pt_idx);    
        return;
    } else if (cell.is_leaf) {
        P2P_cell_pt(cell, pt_idx);
        return;
    }
    treecode_helper(cell, pt_idx);
}

void FMMInfo::treecode_helper(const OctreeCell& cell, int pt_idx) {
    for (int c = 0; c < 8; c++) {
        const int child_idx = cell.children[c];
        if (child_idx == -1) {
            continue;
        }
        const auto child = src_oct.cells[child_idx];
        treecode_process_cell(child, child_idx, pt_idx);
    }
}

void FMMInfo::treecode() {
    const auto pts = obs_oct.elements;
    const int root_idx = src_oct.get_root_index();
    const auto root = src_oct.cells[root_idx];
#pragma omp parallel for
    for(unsigned int i = 0; i < pts.size(); i++) {
        treecode_process_cell(root, root_idx, i);
    }
}

void FMMInfo::M2L_cell_cell(const Box& m_cell_bounds, int m_cell_idx, 
                            const Box& l_cell_bounds, int l_cell_idx) {
    int m_cell_start_idx = m_cell_idx * nodes.size();
    int l_cell_start_idx = l_cell_idx * nodes.size();
    for(unsigned int i = 0; i < nodes.size(); i++) {
        int l_node_idx = l_cell_start_idx + i;
        std::array<double,3> obs_node_loc;
        for(int d = 0; d < 3; d++) {
            obs_node_loc[d] = ref_to_real(nodes[i][d],
                                          l_cell_bounds.min_corner[d],
                                          l_cell_bounds.max_corner[d]);
        }

        for(unsigned int j = 0; j < nodes.size(); j++) {
            int m_node_idx = m_cell_start_idx + j;   
            std::array<double,3> src_node_loc;
            for(int d = 0; d < 3; d++) {
                src_node_loc[d] = ref_to_real(nodes[j][d],
                                              m_cell_bounds.min_corner[d],
                                              m_cell_bounds.max_corner[d]);
            }

            double M2L_kernel = kernel(obs_node_loc, src_node_loc);
            local_weights[l_node_idx] += 
                multipole_weights[m_node_idx] * M2L_kernel;
        }
    }
}

void FMMInfo::L2P_cell_pts(int l_cell_idx) {
    const auto cell = obs_oct.cells[l_cell_idx];
    int cell_start_idx = l_cell_idx * nodes.size();
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        for(unsigned int j = 0; j < nodes.size(); j++) {
            int node_idx = cell_start_idx + j;   
            double L2P_kernel = interp_operator(cell, nodes[j],
                                                obs_oct.elements[i], 
                                                n_exp_pts);
            obs_effect[i] += local_weights[node_idx] * L2P_kernel;
        }
    }
}

void FMMInfo::L2P_helper(int l_cell_idx) {
    const auto cell = obs_oct.cells[l_cell_idx];
    // If no children, do leaf P2M
    if (cell.is_leaf) {
        L2P_cell_pts(l_cell_idx);
        return;
    }

    for (int c = 0; c < 8; c++) {
        int child_idx = cell.children[c];
        if (child_idx == -1) {
            continue;
        }

        // top-down tree traversal, recurse after doing work
        auto child = obs_oct.cells[child_idx];
        //TODO: Extract this function as L2P_cell_cell.
        for(unsigned int i = 0; i < nodes.size(); i++) {
            int child_node_idx = nodes.size() * child_idx + i;
            std::array<double,3> mapped_local_pt;
            for(int d = 0; d < 3; d++) {
                mapped_local_pt[d] = ref_to_real(nodes[i][d],
                                                 child.bounds.min_corner[d],
                                                 child.bounds.max_corner[d]);
            }
            for(unsigned int j = 0; j < nodes.size(); j++) {
                int node_idx = nodes.size() * l_cell_idx + j;   
                double L2P_kernel = interp_operator(cell, nodes[j],
                                            mapped_local_pt, n_exp_pts);
                local_weights[child_node_idx] += 
                    local_weights[node_idx] * L2P_kernel;
            }
        }
        
        L2P_helper(child_idx);
    }
}

void FMMInfo::L2P() {
    L2P_helper(obs_oct.get_root_index()); 
}

void FMMInfo::fmm_process_cell_pair(const OctreeCell& m_cell, int m_cell_idx,
                           const OctreeCell& l_cell, int l_cell_idx) {
    const double dist_squared = dist2<3>(l_cell.bounds.center, 
                                         m_cell.bounds.center);
    const double m_radius_squared = hypot2(m_cell.bounds.half_width); 
    const double l_radius_squared = hypot2(l_cell.bounds.half_width); 
    if (dist_squared > mac2 * (m_radius_squared + l_radius_squared) / 2.0) {
        // M2L_cell_cell(m_cell.bounds, m_cell_idx, l_cell.bounds, l_cell_idx); 
        m2l_jobs.push_back({m_cell_idx, l_cell_idx});
        return;
    } else if (m_cell.is_leaf && l_cell.is_leaf) {
        // P2P_cell_cell(m_cell, l_cell);
        p2p_jobs.push_back({m_cell_idx, l_cell_idx});
        return;
    }
    fmm_process_children(m_cell, m_cell_idx, l_cell, l_cell_idx);
}

void FMMInfo::fmm_process_children(const OctreeCell& m_cell, int m_cell_idx,
                                   const OctreeCell& l_cell, int l_cell_idx) {
    if ((m_cell.level > l_cell.level && !l_cell.is_leaf) || m_cell.is_leaf) {
        //refine l_cell
        assert(!l_cell.is_leaf);
        for (int c = 0; c < 8; c++) {
            const int l_child_idx = l_cell.children[c];
            if (l_child_idx == -1) {
                continue;
            }
            const auto l_child = obs_oct.cells[l_child_idx];
            fmm_process_cell_pair(m_cell, m_cell_idx, l_child, l_child_idx);
        }
    } else {
        //refine m_cell
        assert(!m_cell.is_leaf);
        for (int c = 0; c < 8; c++) {
            const int m_child_idx = m_cell.children[c];
            if (m_child_idx == -1) {
                continue;
            }
            const auto m_child = src_oct.cells[m_child_idx];
            fmm_process_cell_pair(m_child, m_child_idx, l_cell, l_cell_idx);
        }
    }
}

void FMMInfo::fmm() {
    const int src_root_idx = src_oct.get_root_index();
    const auto src_root = src_oct.cells[src_root_idx];
    const int obs_root_idx = obs_oct.get_root_index();
    const auto obs_root = obs_oct.cells[obs_root_idx];
    fmm_process_cell_pair(src_root, src_root_idx, obs_root, obs_root_idx);
    fmm_exec_jobs();
}

void FMMInfo::fmm_exec_jobs() {
    std::cout << "P2P Jobs: " << p2p_jobs.size() << std::endl;
    std::cout << "M2L Jobs: " << m2l_jobs.size() << std::endl;
    std::sort(p2p_jobs.begin(), p2p_jobs.end(),
              [] (std::array<int,2> a, std::array<int,2> b) {
                if(a[0] == b[0]) {
                    return a[1] < b[1];
                }
                return a[0] < b[0];
              });

    std::sort(m2l_jobs.begin(), m2l_jobs.end(),
              [] (std::array<int,2> a, std::array<int,2> b) {
                if(a[0] == b[0]) {
                    return a[1] < b[1];
                }
                return a[0] < b[0];
              });
// #pragma omp parallel for
    for (unsigned int i = 0; i < p2p_jobs.size(); i++) {
        auto p2p = p2p_jobs[i];
        std::cout << p2p[0] << " " << p2p[1] << std::endl;
        auto m_cell = src_oct.cells[p2p[0]];
        auto l_cell = obs_oct.cells[p2p[1]];
        P2P_cell_cell(m_cell, l_cell);
    }

// #pragma omp parallel for
    for (unsigned int i = 0; i < m2l_jobs.size(); i++) {
        auto m2l = m2l_jobs[i];
        auto m_cell = src_oct.cells[m2l[0]];
        auto l_cell = obs_oct.cells[m2l[1]];
        M2L_cell_cell(m_cell.bounds, m2l[0], l_cell.bounds, m2l[1]);
    }
}
