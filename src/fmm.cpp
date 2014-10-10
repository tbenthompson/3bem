#include "fmm.h"
#include "octree.h"
#include "numerics.h"
#include "direct.h"
#include <cassert>

std::array<std::vector<double>,3> get_3d_expansion_nodes(int n_exp_pts) {
    auto nodes = cheb_pts_first_kind(n_exp_pts);
    std::array<std::vector<double>,3> nodes_3d;
    for(int d0 = 0; d0 < n_exp_pts; d0++) {
        for(int d1 = 0; d1 < n_exp_pts; d1++) {
            for(int d2 = 0; d2 < n_exp_pts; d2++) {
                nodes_3d[0].push_back(nodes[d0]);
                nodes_3d[1].push_back(nodes[d1]);
                nodes_3d[2].push_back(nodes[d2]);
            }
        }
    }
    return nodes_3d;
}

double interp_operator(const Box& bounds,
                       double nodex, double nodey, double nodez,
                       double ptx, double pty, double ptz,
                       int n_exp_pts) {
    double effect = 1.0;
    // X interp
    double x_hat = real_to_ref(ptx, bounds.min_corner[0], bounds.max_corner[0]);
    effect *= s_n_fast(nodex, x_hat, n_exp_pts);
    // Y interp
    x_hat = real_to_ref(pty, bounds.min_corner[1], bounds.max_corner[1]);
    effect *= s_n_fast(nodey, x_hat, n_exp_pts);
    // Z interp
    x_hat = real_to_ref(ptz, bounds.min_corner[2], bounds.max_corner[2]);
    effect *= s_n_fast(nodez, x_hat, n_exp_pts);
    return effect;
}


FMMInfo::FMMInfo(Kernel kernel, const Octree& src, std::vector<double>& values, 
                 const Octree& obs, int n_exp_pts, double mac2):
    mac2(mac2),
    n_exp_pts(n_exp_pts),
    nodes(get_3d_expansion_nodes(n_exp_pts)),
    kernel(kernel),
    src_oct(src),
    multipole_weights(src_oct.cells.size() * nodes[0].size(), 0.0),
    values(values),
    obs_oct(obs),
    local_weights(obs_oct.cells.size() * nodes[0].size(), 0.0),
    obs_effect(obs_oct.elements[0].size(), 0.0),
    p2p_jobs(obs_oct.cells.size()),
    m2l_jobs(obs_oct.cells.size())
{}

void FMMInfo::P2M_pts_cell(int m_cell_idx) {
    const auto cell = src_oct.cells[m_cell_idx];
    int cell_start_idx = m_cell_idx * nodes[0].size();
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        for(unsigned int j = 0; j < nodes[0].size(); j++) {
            int node_idx = cell_start_idx + j;   
            double P2M_kernel = interp_operator(cell.bounds, nodes[0][j],
                                                nodes[1][j], nodes[2][j],
                                                src_oct.elements[0][i], 
                                                src_oct.elements[1][i],
                                                src_oct.elements[2][i],
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
        for(unsigned int i = 0; i < nodes[0].size(); i++) {
            int child_node_idx = nodes[0].size() * child_idx + i;
            std::array<double,3> mapped_pt_src;
            for(int d = 0; d < 3; d++) {
                mapped_pt_src[d] = ref_to_real(nodes[d][i],
                                               child.bounds.min_corner[d],
                                               child.bounds.max_corner[d]);
            }
            for(unsigned int j = 0; j < nodes[0].size(); j++) {
                int node_idx = nodes[0].size() * m_cell_idx + j;   
                double P2M_kernel = interp_operator(cell.bounds,
                                                    nodes[0][j],
                                                    nodes[1][j],
                                                    nodes[2][j],
                                                    mapped_pt_src[0],
                                                    mapped_pt_src[1],
                                                    mapped_pt_src[2],
                                                    n_exp_pts);
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
    for(unsigned int i = m_cell.begin; i < m_cell.end; i++) {
        const double P2P_kernel = kernel(obs_oct.elements[0][pt_idx],
                                         obs_oct.elements[1][pt_idx],
                                         obs_oct.elements[2][pt_idx],
                                         src_oct.elements[0][i],
                                         src_oct.elements[1][i],
                                         src_oct.elements[2][i]
                                         );
        obs_effect[pt_idx] += values[i] * P2P_kernel;
    }
}

void FMMInfo::P2P_cell_cell(const OctreeCell& m_cell, const OctreeCell& l_cell) {
    // for(unsigned int i = l_cell.begin; i < l_cell.end; i++) {
    //     P2P_cell_pt(m_cell, i);
    // }
    vec_direct_n_body(src_oct.elements, obs_oct.elements,
                      m_cell.begin, m_cell.end, 
                      l_cell.begin, l_cell.end, values);
}

void FMMInfo::M2P_cell_pt(const Box& m_cell_bounds,
                 int m_cell_idx, int pt_idx) {
    int cell_start_idx = m_cell_idx * nodes[0].size();
    for(unsigned int j = 0; j < nodes[0].size(); j++) {
        int node_idx = cell_start_idx + j;   
        std::array<double,3> src_node_loc;
        for(int d = 0; d < 3; d++) {
            src_node_loc[d] = ref_to_real(nodes[d][j],
                                          m_cell_bounds.min_corner[d],
                                          m_cell_bounds.max_corner[d]);
        }
        const double M2P_kernel = kernel(obs_oct.elements[0][pt_idx],
                                         obs_oct.elements[1][pt_idx],
                                         obs_oct.elements[2][pt_idx],
                                         src_node_loc[0],
                                         src_node_loc[1],
                                         src_node_loc[2]
                                         );
        obs_effect[pt_idx] += multipole_weights[node_idx] * M2P_kernel;
    }
}

void FMMInfo::treecode_process_cell(const OctreeCell& cell, int cell_idx, int pt_idx) {
    const double dist_squared = dist2(obs_oct.elements[0][pt_idx],
                                       obs_oct.elements[1][pt_idx],
                                       obs_oct.elements[2][pt_idx],
                                       cell.bounds.center[0],
                                       cell.bounds.center[1],
                                       cell.bounds.center[2]);
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
    const int root_idx = src_oct.get_root_index();
    const auto root = src_oct.cells[root_idx];
#pragma omp parallel for
    for(unsigned int i = 0; i < obs_oct.n_elements(); i++) {
        treecode_process_cell(root, root_idx, i);
    }
}

void FMMInfo::M2L_cell_cell(const Box& m_cell_bounds, int m_cell_idx, 
                            const Box& l_cell_bounds, int l_cell_idx) {
    int m_cell_start_idx = m_cell_idx * nodes[0].size();
    int l_cell_start_idx = l_cell_idx * nodes[0].size();
    for(unsigned int i = 0; i < nodes[0].size(); i++) {
        int l_node_idx = l_cell_start_idx + i;
        std::array<double,3> obs_node_loc;
        for(int d = 0; d < 3; d++) {
            obs_node_loc[d] = ref_to_real(nodes[d][i],
                                          l_cell_bounds.min_corner[d],
                                          l_cell_bounds.max_corner[d]);
        }

        for(unsigned int j = 0; j < nodes[0].size(); j++) {
            int m_node_idx = m_cell_start_idx + j;   
            std::array<double,3> src_node_loc;
            for(int d = 0; d < 3; d++) {
                src_node_loc[d] = ref_to_real(nodes[d][j],
                                              m_cell_bounds.min_corner[d],
                                              m_cell_bounds.max_corner[d]);
            }

            double M2L_kernel = kernel(obs_node_loc[0], obs_node_loc[1],
                                       obs_node_loc[2], src_node_loc[0],
                                       src_node_loc[1], src_node_loc[2]);
            local_weights[l_node_idx] += 
                multipole_weights[m_node_idx] * M2L_kernel;
        }
    }
}

void FMMInfo::L2P_cell_pts(int l_cell_idx) {
    const auto cell = obs_oct.cells[l_cell_idx];
    int cell_start_idx = l_cell_idx * nodes[0].size();
    for(unsigned int i = cell.begin; i < cell.end; i++) {
        for(unsigned int j = 0; j < nodes[0].size(); j++) {
            int node_idx = cell_start_idx + j;   
            double L2P_kernel = interp_operator(cell.bounds,
                                                nodes[0][j],
                                                nodes[1][j],
                                                nodes[2][j],
                                                obs_oct.elements[0][i],
                                                obs_oct.elements[1][i],
                                                obs_oct.elements[2][i],
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
        for(unsigned int i = 0; i < nodes[0].size(); i++) {
            int child_node_idx = nodes[0].size() * child_idx + i;
            std::array<double,3> mapped_local_pt;
            for(int d = 0; d < 3; d++) {
                mapped_local_pt[d] = ref_to_real(nodes[d][i],
                                                 child.bounds.min_corner[d],
                                                 child.bounds.max_corner[d]);
            }
            for(unsigned int j = 0; j < nodes[0].size(); j++) {
                int node_idx = nodes[0].size() * l_cell_idx + j;   
                double L2P_kernel = interp_operator(cell.bounds,
                                                    nodes[0][j],
                                                    nodes[1][j],
                                                    nodes[2][j],
                                                    mapped_local_pt[0],
                                                    mapped_local_pt[1],
                                                    mapped_local_pt[2],
                                                    n_exp_pts);
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
        m2l_jobs[l_cell_idx].push_back(m_cell_idx);
        return;
    } else if (m_cell.is_leaf && l_cell.is_leaf) {
        // P2P_cell_cell(m_cell, l_cell);
        p2p_jobs[l_cell_idx].push_back(m_cell_idx);
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

bool fmm_compare(std::array<int,2> a, std::array<int,2> b) {
    if(a[1] == b[1]) {
        return a[0] < b[0];
    }
    return a[1] < b[1];
}

void FMMInfo::fmm_exec_jobs() {
    //TODO: make this a templated function and make the P2P_cell_cell
    //and M2L_cell_cell interfaces uniform
    std::vector<std::vector<int>>& job_set = p2p_jobs;
#pragma omp parallel for
    for (unsigned int l_idx = 0; l_idx < job_set.size(); l_idx++) {
        for (unsigned int j = 0; j < job_set[l_idx].size(); j++) {
            int m_idx = job_set[l_idx][j];
            P2P_cell_cell(src_oct.cells[m_idx], obs_oct.cells[l_idx]);
        }
    }

    job_set = m2l_jobs;
#pragma omp parallel for
    for (unsigned int l_idx = 0; l_idx < job_set.size(); l_idx++) {
        for (unsigned int j = 0; j < job_set[l_idx].size(); j++) {
            int m_idx = job_set[l_idx][j];
            M2L_cell_cell(src_oct.cells[m_idx].bounds, m_idx,
                          obs_oct.cells[l_idx].bounds, l_idx);
        }
    }
}
