#ifndef __FMM_H
#define __FMM_H

#include <array>
#include <vector>
#include <functional>

class OctreeCell;
class Octree;
    
typedef std::function<double (std::array<double,3>, std::array<double,3>)> Kernel;

class FMMInfo {
public:
    FMMInfo(Kernel kernel, const Octree& src, std::vector<double>& values,
            const Octree& obs, int n_exp_pts);
    // The n_exp_pts of interpolation.
    int n_exp_pts;

    // The interpolation expansion nodes in the reference cell.
    std::vector<std::array<double,3>> nodes; 

    Kernel kernel;

    // The source points and tree structure
    const Octree& src_oct;
    // The multipole expansion weights in the source cells.
    std::vector<double> multipole_weights;
    // The input point strengths.
    const std::vector<double>& values;

    // The observation points and tree structure
    const Octree& obs_oct;
    // These weight vectors should be of length N_cells * (n_exp_pts^3).
    std::vector<double> local_weights;
    // The observation effect strengths
    // obs_effect is the main output at the end of a FMM processing pass
    std::vector<double> obs_effect;
};

FMMInfo setup_fmm_info(Octree& oct, int n_exp_pts, std::vector<double> values);

//Particle to multipole
double P2M_pt_node(const OctreeCell& cell,
                   const std::array<double,3>& node,
                   const std::array<double,3>& pt,
                   int n_exp_pts);

void P2M_pts_cell(FMMInfo& fmm_info, int m_cell_idx);

// Particle to multipole -- whole upwards pass.
void P2M(FMMInfo& fmm_info);

//Point to point. 
void P2P_cell_pt(FMMInfo& fmm_info, int m_cell_idx, int pt_idx);

//Multipole to point
void M2P_cell_pt(FMMInfo& fmm_info, int m_cell_idx, int pt_idx);

void treecode_eval_helper(FMMInfo& fmm_info, int m_cell_idx, int pt_idx);
void treecode_eval(FMMInfo& fmm_info);

void M2L(FMMInfo& fmm_info, int m_cell_idx, int l_cell_idx);

void L2P(FMMInfo& fmm_info, int l_cell_idx, int pt_idx);


#endif
