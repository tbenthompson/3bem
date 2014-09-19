#include "fmm.h"
#include "octree.h"
#include "numerics.h"

std::vector<double> P2M(Octree& oct, int n_nodes, std::vector<double> p_values)
{
    // int level = oct.max_depth - 1;
    auto cheb_nodes = cheb_pts_first_kind(n_nodes);
    return cheb_nodes;

    // int total_nodes = (int)(pow(oct.levels[level].n_cells_1d, 3) * pow(order, 3));
    // std::vector<double> weights(total_nodes);
    // for (int d1 = 0; d1 < oct.levels[level].n_cells_1d; d1++) {
    // for (int d2 = 0; d2 < oct.levels[level].n_cells_1d; d2++) {
    // for (int d3 = 0; d3 < oct.levels[level].n_cells_1d; d3++) {
    //     int cell_start = threed_to_1d(d1, d2, d3, oct.levels[level].n_cells_1d);
    //     for (int n1 = 0; n1 < order; n1++) {
    //     for (int n2 = 0; n2 < order; n2++) {
    //     for (int n3 = 0; n3 < order; n3++) {
    //         int abc = 0;
    //     }
    //     }
    //     }
    // }
    // }
    // }
    
}
