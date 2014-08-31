#ifndef _OCTREE2_H
#define _OCTREE2_H

#include <memory>
#include <vector>
#include <iostream>
#include <limits>
#include "geom.h"


inline int to_octree_space(double x, double center, 
                           double half_width, int leaves) {
    int res = std::floor(((x - center) / (2 * half_width) + 0.5) * leaves);
    return res;
}

/* One quirk to the behavior of this octree implementation. 
 * All points must be on the interior of the octree, they cannot be on
 * the boundaries. This allows the "to_octree_space" function to ignore the
 * edge cases involving the boundaries.
 */
template <typename T, int dim> 
class Octree2 {
public:
    Octree2(std::unique_ptr<std::vector<T>> &elements,
            int depth);
    
    const static int n_octants = (int)pow(2, dim);
    const int depth;
    const int n_leaves_1d;
    const std::unique_ptr<std::vector<T>> elements;
    Box<dim> bounds;
};

template <typename T, int dim> 
Octree2<T, dim>::Octree2(std::unique_ptr<std::vector<T>> &p_elements,
                         int depth):
    depth(depth),
    n_leaves_1d((int)pow(n_octants, depth)),
    elements(std::move(p_elements))
{ 
    bounds = bounding_box(*elements);
    bounds.half_width.loc[0] *= 1.01;

    std::vector<std::array<int, dim>> leaf_indices(elements->size());
    for (int i = 0; i < elements->size(); ++i) {
        for (int d = 0; d < dim; d++) {
            leaf_indices[i][d] = to_octree_space((*elements)[i].loc[d],
                                         bounds.center.loc[d],
                                         bounds.half_width.loc[d],
                                         n_leaves_1d);
        }
    }
    
    std::vector<int>(pow(n_leaves_1d, dim)
}

// This can be removed as it is clearly completely implicit.
// However, it's nice to have in order to describe the structure 
// in 
// std::vector<std::vector<std::vector<int>>> levels(depth);
// for (int i = 0; i < depth; ++i) {
//     int sum = 0;
//     const int n_branches = (int)pow(2, i);
//     levels[i].resize(n_branches);
//     for (int j = 0; j < n_branches; ++j) {
//         levels[i][j] = naturals(sum, n_octants);
//         sum += n_octants; 
//     }
//     for (auto l: levels[i]) {
//         for (auto c: l) {
//             std::cout << c << std::endl;
//         }
//     }
// }

#endif
