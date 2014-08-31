#include "octree.h"

Box bounding_box(const std::array<std::vector<double>,3>& x)
{
    Vec3 min_corner = {x[0][0], x[1][0], x[2][0]};
    Vec3 max_corner = {x[0][0], x[1][0], x[2][0]};
    for (unsigned int i = 1; i < x.size(); ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            min_corner.loc[d] = std::min(min_corner.loc[d], x[d][i]);
            max_corner.loc[d] = std::max(max_corner.loc[d], x[d][i]);
        }
    }

    Box ext;
    ext.center = (max_corner + min_corner) / 2.0;
    ext.half_width = (max_corner - min_corner) / 2.0;
    return ext;
}

std::vector<int> naturals(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
std::vector<int> naturals(int max) {
    return naturals(0, max);
}

int to_octree_space(double x, double center, 
                    double half_width, int leaves) {
    int res = std::floor(((x - center) / (2 * half_width) + 0.5) * leaves);
    return res;
}

Octree::Octree(std::array<std::vector<double>,3>& p_elements,
               int depth):
    depth(depth),
    n_leaves_1d((int)pow(n_octants, depth)),
    elements(std::move(p_elements)),
    leaf_indices(elements[0].size())
{ 
    bounds = bounding_box(this->elements);
    bounds.half_width.loc[0] *= 1.01;

    for (unsigned int i = 0; i < elements[0].size(); ++i) {
        for (unsigned int d = 0; d < 3; d++) {
            leaf_indices[i][d] = to_octree_space(elements[d][i],
                                         bounds.center.loc[d],
                                         bounds.half_width.loc[d],
                                         n_leaves_1d);
        }
    }
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

