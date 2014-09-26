#include "octree.h"
#include <assert.h>

Box bounding_box(const std::array<std::vector<double>,3>& x)
{
    std::array<double,3> min_corner = {x[0][0], x[1][0], x[2][0]};
    std::array<double,3> max_corner = {x[0][0], x[1][0], x[2][0]};
    for (unsigned int d = 0; d < 3; ++d) {
        for (unsigned int i = 1; i < x[0].size(); ++i) {
            min_corner[d] = std::min(min_corner[d], x[d][i]);
            max_corner[d] = std::max(max_corner[d], x[d][i]);
        }
    }


    Box ext;
    for (int d = 0; d < 3; ++d) {
        ext.center[d] = (max_corner[d] + min_corner[d]) / 2.0;
        // Fudge the width to be slightly wider than necessary in order to make
        // sure that no points on the bounding box boundary.
        ext.half_width[d] = 1.001 * (max_corner[d] - min_corner[d]) / 2.0;
    }
    /* std::cout << max_corner << " " << min_corner << std::endl; */
    return ext;
}

Octree::Octree(std::array<std::vector<double>,3>& p_elements,
               int max_depth):
    max_depth(max_depth),
    elements(std::move(p_elements))
{ 
    // std::cout << "Building octree on " << elements[0].size() 
    //           << " points." << std::endl;

    bounds = bounding_box(this->elements);

    //2^d cells at the leaf level.
    int cells_1d_leaf = 1 << max_depth;

    std::vector<uint64_t> morton_codes(this->elements[0].size());
    std::map<uint64_t, std::vector<
    //insert the morton codes into the balanced binary search tree in
    //std::map. sorting and level adaptivity both achieved!
    for (unsigned int i = 0; i < this->elements[0].size(); i++) {
        unsigned int x = to_octree_space(this->elements[0][i], 
                                         bounds.center[0],
                                         bounds.half_width[0],
                                         cells_1d_leaf);
        unsigned int y = to_octree_space(this->elements[1][i], 
                                         bounds.center[1],
                                         bounds.half_width[1],
                                         cells_1d_leaf);
        unsigned int z = to_octree_space(this->elements[2][i], 
                                         bounds.center[2],
                                         bounds.half_width[2],
                                         cells_1d_leaf);
        morton_codes[i] = morton_encode(x, y, z);
    }

    // std::cout << morton_codes[0] << std::endl;

    uint64_t center = morton_encode(0, 0, 2);
    // std::cout << "Center = " << center << std::endl;
    int n_less = 0;
    for (unsigned int i = 0; i < this->elements[0].size(); i++) {
        if (morton_codes[i] < center) {
            n_less++;
        }
    }
    // std::cout << n_less << std::endl;
// 
//     for (int d = 0; d < max_depth; d++) {
//         levels.emplace_back(this->elements, d, bounds);
//     }
}
