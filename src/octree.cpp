#include "octree.h"
#include <assert.h>

Box bounding_box(const std::array<std::vector<double>,3>& x)
{
    Vec3 min_corner = {x[0][0], x[1][0], x[2][0]};
    Vec3 max_corner = {x[0][0], x[1][0], x[2][0]};
    for (unsigned int i = 1; i < x[0].size(); ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            min_corner.loc[d] = std::min(min_corner.loc[d], x[d][i]);
            max_corner.loc[d] = std::max(max_corner.loc[d], x[d][i]);
        }
    }

    Box ext;
    /* std::cout << max_corner << " " << min_corner << std::endl; */
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

Level::Level(const std::array<std::vector<double>,3>& p_elements,
             int depth,
             Box& bounds):
    depth(depth),
    n_cells_1d((int)pow(2, depth)),
    indices(p_elements[0].size()),
    buckets((int)pow(8, depth))
{
    init_indices(p_elements, bounds);
}

void Level::init_indices(const std::array<std::vector<double>,3>& p_elements,
                         Box& bounds) 
{
    for (unsigned int i = 0; i < p_elements[0].size(); ++i) {
        for (unsigned int d = 0; d < 3; d++) {
            indices[i][d] = to_octree_space(p_elements[d][i],
                                            bounds.center.loc[d],
                                            bounds.half_width.loc[d],
                                            n_cells_1d);
            // std::cout << p_elements[d][i] << " " << bounds.center.loc[d] <<
            //              " " << bounds.half_width.loc[d] << " " << n_cells_1d <<
            //              " " << indices[i][d] << std::endl;
            // assert(indices[i][d] >= 0);
            // assert(indices[i][d] < n_cells_1d);
        }
        /* std::cout << "Indices: " << indices[i][0] << 
         *           " " << indices[i][1] << " " << indices[i][2] << std::endl; */
        int loc = threed_to_1d(indices[i][0], indices[i][1], indices[i][2], n_cells_1d);
        /* std::cout << depth << " " << buckets.size() << " " << loc << 
         *              " " << (int)pow(8, depth) << std::endl; */
        // std::cout << buckets[0][0] << std::endl;
        buckets[loc].push_back(i);
    }
}

Octree::Octree(std::array<std::vector<double>,3>& p_elements,
               int max_depth):
    max_depth(max_depth),
    elements(std::move(p_elements))
{ 
    // std::cout << "Building octree on " << elements[0].size() << " points." << std::endl;
    bounds = bounding_box(this->elements);
    // Fudge the width to be slightly wider than necessary in order to make
    // sure that no points on the bounding box boundary.
    bounds.half_width.loc[0] *= 1.001;
    bounds.half_width.loc[1] *= 1.001;
    bounds.half_width.loc[2] *= 1.001;

    std::vector<uint64_t> morton_codes(this->elements[0].size());
    for (unsigned int i = 0; i < this->elements[0].size(); i++) {
        unsigned int x = to_octree_space(this->elements[0][i], 
                                         bounds.center.loc[0],
                                         bounds.half_width.loc[0],
                                         4);
        unsigned int y = to_octree_space(this->elements[1][i], 
                                         bounds.center.loc[1],
                                         bounds.half_width.loc[1],
                                         4);
        unsigned int z = to_octree_space(this->elements[2][i], 
                                         bounds.center.loc[2],
                                         bounds.half_width.loc[2],
                                         4);
        morton_codes[i] = morton_encode(x, y, z);
    }
    uint64_t center = morton_encode(0, 0, 2);
    std::cout << "Center = " << center << std::endl;
    int n_less = 0;
    for (unsigned int i = 0; i < this->elements[0].size(); i++) {
        if (morton_codes[i] < center) {
            n_less++;
        }
    }
    std::cout << n_less << std::endl;
    std::sort(morton_codes.begin(), morton_codes.end());
    std::cout << morton_codes[0] << std::endl;
// 
//     for (int d = 0; d < max_depth; d++) {
//         levels.emplace_back(this->elements, d, bounds);
//     }
}
