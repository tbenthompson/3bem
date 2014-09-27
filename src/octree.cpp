#include "octree.h"
#include "algs.h"
#include <assert.h>

Box bounding_box(const std::vector<std::array<double,3>>& x)
{
    Box ext;
    ext.min_corner = {x[0][0], x[0][1], x[0][2]};
    ext.max_corner = {x[0][0], x[0][1], x[0][2]};
    for (unsigned int i = 1; i < x.size(); ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            ext.min_corner[d] = std::min(ext.min_corner[d], x[i][d]);
            ext.max_corner[d] = std::max(ext.max_corner[d], x[i][d]);
        }
    }


    for (int d = 0; d < 3; ++d) {
        ext.center[d] = (ext.max_corner[d] + ext.min_corner[d]) / 2.0;
        // Fudge the width to be slightly wider than necessary in order to make
        // sure that no points on the bounding box boundary.
        ext.half_width[d] = (ext.max_corner[d] - ext.min_corner[d]) / 2.0;
    }
    /* std::cout << max_corner << " " << min_corner << std::endl; */
    return ext;
}

Box get_child_box(std::array<int,3> which, Box& parent_bounds) {
    Box child_box;
    for (int d = 0; d < 3; d++) {
        child_box.min_corner[d] = parent_bounds.min_corner[d] + 
                               which[d] * parent_bounds.half_width[d];
        child_box.max_corner[d] = parent_bounds.center[d] + 
                               which[d] * parent_bounds.half_width[d];
        child_box.half_width[d] = parent_bounds.half_width[d] / 2.0;
        child_box.center[d] += (child_box.min_corner[d] + 
                                child_box.max_corner[d]) / 2.0;
    }
    return child_box;
}

Octree::Octree(std::vector<std::array<double,3>>& p_elements,
               unsigned int max_elements_per_cell):
    max_elements_per_cell(max_elements_per_cell),
    elements(std::move(p_elements)),
    morton_codes(elements.size())
{ 
    Box root_bounds = bounding_box(elements);
    for (int d = 0; d < 3; d++) {
        root_bounds.half_width[d] *= 1.001;
    }
    for (unsigned int i = 0; i < elements.size(); i++) {
        morton_codes[i] = compute_morton(elements[i], root_bounds);
    }
    sort_elements();

    OctreeCell root = {0, root_bounds, {0,0,0}, 
                       0, (unsigned int)morton_codes.size() - 1,
                       compute_morton(root_bounds.min_corner, root_bounds),
                       compute_morton(root_bounds.max_corner, root_bounds)};
    cells.push_back(root);
    build_children(root);
}

uint64_t Octree::compute_morton(std::array<double,3> pt,
                                const Box& bounds) {
    unsigned int x = to_octree_space(pt[0], bounds.center[0],
                                     bounds.half_width[0], deepest);
    unsigned int y = to_octree_space(pt[1], bounds.center[1],
                                     bounds.half_width[1], deepest);
    unsigned int z = to_octree_space(pt[2], bounds.center[2],
                                     bounds.half_width[2], deepest);
    return morton_encode(z, y, x);
}

void Octree::sort_elements() {
    //replace with a parallel merge sort
    auto p = sort_permutation(morton_codes, [](uint64_t a, uint64_t b) {return a < b;});
    morton_codes = apply_permutation(morton_codes, p);
    elements = apply_permutation(elements, p);
}

void Octree::build_child(OctreeCell& cur_cell, int i, int j, int k) {
    int child_idx = 4 * i + 2 * j + k;

    unsigned int level = cur_cell.level + 1;

    Box child_box = get_child_box({i, j, k}, cur_cell.bounds);

    auto loc = cur_cell.loc;
    loc[0] = loc[0] * 2 + i;
    loc[1] = loc[1] * 2 + j;
    loc[2] = loc[2] * 2 + k;

    //Compute morton code of each corner.
    auto morton_steps = ((cur_cell.max_code - cur_cell.min_code) / 8);
    auto min_code = cur_cell.min_code + child_idx * morton_steps;
    auto max_code = cur_cell.max_code - (7 - child_idx) * morton_steps;
    std::cout << min_code << " " << max_code << std::endl;
    assert(min_code < max_code);

    unsigned int low = std::upper_bound(morton_codes.begin() + cur_cell.begin,
                                morton_codes.begin() + cur_cell.end, min_code)
             - morton_codes.begin();
    unsigned int up = std::lower_bound(morton_codes.begin() + cur_cell.begin,
                               morton_codes.begin() + cur_cell.end, max_code) 
             - morton_codes.begin();
    std::cout << up - low << std::endl;
    std::cout << low << " " << up << std::endl;

    if (up - low == 0) {
        std::cout << "EMPTY CELL" << std::endl; 
    }
    

    unsigned int new_cell_idx = cells.size();
    cells.push_back({level, child_box, loc, low, up});
    
    std::cout << child_idx << " " << cells.size() << std::endl;
    //build children recursively -- depth first
    if (up - low > max_elements_per_cell) {
        build_children(cells[new_cell_idx]);
    }
    cur_cell.children[child_idx] = new_cell_idx;
}

void Octree::build_children(OctreeCell& cur_cell) {
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            for(int k = 0; k < 2; k++) {
                build_child(cur_cell, i, j, k);
            }
        }
    }
}
