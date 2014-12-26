#include "octree.h"
#include "numerics.h"
#include "algs.h"
#include "vec.h"
#include <assert.h>

namespace tbem {

Box box_from_min_max(const Vec3<double> min_corner,
                     const Vec3<double> max_corner) {
    Box ext;
    ext.min_corner = min_corner;
    ext.max_corner = max_corner;
    for (int d = 0; d < 3; ++d) {
        ext.center[d] = (ext.max_corner[d] + ext.min_corner[d]) / 2.0;
        ext.half_width[d] = (ext.max_corner[d] - ext.min_corner[d]) / 2.0;
    }
    ext.radius2 = hypot2(ext.half_width);
    return ext;
}

Box box_from_center_half_width(const Vec3<double>& center,
                               const Vec3<double>& half_width) {
    Box ext;
    ext.center = center;
    ext.half_width = half_width;
    for (int d = 0; d < 3; d++) {
        ext.min_corner[d] = ext.center[d] - ext.half_width[d];
        ext.max_corner[d] = ext.center[d] + ext.half_width[d];
    }
    ext.radius2 = hypot2(ext.half_width);
    return ext;
}

Box bounding_box(const std::array<std::vector<double>,3>& x, int begin, int end)
{
    assert(begin < end);
    Vec3<double> min_corner = {x[0][begin], x[1][begin], x[2][begin]};
    Vec3<double> max_corner = {x[0][begin], x[1][begin], x[2][begin]};
    for (int i = begin + 1; i < end; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            min_corner[d] = std::min(min_corner[d], x[d][i]);
            max_corner[d] = std::max(max_corner[d], x[d][i]);
        }
    }
    return box_from_min_max(min_corner, max_corner);
}

Octree::Octree(std::array<std::vector<double>,3>& p_elements,
               unsigned int max_elements_per_cell):
    max_elements_per_cell(max_elements_per_cell),
    elements(std::move(p_elements)),
    morton_codes(elements[0].size())
{ 
    bounds = bounding_box(elements, 0, n_elements());

    // Fudge the width for the box used to calculate morton codes
    // so that no point is exactly on the outer boundary of the box.
    // This allows ignoring the edge case.
    morton_bounds = bounds;
    for (int d = 0; d < 3; d++) {
        morton_bounds.half_width[d] *= 1.001;
    }

    for (unsigned int i = 0; i < n_elements(); i++) {
        morton_codes[i] = compute_morton(elements[0][i],
                                         elements[1][i],
                                         elements[2][i]);
    }
    sort_elements();

    OctreeCell root = {0, {0,0,0}, bounds,
                       0, (unsigned int)morton_codes.size(),
                       compute_morton(bounds.min_corner),
                       compute_morton(bounds.max_corner),
                       {-1, -1, -1, -1, -1, -1, -1, -1},
                       false};
    if (n_elements() > max_elements_per_cell) {
        root.children = build_children(root);
    } else {
        root.is_leaf = true;
    }
    cells.push_back(root);
}

uint64_t Octree::compute_morton(Vec3<double> pt) {
    return compute_morton(pt[0], pt[1], pt[2]);
}

uint64_t Octree::compute_morton(double ptx, double pty, double ptz) {
    unsigned int x = to_octree_space(ptx, morton_bounds.center[0],
                                     morton_bounds.half_width[0], deepest);
    unsigned int y = to_octree_space(pty, morton_bounds.center[1],
                                     morton_bounds.half_width[1], deepest);
    unsigned int z = to_octree_space(ptz, morton_bounds.center[2],
                                     morton_bounds.half_width[2], deepest);
    return morton_encode(z, y, x);
}

void Octree::sort_elements() {
    permutation = sort_permutation(morton_codes, 
                              [](uint64_t a, uint64_t b) {return a < b;});
    morton_codes = apply_permutation(morton_codes, permutation);
    for (int d = 0; d < 3; d++) {
        elements[d] = apply_permutation(elements[d], permutation);
    }
}

int Octree::build_child(OctreeCell& cur_cell, int i, int j, int k) {
    int child_idx = 4 * i + 2 * j + k;


    //Compute morton code of each corner.
    auto morton_steps = ((cur_cell.max_code - cur_cell.min_code) / 8);
    auto min_code = cur_cell.min_code + child_idx * morton_steps;
    auto max_code = cur_cell.max_code - (7 - child_idx) * morton_steps;
    assert(max_code >= min_code);


    unsigned int begin = std::lower_bound(morton_codes.begin() + cur_cell.begin,
                                morton_codes.begin() + cur_cell.end, min_code)
             - morton_codes.begin();
    unsigned int end = std::upper_bound(morton_codes.begin() + cur_cell.begin,
                               morton_codes.begin() + cur_cell.end, max_code) 
             - morton_codes.begin();
    assert(end >= begin);
    // std::cout << cur_cell.begin << " " << cur_cell.end << " " << begin << " " << end << std::endl;

    // Every cell must have at least two points in it.
    if (end - begin == 0) {
        return -1;
    }

    Box box;
    if (end - begin == 1) {

        Vec3<double> center;
        Vec3<double> half_width;
        for (int d = 0; d < 3; d++) {
            center[d] = (elements[d][begin] + cur_cell.bounds.center[d]) / 2.0;
            half_width[d] = std::fabs(elements[d][begin] - 
                                      cur_cell.bounds.center[d]) / 2.5;
        }
        box = box_from_center_half_width(center, half_width);
    } else {
        box = bounding_box(elements, begin, end);
    }

    auto level = cur_cell.level + 1;

    auto loc = cur_cell.loc;
    loc[0] = loc[0] * 2 + i;
    loc[1] = loc[1] * 2 + j;
    loc[2] = loc[2] * 2 + k;

    OctreeCell new_cell = {level, loc, box, begin, end, min_code, max_code,
                          {-1, -1, -1, -1, -1, -1, -1, -1}, false};
    
    //build children recursively -- depth first
    if (end - begin > max_elements_per_cell) {
        new_cell.children = build_children(new_cell);
    } else {
        new_cell.is_leaf = true;
    }
    cells.push_back(new_cell);
    return cells.size() - 1;
}

std::array<int,8> Octree::build_children(OctreeCell& cur_cell) {
    std::array<int,8> child_indices;
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            for(int k = 0; k < 2; k++) {
                child_indices[4 * i + 2 * j + k] = build_child(cur_cell, i, j, k);
            }
        }
    }
    return child_indices;
}

OctreeCell& Octree::get_root() {
    return cells[get_root_index()];
}

int Octree::get_root_index() const {
    return cells.size() - 1;
}

} // END namespace tbem
