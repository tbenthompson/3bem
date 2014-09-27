#ifndef _OCTREE2_H
#define _OCTREE2_H

#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>
#include <limits>


/* A box class defined by its center and half_width. These are using as
 * bounding boxes for the nodes in the octree hierarchy.
 */
class Box {
public:
    std::array<double,3> center;
    std::array<double,3> half_width;
    std::array<double,3> min_corner;
    std::array<double,3> max_corner;

    friend std::ostream& operator<<(std::ostream& os, const Box& obj)
    {
        os << "{Box center={";
        os << obj.center[0] << ",";
        os << obj.center[1] << ",";
        os << obj.center[2] << "}";
        os << ", half_width={";
        os << obj.half_width[0] << ",";
        os << obj.half_width[1] << ",";
        os << obj.half_width[2] << "}}";
        return os;
    }
};

Box bounding_box(const std::vector<std::array<double,3>>& x);

Box get_child_box(std::array<int,3> which, Box& parent_bounds);

inline int to_octree_space(double x, double center, 
                           double half_width, int leaves) {
    int res = std::floor(((x - center) / (2 * half_width) + 0.5) * leaves);
    return res;
}


class OctreeCell {
public:
    // Vertical depth of the cell. The root is level 0.
    unsigned int level;

    Box bounds;

    // The location of the cell within this level. For example,
    // a level 2 cell would have a location ranging from 0 to 3
    // in each index, for a total of 8^2 = 64 cells on level 2.
    // TODO: May not be necessary to store this.
    std::array<int,3> loc;

    // The beginning and end indices of the elements in this cell.
    // The index range is [begin, end] (inclusive of end indices);
    unsigned int begin;
    unsigned int end;

    uint64_t min_code;
    uint64_t max_code;

    // The cell index of the children of this cell.
    // This could be replaced by ensuring a specific ordering of cells.
    // For example, the next cells until the level goes back up are all
    // children of this cell.
    // cell index of -1 means no child.
    std::array<int,8> children;
};

/* One quirk to the behavior of this octree implementation. 
 * All points must be on the interior of the octree, they cannot be on
 * the boundaries. This allows the "to_octree_space" function to ignore the
 * edge cases involving the boundaries.
 */
class Octree {
public:
    Octree(std::vector<std::array<double,3>>& elements,
           unsigned int max_elements_per_cell);

    uint64_t compute_morton(std::array<double,3> pt, const Box& bounds);
    void sort_elements();
    void build_children(OctreeCell& cell);
    void build_child(OctreeCell& cur_cell, int i, int j, int k);

    const unsigned int max_elements_per_cell;
    // The elements and z-order morton location codes of the points in the
    // octree. These could be made into a simple POD class.
    std::vector<std::array<double,3>> elements;
    std::vector<uint64_t> morton_codes;

    std::vector<OctreeCell> cells;
    
    //The deepest possible depth
    static const int deepest = 1 << 20;
};
#endif
