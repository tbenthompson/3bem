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

Box bounding_box(const std::array<std::vector<double>,3>& x);

inline int to_octree_space(double x, double center, 
                    double half_width, int leaves) {
    int res = std::floor(((x - center) / (2 * half_width) + 0.5) * leaves);
    return res;
}
//
// method to seperate bits from a given integer 3 positions apart
inline uint64_t split_by_3(unsigned int a){
    uint64_t x = a & 0x1fffff; // we only look at the first 21 bits
    // shift left 32 bits, OR with self, and 
    // 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 32) & 0x1f00000000ffff;  
    // shift left 32 bits, OR with self, and 
    // 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 16) & 0x1f0000ff0000ff;  
    // shift left 32 bits, OR with self, and
    // 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 8) & 0x100f00f00f00f00f;
    // shift left 32 bits, OR with self, and 
    // 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; 
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

inline uint64_t morton_encode(unsigned int x, unsigned int y, unsigned int z){
    uint64_t answer = 0 | split_by_3(x) | split_by_3(y) << 1 | split_by_3(z) << 2;
    return answer;
}

/* One quirk to the behavior of this octree implementation. 
 * All points must be on the interior of the octree, they cannot be on
 * the boundaries. This allows the "to_octree_space" function to ignore the
 * edge cases involving the boundaries.
 */
class Octree {
public:
    Octree(std::array<std::vector<double>,3>& elements,
           int max_depth);
    
    const int max_depth;
    const std::array<std::vector<double>, 3> elements;

    Box bounds;
};
#endif
