#ifndef _OCTREE2_H
#define _OCTREE2_H

#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>
#include <limits>

/* A templated vector class. Can be used to represent a 2D or 3D
 * vector. The class is a thin wrapper over an array.
 */
class Vec3 {
public:
    double loc[3];

    inline bool operator==(const Vec3& rhs) const {
        bool retval = true;
        for(int i = 0; i < 3; ++i) {
            retval &= (loc[i] == rhs.loc[i]);
        }
        return retval;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3& obj)
    {
        os << "{";
        for (int d = 0; d < 3 - 1; ++d) {
            os << obj.loc[d] << ", ";
        }
        return os << obj.loc[3 - 1] << "}";
    }

    friend Vec3 operator+(Vec3 lhs, const Vec3& rhs) {
        for(int i = 0; i < 3; ++i) {
            lhs.loc[i] += rhs.loc[i];
        }
        return lhs;
    }

    friend Vec3 operator-(Vec3 lhs, const Vec3& rhs) {
        for(int i = 0; i < 3; ++i) {
            lhs.loc[i] -= rhs.loc[i];
        }
        return lhs;
    }

    friend Vec3 operator*(Vec3 lhs, const double& rhs) {
        for(int i = 0; i < 3; ++i) {
            lhs.loc[i] *= rhs;
        }
        return lhs;
    }

    friend Vec3 operator/(Vec3 lhs, const double& rhs) {
        for(int i = 0; i < 3; ++i) {
            lhs.loc[i] /= rhs;
        }
        return lhs;
    }
};


/* A box class defined by its center and half_width. These are using as
 * bounding boxes for the nodes in the octree hierarchy.
 */
class Box {
public:
    Vec3 center;
    Vec3 half_width;

    friend std::ostream& operator<<(std::ostream& os, const Box& obj)
    {
        os << "{Box center=" << obj.center << 
               ", half_width=" << obj.half_width << "}";
        return os;
    }
};

Box bounding_box(const std::array<std::vector<double>,3>& x);
std::vector<int> naturals(int min, int max);
std::vector<int> naturals(int max);

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
