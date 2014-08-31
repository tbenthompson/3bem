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


int to_octree_space(double x, double center, double half_width, int leaves);

/* One quirk to the behavior of this octree implementation. 
 * All points must be on the interior of the octree, they cannot be on
 * the boundaries. This allows the "to_octree_space" function to ignore the
 * edge cases involving the boundaries.
 */
class Octree {
public:
    Octree(std::array<std::vector<double>,3>& elements,
           int depth);
    
    const static int n_octants = 8;
    const int depth;
    const int n_leaves_1d;
    const std::array<std::vector<double>, 3> elements;
    std::vector<std::array<int, 3>> leaf_indices;
    Box bounds;
};
#endif
