#ifndef __GEOM_H
#define __GEOM_H

#include <vector>
#include <algorithm>
/* A templated vector class. Can be used to represent a 2D or 3D
 * vector. The class is a thin wrapper over an array.
 */
template <int dim>
class Vec {
public:
    double loc[dim];

    inline bool operator==(const Vec<dim>& rhs) const {
        bool retval = true;
        for(int i = 0; i < dim; ++i) {
            retval &= (loc[i] == rhs.loc[i]);
        }
        return retval;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec<dim>& obj)
    {
        os << "{";
        for (int d = 0; d < dim - 1; ++d) {
            os << obj.loc[d] << ", ";
        }
        return os << obj.loc[dim - 1] << "}";
    }

    friend Vec<dim> operator+(Vec<dim> lhs, const Vec<dim>& rhs) {
        for(int i = 0; i < dim; ++i) {
            lhs.loc[i] += rhs.loc[i];
        }
        return lhs;
    }

    friend Vec<dim> operator-(Vec<dim> lhs, const Vec<dim>& rhs) {
        for(int i = 0; i < dim; ++i) {
            lhs.loc[i] -= rhs.loc[i];
        }
        return lhs;
    }

    friend Vec<dim> operator*(Vec<dim> lhs, const double& rhs) {
        for(int i = 0; i < dim; ++i) {
            lhs.loc[i] *= rhs;
        }
        return lhs;
    }

    friend Vec<dim> operator/(Vec<dim> lhs, const double& rhs) {
        for(int i = 0; i < dim; ++i) {
            lhs.loc[i] /= rhs;
        }
        return lhs;
    }
};


/* A box class defined by its center and half_width. These are using as
 * bounding boxes for the nodes in the octree hierarchy.
 */
template <int dim>
class Box {
public:
    Vec<dim> center;
    Vec<dim> half_width;

    friend std::ostream& operator<<(std::ostream& os, const Box<dim>& obj)
    {
        os << "{Box center=" << obj.center << 
               ", half_width=" << obj.half_width << "}";
        return os;
    }
};

template <int dim>
Box<dim> bounding_box(const std::vector<Vec<dim>> &elements)
{
    Vec<dim> min_corner = elements[0];
    Vec<dim> max_corner = elements[0];
    for (int i = 1; i < elements.size(); ++i) {
        for (int d = 0; d < dim; ++d) {
            min_corner.loc[d] = std::min(min_corner.loc[d],
                                         elements[i].loc[d]);
            max_corner.loc[d] = std::max(max_corner.loc[d],
                                         elements[i].loc[d]);
        }
    }

    Box<dim> ext;
    ext.center = (max_corner + min_corner) / 2.0;
    ext.half_width = (max_corner - min_corner) / 2.0;
    return ext;
}

/* Compute the bounding box of a set of vectors.
 */
template <int dim>
Box<dim> bounding_box_subset(const std::vector<Vec<dim>> &elements,
                             const std::vector<int>& indices) 
{
    Vec<dim> min_corner = elements[indices[0]];
    Vec<dim> max_corner = elements[indices[0]];
    for (int i = 1; i < indices.size(); ++i) {
        for (int d = 0; d < dim; ++d) {
            min_corner.loc[d] = std::min(min_corner.loc[d],
                                         elements[indices[i]].loc[d]);
            max_corner.loc[d] = std::max(max_corner.loc[d],
                                         elements[indices[i]].loc[d]);
        }
    }

    Box<dim> ext;
    ext.center = (max_corner + min_corner) / 2.0;
    ext.half_width = (max_corner - min_corner) / 2.0;
    return ext;
}

inline std::vector<int> naturals(int min, int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), min);
    return indices;
}

/* equivalent to range(0, max) in python */
inline std::vector<int> naturals(int max) {
    return naturals(0, max);
}

#endif
