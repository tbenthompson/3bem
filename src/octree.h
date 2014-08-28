#ifndef __OCTREE_H
#define __OCTREE_H

#include <algorithm>
#include <vector>
#include <array>
#include <iostream>
#include <memory>

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

/* Compute the bounding box of a set of vectors.
 */
template <int dim>
Box<dim> bounding_box(const std::vector<Vec<dim> > &elements,
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

/* equivalent to range(0, max) in python */
std::vector<int> naturals(int max) {
    std::vector<int> indices(max);
    std::iota(indices.begin(), indices.end(), 0);
    return indices;
}

/* A shared static class to define some compile time information about
 * octrees.
 */
template <int dim>
class OctreeInfo {
public:
    static const int n_octants = (int)pow(2, dim);
};

template <int dim>
class OctreeNode {
public:
    OctreeNode(std::vector<Vec<dim> >& elements,
               int leaf_elements,
               const std::vector<int> &indices);
    Box<dim> extents;
    std::vector<OctreeNode<dim> > children;
    std::vector<int> indices;

    friend std::ostream& operator<<(std::ostream& os,
                                    const OctreeNode<dim>& obj)
    {
        os << "{OctreeNode extents=" << obj.extents << ", children={";
        for (auto c: obj.children) {
            os << c; 
        }
        os << "}, indices={";
        for (auto i: obj.indices) {
            os << i; 
        }
        os << "}}";
        return os;
    }
};

template <int dim>
OctreeNode<dim>::OctreeNode(std::vector<Vec<dim> >& elements,
                            int leaf_elements,
                            const std::vector<int> &indices)
{
    if (indices.size() == 0) {
        return;
    }

    extents = bounding_box(elements, indices);

    if (indices.size() <= leaf_elements) {
        this->indices = indices;
        return;
    }

    const int n_oct = OctreeInfo<dim>::n_octants;
    children.reserve(n_oct);
    const auto octant_indices = partition(elements, indices, extents.center);
    for(int i = 0; i < n_oct; ++i) {
        children.push_back(OctreeNode<dim>(elements, leaf_elements, octant_indices[i]));
    }
}

template <int dim>
std::array<std::vector<int>, OctreeInfo<dim>::n_octants> partition(
                                    const std::vector<Vec<dim> > &elements,
                                    const std::vector<int>& indices,
                                    const Vec<dim>& center) {
    const int n_oct = OctreeInfo<dim>::n_octants;
    std::vector<int> octant(indices.size());
    std::array<int, n_oct> counts;
    for (int o = 0; o < n_oct; o++) {
        counts[o] = 0;
    }

    for (int i = 0; i < indices.size(); i++) {
        octant[i] = find_octant(elements[indices[i]], center);
        counts[octant[i]] += 1;
    }
    std::array<std::vector<int>, n_oct> octants;
    for (int o = 0; o < n_oct; o++) {
        octants[o].resize(counts[o]);
        counts[o] = 0;
    }

    for (int i = 0; i < indices.size(); i++) {
        const int which_octant = octant[i];
        octants[which_octant][counts[which_octant]] = indices[i];
        counts[which_octant] += 1;
    }
    return octants;
}


template <int dim>
int find_octant(const Vec<dim>& point, const Vec<dim>& center) {
    int which_octant = 0;
    for(int d = 0; d < dim; ++d) {
        bool sign = point.loc[d] < center.loc[d];
        which_octant |= sign * (1 << d);
    }
    return which_octant;
}


template <int dim> 
class Octree {
public:
    Octree(std::unique_ptr<std::vector<Vec<dim> > > &elements,
           int leaf_elements);
    
    const int leaf_elements;
    const std::unique_ptr<std::vector<Vec<dim> > > elements;
    const OctreeNode<dim> root;
};

template <int dim>
Octree<dim>::Octree(std::unique_ptr<std::vector<Vec<dim> > > &elements,
                    int leaf_elements): 
    leaf_elements(leaf_elements),
    elements(std::move(elements)),
    root(*this->elements, leaf_elements, naturals(this->elements->size()))
{}

#endif
