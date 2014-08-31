#ifndef __OCTREE_H
#define __OCTREE_H

#include <algorithm>
#include <vector>
#include <array>
#include <iostream>
#include <memory>
#include <omp.h>
#include "geom.h"

/* A shared static class to define some compile time information about
 * octrees.
 */
template <int dim>
class OctreeInfo {
public:
    static const int n_octants = (int)pow(2, dim);
};

template <typename T, int dim>
class OctreeNode {
public:
    OctreeNode(std::vector<T>& elements,
               int leaf_elements,
               const std::vector<int> &indices);

    Box<dim> extents;
    std::vector<OctreeNode<T, dim>> children;
    std::vector<int> indices;

    friend std::ostream& operator<<(std::ostream& os,
                                    const OctreeNode<T, dim>& obj)
    {
        os << "{OctreeNode extents=" << obj.extents << ", children={";
        for (auto c: obj.children) {
            os << c << ","; 
        }
        os << "}, indices={";
        for (auto i: obj.indices) {
            os << i << ","; 
        }
        os << "}}";
        return os;
    }
};

template <typename T, int dim>
OctreeNode<T, dim>::OctreeNode(std::vector<T>& elements,
                            int leaf_elements,
                            const std::vector<int> &indices)
{
    if (indices.size() == 0) {
        return;
    }

    extents = bounding_box_subset(elements, indices);

    if (indices.size() <= leaf_elements) {
        this->indices = indices;
        return;
    }

    const int n_oct = OctreeInfo<dim>::n_octants;
    children.reserve(n_oct);
    const auto octant_indices = partition(elements, indices, extents.center);

    for(int i = 0; i < n_oct; ++i) {
        children.push_back(OctreeNode<T, dim>(elements, leaf_elements, octant_indices[i]));
    }
}

template <typename T, int dim>
std::array<std::vector<int>, OctreeInfo<dim>::n_octants> partition(
                                    const std::vector<T> &elements,
                                    const std::vector<int>& indices,
                                    const Vec<dim>& center) {
    // Partitioning happens in two steps. First, the octant for each element
    // is calculated. Then, afterwards, the elements are grouped into an array
    // where each list contains all of that octants indices.
    const int n_oct = OctreeInfo<dim>::n_octants;
    std::vector<int> octant(indices.size());
    std::array<int, n_oct> counts;
    for (int o = 0; o < n_oct; ++o) {
        counts[o] = 0;
    }

    for (int i = 0; i < indices.size(); ++i) {
        octant[i] = find_octant(elements[indices[i]], center);
        counts[octant[i]] += 1;
    }
    std::array<std::vector<int>, n_oct> octants;
    for (int o = 0; o < n_oct; ++o) {
        octants[o].resize(counts[o]);
        counts[o] = 0;
    }

    for (int i = 0; i < indices.size(); ++i) {
        const int which_octant = octant[i];
        octants[which_octant][counts[which_octant]] = indices[i];
        counts[which_octant] += 1;
    }
    return octants;
}


template <typename T, int dim>
int find_octant(const T& point, const Vec<dim>& center) {
    int which_octant = 0;
    for(int d = 0; d < dim; ++d) {
        bool sign = point.loc[d] < center.loc[d];
        which_octant |= sign * (1 << d);
    }
    return which_octant;
}


template <typename T, int dim> 
class Octree {
public:
    Octree(std::unique_ptr<std::vector<T>> &elements,
           int leaf_elements);
    
    const int leaf_elements;
    const std::unique_ptr<std::vector<T>> elements;
    const OctreeNode<T, dim> root;
};

template <typename T, int dim>
Octree<T, dim>::Octree(std::unique_ptr<std::vector<T>> &elements,
                    int leaf_elements): 
    leaf_elements(leaf_elements),
    elements(std::move(elements)),
    root(*this->elements, leaf_elements, naturals(this->elements->size()))
{}

#endif
