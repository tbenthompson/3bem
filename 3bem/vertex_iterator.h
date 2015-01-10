#ifndef __QWEQASDKJVNS_VERTEX_ITERATOR_H
#define __QWEQASDKJVNS_VERTEX_ITERATOR_H
#include <cassert>

#include "mesh.h"

namespace tbem {

/* A barebones iterator for the vertices of a MeshField 
 * If this becomes heavily used, it might be worth becoming a more standard 
 * iterator design patterns in c++
 * http://stackoverflow.com/questions/8054273/how-to-implement-an-stl-style-iterator-and-avoid-common-pitfalls
 */
template <typename T, size_t dim>
struct FacetCornerIterator {
    typedef FacetCornerIterator<T,dim> iterator;

    const MeshField<T,dim>& mesh;
    size_t facet_idx;    
    size_t vertex_idx;

    FacetCornerIterator(const MeshField<T,dim>& mesh, int facet_idx, int vertex_idx):
        mesh(mesh),
        facet_idx(facet_idx),
        vertex_idx(vertex_idx)
    {} 

    friend std::ostream& operator<<(std::ostream& os, const iterator& it) {
        os << "(" << it.facet_idx << ", " << it.vertex_idx << ")";
        return os;
    }

    iterator& operator++() {
        *this += 1;
        return *this;
    }

    iterator& operator+=(size_t step) {
        int vertex_idx_too_far = vertex_idx + step; 
        facet_idx += vertex_idx_too_far / dim;
        vertex_idx = vertex_idx_too_far % dim;
        assert(vertex_idx < dim);
        return *this;
    }

    friend iterator operator+(const iterator& it, size_t step) {
        auto out_it = it;
        out_it += step;
        return out_it;
    }

    bool is_end() const {
        return *this == mesh.end();
    }

    int absolute_index() const {
        return facet_idx * dim + vertex_idx; 
    }

    const T& operator*() const {
        return get_vertex();
    }

    const FacetField<T,dim>& get_facet() const {
        return mesh.facets[facet_idx];
    }

    const T& get_vertex() const {
        return mesh.facets[facet_idx].vertices[vertex_idx];
    }
    
    friend bool operator==(const iterator& a, const iterator& b) {
        return (a.absolute_index() == b.absolute_index()) && 
            ((&a.mesh) == (&b.mesh));
    }

    friend bool operator!=(const iterator& a, const iterator& b) {
        return !(a == b);
    }
};

template <size_t dim>
using VertexIterator = FacetCornerIterator<Vec<double,dim>,dim>;

/* This is useful to put VertexIterator as the key in an unordered_map
 */
template <size_t dim>
struct HashVertexIterator {
    size_t operator()(const VertexIterator<dim>& it) const {
        return std::hash<size_t>()(it.absolute_index());
    }
};

} // END namespace tbem
#endif
