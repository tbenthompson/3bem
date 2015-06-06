#ifndef TBEMQWEQASDKJVNS_VERTEX_ITERATOR_H
#define TBEMQWEQASDKJVNS_VERTEX_ITERATOR_H
#include <cassert>
#include <iostream>

#include "mesh.h"

namespace tbem {

/* A barebones iterator for the vertices of a Mesh 
 * If this becomes heavily used, it might be worth becoming a more standard 
 * iterator design patterns in c++
 * http://stackoverflow.com/questions/8054273/how-to-implement-an-stl-style-iterator-and-avoid-common-pitfalls
 */
template <size_t dim>
struct VertexIterator {
    typedef VertexIterator<dim> iterator;

    const Mesh<dim>& mesh;
    size_t facet_idx;    
    size_t vertex_idx;

    VertexIterator(const Mesh<dim>& mesh, int facet_idx, int vertex_idx):
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

    const Vec<double,dim>& operator*() const {
        return get_vertex();
    }

    const Facet<dim>& get_facet() const {
        return mesh.facets[facet_idx];
    }

    const Vec<double,dim>& get_vertex() const {
        return mesh.get_vertex(facet_idx, vertex_idx);
    }
    
    friend bool operator==(const iterator& a, const iterator& b) {
        return (a.absolute_index() == b.absolute_index()) && 
            ((&a.mesh) == (&b.mesh));
    }

    friend bool operator!=(const iterator& a, const iterator& b) {
        return !(a == b);
    }
};

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
