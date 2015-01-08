#ifndef __QWEQASDKJVNS_VERTEX_ITERATOR_H
#define __QWEQASDKJVNS_VERTEX_ITERATOR_H

namespace tbem {

/* A barebones iterator for the vertices of a MeshField 
 * If this becomes heavily used, it would be worth conforming to standard 
 * iterator design patterns in c++
 * http://stackoverflow.com/questions/8054273/how-to-implement-an-stl-style-iterator-and-avoid-common-pitfalls
 */
template <typename T, int dim>
struct VertexIterator {
    const MeshField<T,dim>& mesh;
    size_t facet_idx;    
    size_t vertex_idx;

    VertexIterator(const MeshField<T,dim>& mesh, int facet_idx, int vertex_idx):
        mesh(mesh),
        facet_idx(facet_idx),
        vertex_idx(vertex_idx)
    {} 

    VertexIterator<T,dim>& operator++() {
        vertex_idx++;
        if (vertex_idx == dim) {
            vertex_idx = 0;
            facet_idx++;
        }
        return *this;
    }

    const T& operator*() const {
        return mesh.facets[facet_idx].vertices[vertex_idx];
    }

    friend bool operator==(const VertexIterator<T,dim>& a,
                           const VertexIterator<T,dim>& b) {
        return (a.facet_idx == b.facet_idx) && (a.vertex_idx == b.vertex_idx);
    }

    friend bool operator!=(const VertexIterator<T,dim>& a,
                           const VertexIterator<T,dim>& b) {
        return !(a == b);
    }
};

} // END namespace tbem
#endif
