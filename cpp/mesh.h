#ifndef TBEMGGggGGggTTFDSSDf_MESH_H
#define TBEMGGggGGggTTFDSSDf_MESH_H
#include <vector>
#include "vec.h"

namespace tbem {

template <size_t dim>
struct VertexIterator;

template <size_t dim>
using Facet = Vec<Vec<double,dim>,dim>;

template <size_t dim>
struct Mesh {
    const std::vector<Facet<dim>> facets;

    Mesh(const std::vector<Facet<dim>>& facets):
        facets(facets)
    {}

    const Vec<double,dim>& get_vertex(size_t facet_idx, size_t vertex_idx) const;
    const Vec<double,dim>& get_vertex_from_dof(size_t absolute_index) const;
    size_t n_facets() const;
    size_t n_dofs() const;

    VertexIterator<dim> begin() const;
    VertexIterator<dim> end() const;

    Mesh<dim> refine(const std::vector<size_t>& refine_these) const;
    Mesh<dim> refine_once() const;
    Mesh<dim> refine_repeatedly(size_t times) const;

    Mesh<dim> remove_facets(const std::vector<size_t>& remove_these) const;

    static Mesh<dim> create_union(const std::vector<Mesh<dim>>& others);

    static Mesh<dim> from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
        const std::vector<std::array<size_t,dim>>& facets_by_vert_idx);
};

} //END NAMESPACE tbem

#endif
