#include "mesh.h"
#include "vertex_iterator.h"
#include "vec_ops.h"
#include <algorithm>

namespace tbem {

template <size_t dim>
const Vec<double,dim>& Mesh<dim>::get_vertex(size_t facet_idx, size_t vertex_idx) const
{
    return facets[facet_idx][vertex_idx];
}

template <size_t dim>
const Vec<double,dim>& Mesh<dim>::get_vertex_from_dof(size_t absolute_index) const
{
    auto v_idx = absolute_index % dim;
    auto f_idx = (absolute_index - v_idx) / dim;
    return get_vertex(f_idx, v_idx);
}

template <size_t dim>
size_t Mesh<dim>::n_facets() const 
{
    return facets.size();
}

template <size_t dim>
size_t Mesh<dim>::n_dofs() const 
{
    return dim * n_facets();
}

template <size_t dim>
VertexIterator<dim> Mesh<dim>::begin() const 
{
    return VertexIterator<dim>(*this, 0, 0);
}

template <size_t dim>
VertexIterator<dim> Mesh<dim>::end() const 
{
    return VertexIterator<dim>(*this, facets.size(), 0);
}

/* Produces 2 new segments by splitting the current segment 
 * in half. 
 */
std::array<Facet<2>,2> refine_facet(const Facet<2>& f) 
{
    auto midpt = (f[0] + f[1]) / 2.0;
    return {{
        {{f[0], midpt}},
        {{midpt, f[1]}}
    }};
}

/* Produces 4 new triangles from an initial triangle by adding a
 * vertex to each edge of the triangle
 *         /\          /\
 *        /  \   ->   /__\
 *       /    \      /\  /\
 *      /______\    /__\/__\
 *      (My first ever ASCII art, a masterpiece that will
 *      stand the test of time!)
 *      Approximate date of creation: October 2014
 */
std::array<Facet<3>,4> refine_facet(const Facet<3>& f) 
{
    auto midpt01 = (f[0] + f[1]) / 2.0;
    auto midpt12 = (f[1] + f[2]) / 2.0;
    auto midpt20 = (f[2] + f[0]) / 2.0;

    //Maintain the orientation. Since vertex 1 is "next" after vertex 0 
    //in the original triangle, midpt01 should be "next" after vertex 0
    //in the refined triangle. Following this principle for all the triangles
    //gives this set of new faces
    return {{
        {{f[0], midpt01, midpt20}},
        {{f[1], midpt12, midpt01}},
        {{f[2], midpt20, midpt12}},
        {{midpt01, midpt12, midpt20}}
    }};
}

template <size_t dim>
Mesh<dim> 
Mesh<dim>::refine(const std::vector<size_t>& refine_these) const 
{
    if (refine_these.empty()) {
        return *this;
    }

    std::vector<Facet<dim>> out_facets;

    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    auto sorted_refines = refine_these;
    std::sort(sorted_refines.begin(), sorted_refines.end());

    // The next index of sorted_refines.
    size_t current = 0;

    for (size_t i = 0; i < facets.size(); i++) {
        if (current == sorted_refines.size()) {
            out_facets.push_back(facets[i]);
            continue;
        }
        if (i == sorted_refines[current]) {
            auto refined = refine_facet(facets[i]);
            for (auto r: refined) {
                out_facets.push_back(r);
            }
            current++;
        } else {
            out_facets.push_back(facets[i]);
        }
    }
    return Mesh<dim>{out_facets};
}

/* A helper function to refine all the facets. */
template <size_t dim>
Mesh<dim> 
Mesh<dim>::refine_once() const 
{
    return refine(range(facets.size()));
}

/* A helper function to refine all the facets multiple times. */
template <size_t dim>
Mesh<dim> 
Mesh<dim>::refine_repeatedly(size_t times) const 
{
    if (times == 0) {
        return *this;
    }
    return refine_repeatedly(times - 1).refine_once();
}

template <size_t dim>
Mesh<dim> 
Mesh<dim>::create_union(const std::vector<Mesh<dim>>& meshes) 
{

    std::vector<Facet<dim>> new_facets;
    for (size_t i = 0; i < meshes.size(); i++) {
        for (auto f: meshes[i].facets) {
            new_facets.push_back(f);
        }
    }

    return Mesh<dim>{new_facets};
}

/* Given a list of vertices and a list of faces that references the vertex
 * list, this will construct a mesh object.
 */
template <size_t dim>
Mesh<dim>
Mesh<dim>::from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<size_t,dim>>& facets_by_vert_idx) 
{
    std::vector<Facet<dim>> facets;
    for (auto in_facet: facets_by_vert_idx) { 
        Facet<dim> out_verts;
        for (size_t d = 0; d < dim; d++) {
            out_verts[d] = vertices[in_facet[d]];
        }
        facets.push_back(out_verts);
    }
    return Mesh<dim>{facets};
}

template struct Mesh<2>;
template struct Mesh<3>;

} //END NAMESPACE tbem
