#include "mesh.h"
#include "vertex_iterator.h"
#include "vec_ops.h"
#include <algorithm>

namespace tbem {


template <typename T, size_t dim>
size_t Function<T,dim>::n_facets() const {
    return facets.size();
}

template <typename T, size_t dim>
size_t Function<T,dim>::n_dofs() const {
    return dim * n_facets();
}

template <typename T, size_t dim>
FunctionDOFIterator<T,dim> Function<T,dim>::begin() const {
    return FunctionDOFIterator<T,dim>(*this, 0, 0);
}

template <typename T, size_t dim>
FunctionDOFIterator<T,dim> Function<T,dim>::end() const {
    return FunctionDOFIterator<T,dim>(*this, facets.size(), 0);
}

/* Produces 2 new segments by splitting the current segment 
 * in half. 
 */
template <typename T>
std::array<Vec<T,2>,2> refine_facet(const Vec<T,2>& f) {
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
 */
template <typename T>
std::array<Vec<T,3>,4> refine_facet(const Vec<T,3>& f) {
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

template <typename T, size_t dim>
Function<T,dim> 
Function<T,dim>::refine(const std::vector<int>& refine_these) const {
    if (refine_these.empty()) {
        return *this;
    }

    std::vector<Vec<T,dim>> out_facets;

    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::vector<int> sorted_refines = refine_these;
    std::sort(sorted_refines.begin(), sorted_refines.end());

    // The next index of sorted_refines.
    int current = 0;

    for (unsigned int i = 0; i < facets.size(); i++) {
        if (i == refine_these[current]) {
            auto refined = refine_facet(facets[i]);
            for (const auto& r: refined) {
                out_facets.push_back(r);
            }
            current++;
        } else {
            out_facets.push_back(facets[i]);
        }
    }
    return Function<T,dim>{out_facets};
}

/* A helper function to refine all the facets. */
template <typename T, size_t dim>
Function<T,dim> 
Function<T,dim>::refine() const {
    return refine(range(facets.size()));
}

/* A helper function to refine all the facets multiple times. */
template <typename T, size_t dim>
Function<T,dim> 
Function<T,dim>::refine_repeatedly(unsigned int times) const {
    if (times == 0) {
        return *this;
    }
    return refine_repeatedly(times - 1).refine();
}

template <typename T, size_t dim>
Function<T,dim> 
Function<T,dim>::create_union(const std::vector<Function<T,dim>>& meshes) {

    std::vector<Vec<T,dim>> new_facets;
    for (int i = 0; i < meshes.size(); i++) {
        for (const auto& f: meshes[i].facets) {
            new_facets.push_back(f);
        }
    }

    return Function<T,dim>{new_facets};
}

/* Given a list of vertices and a list of faces that references the vertex
 * list, this will construct a mesh object.
 */
template <typename T, size_t dim>
Function<T,dim>
Function<T,dim>::from_vertices_faces(const std::vector<T>& vertices,
                         const std::vector<std::array<int,dim>>& facets_by_vert_idx) {
    std::vector<Vec<T,dim>> facets;
    for (const auto& in_facet: facets_by_vert_idx) { 
        Vec<T,dim> out_verts;
        for (int d = 0; d < dim; d++) {
            out_verts[d] = vertices[in_facet[d]];
        }
        facets.push_back(out_verts);
    }
    return Function<T,dim>{facets};
}

template class Function<double,2>;
template class Function<Vec<double,2>,2>;
template class Function<double,3>;
template class Function<Vec<double,3>,3>;

} //END NAMESPACE tbem
