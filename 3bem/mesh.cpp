#include "mesh.h"
#include "util.h"


namespace tbem {

/* Produces 2 new segments by splitting the current segment 
 * in half. 
 */
template <typename T>
std::array<FacetField<T,2>,2> refine_facet(const FacetField<T,2>& f) {
    auto midpt = (f.vertices[0] + f.vertices[1]) / 2.0;
    return {
        FacetField<T,2>{{f.vertices[0], midpt}},
        FacetField<T,2>{{midpt, f.vertices[1]}}
    };
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
std::array<FacetField<T,3>,4> refine_facet(const FacetField<T,3>& f) {
    auto midpt01 = (f.vertices[0] + f.vertices[1]) / 2.0;
    auto midpt12 = (f.vertices[1] + f.vertices[2]) / 2.0;
    auto midpt20 = (f.vertices[2] + f.vertices[0]) / 2.0;

    //Maintain the orientation. Since vertex 1 is "next" after vertex 0 
    //in the original triangle, midpt01 should be "next" after vertex 0
    //in the refined triangle. Following this principle for all the triangles
    //gives this set of new faces
    return {
        FacetField<T,3>{{f.vertices[0], midpt01, midpt20}},
        FacetField<T,3>{{f.vertices[1], midpt12, midpt01}},
        FacetField<T,3>{{f.vertices[2], midpt20, midpt12}},
        FacetField<T,3>{{midpt01, midpt12, midpt20}}
    };
}

/* Applies the Mesh modification function to a refined facet. This
 * can be used to, for example, enforce that all refined vertices in
 * a sphere mesh lie a certain distance from the sphere's center.
 */
template <typename T,int dim>
FacetField<T,dim> refine_modify(const FacetField<T,dim>& f,
                                const MeshField<T,dim>& m) {
    if (!m.has_refine_mod) {
        return f;
    }
    Vec<T,dim> vertices;
    for (int d = 0; d < dim; d++) {
        vertices[d] = m.refine_mod(f.vertices[d]);
    }
    return FacetField<T,dim>{vertices};
}

template <typename T, int dim>
MeshField<T,dim> 
MeshField<T,dim>::refine(const std::vector<int>& refine_these) const {
    std::vector<FacetField<T,dim>> out_facets;

    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::vector<int> sorted_refines = refine_these;
    std::sort(sorted_refines.begin(), sorted_refines.end());

    // The next index of sorted_refines.
    int current = 0;

    for (unsigned int i = 0; i < facets.size(); i++) {
        if (i == refine_these[current]) {
            auto refined = refine_facet(facets[i]);
            for (auto r: refined) {
                auto mod_r = refine_modify<T,dim>(r, *this);
                out_facets.push_back(mod_r);
            }
            current++;
        } else {
            out_facets.push_back(facets[i]);
        }
    }
    return MeshField<T,dim>{out_facets, has_refine_mod, refine_mod};
}

/* A helper function to refine all the facets. */
template <typename T, int dim>
MeshField<T,dim> 
MeshField<T,dim>::refine() const {
    return refine(naturals(facets.size()));
}

/* A helper function to refine all the facets multiple times. */
template <typename T, int dim>
MeshField<T,dim> 
MeshField<T,dim>::refine_repeatedly(unsigned int times) const {
    if (times == 0) {
        return *this;
    }
    return refine_repeatedly(times - 1).refine();
}

template <typename T, int dim>
MeshField<T,dim> 
MeshField<T,dim>::form_union(const std::vector<MeshField<T,dim>>& meshes) {

    std::vector<FacetField<T,dim>> new_facets;
    for (int i = 0; i < meshes.size(); i++) {
        if (meshes[i].has_refine_mod == true) {
            throw std::domain_error("Mesh unions can only be formed from meshes\
                                     with has_refine_mod == false");
        }
        for (auto f: meshes[i].facets) {
            new_facets.push_back(f);
        }
    }

    return MeshField<T,dim>{new_facets, false, nullptr};
}

/* Given a list of vertices and a list of faces that references the vertex
 * list, this will construct a mesh object.
 */
template <typename T, int dim>
MeshField<T,dim>
MeshField<T,dim>::from_vertices_faces(const std::vector<T>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod,
                         const typename MeshField<T,dim>::RefineFnc& refine_mod) {
    std::vector<FacetField<T,dim>> new_facets;
    for (auto in_facet: facets) { 
        Vec<T,dim> out_verts;
        for (int d = 0; d < dim; d++) {
            out_verts[d] = vertices[in_facet[d]];
        }
        auto out_facet = FacetField<T,dim>{out_verts};
        new_facets.push_back(out_facet);
    }
    return MeshField<T,dim>{new_facets, has_refine_mod, refine_mod};
}

template class MeshField<double,2>;
template class MeshField<Vec<double,2>,2>;
template class MeshField<double,3>;
template class MeshField<Vec<double,3>,3>;

} //END NAMESPACE tbem
