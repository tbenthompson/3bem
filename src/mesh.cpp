#include "mesh.h"
#include "util.h"

/* Produces 2 new segments by splitting the current segment 
 * in half.
 */
std::array<Facet<2>,2> refine_facet(const Facet<2>& f) {
    auto midpt = (f.vertices[0] + f.vertices[1]) / 2.0;
    return {
        Facet<2>{{f.vertices[0], midpt}},
        Facet<2>{{midpt, f.vertices[1]}}
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
std::array<Facet<3>,4> refine_facet(const Facet<3>& f) {
    auto midpt01 = (f.vertices[0] + f.vertices[1]) / 2.0;
    auto midpt12 = (f.vertices[1] + f.vertices[2]) / 2.0;
    auto midpt20 = (f.vertices[2] + f.vertices[0]) / 2.0;

    //Maintain the orientation. Since vertex 1 is "next" after vertex 0 
    //in the original triangle, midpt01 should be "next" after vertex 0
    //in the refined triangle. Following this principle for all the triangles
    //gives this set of new faces
    return {
        Facet<3>{{f.vertices[0], midpt01, midpt20}},
        Facet<3>{{f.vertices[1], midpt12, midpt01}},
        Facet<3>{{f.vertices[2], midpt20, midpt12}},
        Facet<3>{{midpt01, midpt12, midpt20}}
    };
}

/* Applies the Mesh modification function to a refined facet. This
 * can be used to, for example, enforce that all refined vertices in
 * a sphere mesh lie a certain distance from the sphere's center.
 */
template <int dim>
Facet<dim> refine_modify(const Facet<dim>& f, const Mesh<dim>& m) {
    if (!m.has_refine_mod) {
        return f;
    }
    Vec<Vec<double,dim>,dim> vertices;
    for (int d = 0; d < dim; d++) {
        vertices[d] = m.refine_mod(f.vertices[d]);
    }
    return Facet<dim>{vertices};
}

template 
Facet<2> refine_modify(const Facet<2>& f, const Mesh<2>& m);
template
Facet<3> refine_modify(const Facet<3>& f, const Mesh<3>& m);

template <int dim>
Mesh<dim> 
Mesh<dim>::refine(const std::vector<int>& refine_these) const {
    std::vector<Facet<dim>> out_facets;

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
                auto mod_r = refine_modify(r, *this);
                out_facets.push_back(mod_r);
            }
            current++;
        } else {
            out_facets.push_back(facets[i]);
        }
    }
    return Mesh<dim>{out_facets, has_refine_mod, refine_mod};
}

/* A helper function to refine all the facets. */
template <int dim>
Mesh<dim> 
Mesh<dim>::refine() const {
    return refine(naturals(facets.size()));
}

/* A helper function to refine all the facets multiple times. */
template <int dim>
Mesh<dim> 
Mesh<dim>::refine_repeatedly(unsigned int times) const {
    if (times == 0) {
        return *this;
    }
    return refine_repeatedly(times - 1).refine();
}

/* Given a list of vertices and a list of faces that references the vertex
 * list, this will construct a mesh object.
 */
template <int dim>
Mesh<dim>
Mesh<dim>::from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod,
                         const typename Mesh<dim>::RefineFnc& refine_mod) {
    std::vector<Facet<dim>> new_facets;
    for (auto in_facet: facets) { 
        Vec<Vec<double,dim>,dim> out_verts;
        for (int d = 0; d < dim; d++) {
            out_verts[d] = vertices[in_facet[d]];
        }
        auto out_facet = Facet<dim>{out_verts};
        new_facets.push_back(out_facet);
    }
    return Mesh<dim>{new_facets, has_refine_mod, refine_mod};
}

template class Mesh<2>;
template class Mesh<3>;
