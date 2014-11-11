#include "mesh.h"
#include "util.h"

std::unordered_map<int, int> find_duplicate_map(const Mesh3D& mesh, double eps) {
    std::unordered_map<int, int> old_to_new;
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        //TODO: FIX THE O(N^2) problem -- blocking, refactor octree to use Vec3
        for (unsigned int j = i + 1; j < mesh.vertices.size(); j++) {
            if (!same_vertex(mesh.vertices[i], mesh.vertices[j], eps)) {
                continue;
            } 
            // Map the deleted vertex to the label of the retained vertex.
            // The retained vertex will always be the last in the set of 
            // identical vertices. So, if indices 1,3,9,11 are all the same
            // point, only vertex 11 will be kept.
            old_to_new[i] = j;
        }
        auto found = old_to_new.find(i);
        if (found == old_to_new.end()) {
            old_to_new[i] = i;
        }
    }

    int next_new_vertex = 0;
    // Traverse the vertices in reverse, so that when we get to a deleted
    // vertex it has already been assigned its new location.
    for (int i = (int)mesh.vertices.size() - 1; i >= 0; i--) {
        if (old_to_new[i] != (int)i) {
            old_to_new[i] = old_to_new[old_to_new[i]];
        } else {
            old_to_new[i] = next_new_vertex;
            next_new_vertex++;
        }
    }

    return old_to_new;
}

int count_unique_vertices(const std::unordered_map<int,int>& old_to_new) {
    int max_index = 0;
    for (auto it = old_to_new.begin(); it != old_to_new.end(); it++) {
        max_index = std::max(max_index, it->second);
    }
    return max_index + 1;
}

std::vector<Vec3<double>> unique_vertices(
        const Mesh3D& mesh, std::unordered_map<int,int>& old_to_new) {
    int n_unique = count_unique_vertices(old_to_new);
    std::vector<Vec3<double>> vertices(n_unique);
    for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        vertices[old_to_new[i]] = mesh.vertices[i]; 
    }
    return vertices;
}
            
Mesh3D clean_mesh(const Mesh3D& mesh, double vertex_smear) {
    //Find duplicate vertices
    auto old_to_new = find_duplicate_map(mesh, vertex_smear);
    auto new_vertices = unique_vertices(mesh, old_to_new);
    // std::cout << mesh.vertices.size() << " " << new_vertices.size() << " " << 
        // count_unique_vertices(old_to_new) << std::endl;

    std::vector<std::array<int,3>> new_faces;
    for (unsigned int i = 0; i < mesh.faces.size(); i++) {
        std::array<int,3> face;
        for (int v = 0; v < 3; v++) {
            face[v] = old_to_new[mesh.faces[i][v]];
        }
        new_faces.push_back(face);
    }

    Mesh3D retval{new_vertices, new_faces, mesh.has_refine_mod, mesh.refine_mod};
    return retval;
}

void refine_face(Mesh3D& new_mesh, std::array<int, 3> face) {
    // Find the new vertex and faces.
    const auto v0 = new_mesh.vertices[face[0]];
    const auto v1 = new_mesh.vertices[face[1]];
    const auto v2 = new_mesh.vertices[face[2]];

    // Calculate the midpoints of each edge of the triangle. These
    // are used to refine the triangle
    auto midpt01 = (v0 + v1) / 2.0;
    auto midpt12 = (v1 + v2) / 2.0;
    auto midpt20 = (v2 + v0) / 2.0;
    if (new_mesh.has_refine_mod) {
        midpt01 = new_mesh.refine_mod(midpt01);
        midpt12 = new_mesh.refine_mod(midpt12);
        midpt20 = new_mesh.refine_mod(midpt20);
    }

    // Add the vertices while grabbing their indices.
    int midpt01_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt01);
    int midpt12_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt12);
    int midpt20_idx = new_mesh.vertices.size();
    new_mesh.vertices.push_back(midpt20);

    //Maintain the orientation. Since vertex 1 is "next" after vertex 0 
    //in the original triangle, midpt01 should be "next" after vertex 0
    //in the refined triangle. Following this principle for all the triangles
    //gives this set of new faces
    new_mesh.faces.push_back({face[0], midpt01_idx, midpt20_idx});
    new_mesh.faces.push_back({face[1], midpt12_idx, midpt01_idx});
    new_mesh.faces.push_back({face[2], midpt20_idx, midpt12_idx});
    new_mesh.faces.push_back({midpt01_idx, midpt12_idx, midpt20_idx});
}

Mesh3D refine_mesh(const Mesh3D& m, std::vector<int> refine_these) {
    Mesh3D new_mesh {m.vertices, {}, m.has_refine_mod, m.refine_mod};
    new_mesh.vertices = m.vertices;

    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::sort(refine_these.begin(), refine_these.end());

    // The next refined edge.
    int current = 0;

    for (int i = 0; i < (int)m.faces.size(); i++) {
        if (i == refine_these[current]) {
            current += 1;
            refine_face(new_mesh, m.faces[i]);
        } else {
            new_mesh.faces.push_back(m.faces[i]);
        }
    }

    return new_mesh;
}

Mesh3D refine_mesh(const Mesh3D& m) {
    return refine_mesh(m, naturals(m.faces.size()));
}

Mesh3D refine_clean(const Mesh3D& m, unsigned int times) {
    if (times == 0) {
        return clean_mesh(m);
    }
    Mesh3D ret = refine_mesh(m);
    for(unsigned int i = 1; i < times; i++) {
        ret = refine_mesh(ret);
    }
    ret = clean_mesh(ret);
    return ret;
}

// Explicitly instantiate the possible mesh templates.
template class Mesh<2>;
template class Mesh<3>;

std::array<Facet<2>,2> refine_facet(const Facet<2>& f) {
    auto midpt = (f.vertices[0] + f.vertices[1]) / 2.0;
    return {
        Facet<2>{{f.vertices[0], midpt}},
        Facet<2>{{midpt, f.vertices[1]}}
    };
}

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

template <int dim>
Facet<dim> refine_modify(const Facet<dim>& f, const NewMesh<dim>& m) {
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
Facet<2> refine_modify(const Facet<2>& f, const NewMesh<2>& m);
template
Facet<3> refine_modify(const Facet<3>& f, const NewMesh<3>& m);

template <int dim>
NewMesh<dim> 
NewMesh<dim>::refine(const std::vector<int>& refine_these) const {
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
    return NewMesh<dim>{out_facets, has_refine_mod, refine_mod};
}

template <int dim>
NewMesh<dim> 
NewMesh<dim>::refine() const {
    return refine(naturals(facets.size()));
}

template <int dim>
NewMesh<dim> 
NewMesh<dim>::refine_repeatedly(unsigned int times) const {
    if (times == 0) {
        return *this;
    }
    return refine_repeatedly(times - 1).refine();
}

template <int dim>
NewMesh<dim>
NewMesh<dim>::from_vertices_faces(const std::vector<Vec<double,dim>>& vertices,
                         const std::vector<std::array<int,dim>>& facets,
                         bool has_refine_mod,
                         const typename NewMesh<dim>::RefineFnc& refine_mod) {
    std::vector<Facet<dim>> new_facets;
    for (auto in_facet: facets) { 
        Vec<Vec<double,dim>,dim> out_verts;
        for (int d = 0; d < dim; d++) {
            out_verts[d] = vertices[in_facet[d]];
        }
        auto out_facet = Facet<dim>{out_verts};
        new_facets.push_back(out_facet);
    }
    return NewMesh<dim>{new_facets, has_refine_mod, refine_mod};
}

template class NewMesh<2>;
template class NewMesh<3>;
