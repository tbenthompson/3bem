#ifndef __ASKDJAWERWFJS_CONTINUITY_BUILDER
#define __ASKDJAWERWFJS_CONTINUITY_BUILDER

#include "constraint.h"
#include "mesh.h"
#include "vertex_iterator.h"
    
namespace tbem {

// template <int dim>
// std::vector<ConstraintEQ> mesh_continuity(const Mesh<dim>& m,
//                                           double eps = 1e-10) {
// 
//     std::vector<ConstraintEQ> constraints;
//     for (auto outer_it = m.begin(); outer_it != m.end(); ++outer_it) {
//         auto outer_pt = *outer_it;
//         int outer_dof = dim * outer_it.facet_idx + outer_it.vertex_idx;
//         for (auto inner_it = m.begin(); inner_it != m.end(); ++inner_it) {
//             auto inner_pt = *inner_it;
//             int inner_dof = dim * inner_it.facet_idx + inner_it.vertex_idx;
//             if (!all(fabs(inner_pt - outer_pt) < eps)) {
//                 continue;
//             } 
//             constraints.push_back(continuity_constraint(outer_dof, inner_dof));
//         }
//     }
//     return constraints;
// }

template <int dim>
std::vector<ConstraintEQ> mesh_continuity(const Mesh<dim>& m,
                                          double eps = 1e-10) {

    std::vector<ConstraintEQ> constraints;
    for (std::size_t i = 0; i < m.facets.size(); i++) {
        for (std::size_t vertex1 = 0; vertex1 < dim; vertex1++) {
            auto i_pt = m.facets[i].vertices[vertex1];
            for (std::size_t j = i + 1; j < m.facets.size(); j++) {
                for (std::size_t vertex2 = 0; vertex2 < dim; vertex2++) {
                    auto j_pt = m.facets[j].vertices[vertex2];
                    if (!all(fabs(i_pt - j_pt) < eps)) {
                        continue;
                    } 
                    constraints.push_back(continuity_constraint(dim * i + vertex1,
                                                                dim * j + vertex2));
                }
            }
        }
    }
    return constraints;
}


template <int dim> 
ConstraintMatrix apply_discontinuities(const Mesh<dim>& surface,
                                       const Mesh<dim>& disc,
                                       const ConstraintMatrix& c_matrix,
                                       double eps = 1e-10) {
    auto out_map = c_matrix.map;
    for (std::size_t i = 0; i < disc.facets.size(); i++) {
        for (std::size_t vertex1 = 0; vertex1 < dim; vertex1++) {
            auto disc_pt = disc.facets[i].vertices[vertex1];
            for (std::size_t j = 0; j < surface.facets.size(); j++) {
                for (std::size_t vertex2 = 0; vertex2 < dim; vertex2++) {
                    auto surf_pt = surface.facets[j].vertices[vertex2];

                    // If the vertices do not overlap, nothing is done.
                    if (!all(fabs(disc_pt - surf_pt) < eps)) {
                        continue;
                    } 

                    // Is this DOF constrained? If not, move on.
                    int dof = dim * j + vertex2;
                    if (!is_constrained(out_map, dof)) {
                        continue;
                    }

                    // Get the other dof for the constraint.
                    int other_dof = out_map.find(dof)->second.terms[0].dof;
                    int other_vert = other_dof % 3;
                    int other_face = (other_dof - other_vert) / 3;
                    
                    // Calculate which side of the disc face the other vertex is on.
                    auto my_side = which_side_facet<dim>(disc.facets[i].vertices, 
                                                surface.facets[j].vertices);
                    auto other_side = which_side_facet<dim>(disc.facets[i].vertices, 
                                                surface.facets[other_face].vertices);
                    if (my_side == other_side) {
                        continue;
                    }
                    out_map.erase(dof);
                }
            }
        }
    }
    return ConstraintMatrix{out_map};
}

} //END namespace tbem

#endif
