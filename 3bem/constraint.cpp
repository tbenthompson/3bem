#include <iostream>
#include <algorithm>
#include <cassert>
#include "constraint.h"
#include "mesh.h"

namespace tbem {

std::ostream& operator<<(std::ostream& os, const Constraint& c) {
    os << "Constraint[[(RHS, " << c.rhs_value << "), ";
    for (auto v: c.dof_constraints) {
        os << "(" << v.first << ", " << v.second << "), ";
    }
    os << "]]";
    return os;
}

ConstraintMatrix ConstraintMatrix::add_constraints(
                                 const std::vector<Constraint>& constraints) {
    std::vector<std::pair<int, Constraint>> entries;
    for (auto it = c_map.begin(); it != c_map.end(); ++it) {
        entries.push_back({it->first, it->second});
    }

    for (std::size_t i = 0; i < constraints.size(); ++i) {
        auto in_constraint = constraints[i];
        auto new_dof_constraints = in_constraint.dof_constraints;
        auto last = std::max_element(new_dof_constraints.begin(), 
                                     new_dof_constraints.end(),
            [] (const DOFWeight<double>& a, const DOFWeight<double>& b) {
                return a.first < b.first; 
            });
        auto last_dof = last->first;
        std::iter_swap(last, new_dof_constraints.end() - 1);
        entries.push_back({last_dof, {new_dof_constraints, in_constraint.rhs_value}});
    }
    MapT new_map(entries.begin(), entries.end());
    return {new_map};
}

bool ConstraintMatrix::is_constrained(int dof) const {
    auto dof_constraint = c_map.find(dof);
    if (dof_constraint == c_map.end()) {
        return false;
    }
    return true;
}

ConstraintMatrix ConstraintMatrix::from_constraints(
        const std::vector<Constraint>& constraints) {
    ConstraintMatrix c;
    return c.add_constraints(constraints);
};
    
std::ostream& operator<<(std::ostream& os, const ConstraintMatrix& cm) {
    os << "ConstraintMatrix[[";
    for (auto it = cm.c_map.begin(); it != cm.c_map.end(); ++it) {
         os << "(" << it->first << ", " << it->second << "), ";
    }
    os << "]]";
    return os;
}
Constraint continuity_constraint(int dof1, int dof2) {
    return {
        {DOFWeight<double>{dof1, 1.0}, DOFWeight<double>{dof2, -1.0}},
        0.0
    };
}

Constraint offset_constraint(int dof1, int dof2, double offset) {
    return {
        continuity_constraint(dof1, dof2).dof_constraints,
        -offset
    };
}

Constraint boundary_condition(int dof, double value) {
    return {
        {DOFWeight<double>{dof, 1.0}},
        value
    };
}

template <int dim>
std::vector<Constraint> mesh_continuity(const Mesh<dim>& m, double eps) {

    std::vector<Constraint> constraints;
    for (std::size_t i = 0; i < m.facets.size(); i++) {
        for (std::size_t d1 = 0; d1 < dim; d1++) {
            auto i_pt = m.facets[i].vertices[d1];
            for (std::size_t j = i + 1; j < m.facets.size(); j++) {
                for (std::size_t d2 = 0; d2 < dim; d2++) {
                    auto j_pt = m.facets[j].vertices[d2];
                    if (!all(fabs(i_pt - j_pt) < eps)) {
                        continue;
                    } 
                    constraints.push_back(continuity_constraint(dim * i + d1,
                                                                dim * j + d2));
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
    auto out_c_map = c_matrix.c_map;
    for (std::size_t i = 0; i < disc.facets.size(); i++) {
        for (std::size_t d1 = 0; d1 < dim; d1++) {
            auto disc_pt = disc.facets[i].vertices[d1];
            for (std::size_t j = 0; j < surface.facets.size(); j++) {
                for (std::size_t d2 = 0; d2 < dim; d2++) {
                    auto surf_pt = surface.facets[j].vertices[d2];

                    // If the vertices do not overlap, nothing is done.
                    if (!all(fabs(disc_pt - surf_pt) < eps)) {
                        continue;
                    } 

                    // Is this DOF constrained? If not, move on.
                    int dof = dim * j + d2;
                    auto dof_and_constraint = out_c_map.find(dof);
                    if (dof_and_constraint == out_c_map.end()) {
                        continue;
                    }

                    // Get the other dof for the constraint.
                    int other_dof = dof_and_constraint->second.dof_constraints[0].first;
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
                    out_c_map.erase(dof);
                }
            }
        }
    }
    return ConstraintMatrix{out_c_map};
}

// INSTANTIATE TEMPLATES:
template 
std::vector<Constraint> mesh_continuity<2>(const Mesh<2>& m, double eps);
template 
std::vector<Constraint> mesh_continuity<3>(const Mesh<3>& m, double eps);

template
ConstraintMatrix apply_discontinuities<2>(const Mesh<2>& surface,
                                          const Mesh<2>& disc,
                                          const ConstraintMatrix& c_mat, double eps);
template
ConstraintMatrix apply_discontinuities<3>(const Mesh<3>& surface,
                                          const Mesh<3>& disc,
                                          const ConstraintMatrix& c_mat, double eps);

} // END namespace tbem
