#include <iostream>
#include <algorithm>
#include <cassert>
#include "constraint.h"
#include "mesh.h"

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
    MapT new_map;
    for (auto it = c_map.begin(); it != c_map.end(); ++it) {
        new_map[it->first] = it->second;
    }

    for (std::size_t i = 0; i < constraints.size(); ++i) {
        auto in_constraint = constraints[i];
        auto last = std::max_element(in_constraint.dof_constraints.begin(), 
                                     in_constraint.dof_constraints.end(),
            [] (const DOFWeight& a, const DOFWeight& b) {
                return a.first < b.first; 
            });
        auto last_dof = last->first;
        //TODO: Make the constrained dof first, not last! simpler
        std::iter_swap(last, in_constraint.dof_constraints.end() - 1);
        new_map[last_dof] = in_constraint;
    }
    return {new_map};
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

std::vector<double> ConstraintMatrix::get_all(const std::vector<double>& in,
                                              int total_dofs) {
    std::vector<double> out(total_dofs); 
    int next_in = 0;
    for (int i = 0; i < total_dofs; i++) {
        auto dof_and_constraint = c_map.find(i);
        if (dof_and_constraint == c_map.end()) {
            out[i] = in[next_in];
            next_in++;
            continue;
        }

        auto dof_constraint = dof_and_constraint->second.dof_constraints;
        double out_val = dof_and_constraint->second.rhs_value;
        for (std::size_t j = 0; j < dof_constraint.size() - 1; j++) {
            out_val -= dof_constraint[j].second * out[dof_constraint[j].first];
        }
        auto this_dof_weight = (dof_constraint.end() - 1)->second;
        out[i] = out_val / this_dof_weight;
    }
    return out;
}

//TODO: Is there an alternate formulation for this that does not require modifying the
//input vec? Should recursive constraints be dealt with upfront? The constrained-dof 
//is last rule should prevent any cyclic constraints.
void ConstraintMatrix::add_vec_with_constraints(DOFWeight entry,
                                                std::vector<double>& vec) {
    auto dof_and_constraint = c_map.find(entry.first);
    if (dof_and_constraint == c_map.end()) {
        vec[entry.first] += entry.second;
        return;
    }

    auto constraint = dof_and_constraint->second;
    auto dof_weights = constraint.dof_constraints;
    int n_dof_weights = dof_weights.size();
    double this_weight = dof_weights[n_dof_weights - 1].second;
    for (int i = 0; i < n_dof_weights - 1; i++) {
        double recurse_weight = -dof_weights[i].second * entry.second / this_weight;
        DOFWeight new_entry{dof_weights[i].first, recurse_weight};
        add_vec_with_constraints(new_entry, vec);
    }
}

std::vector<double> ConstraintMatrix::get_reduced(const std::vector<double>& all) {
    std::vector<double> condensed_dofs(all.size(), 0.0);
    for (std::size_t i = 0; i < all.size(); i++) {
        add_vec_with_constraints(DOFWeight{i, all[i]}, condensed_dofs);
    }

    std::vector<double> out;
    for (std::size_t i = 0; i < all.size(); i++) {
        auto dof_and_constraint = c_map.find(i);
        if (dof_and_constraint == c_map.end()) {
            out.push_back(condensed_dofs[i]);
            continue;
        }
    }

    return out;
}

Constraint continuity_constraint(int dof1, int dof2) {
    return {
        {DOFWeight{dof1, 1.0}, DOFWeight{dof2, -1.0}},
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
        {DOFWeight{dof, 1.0}},
        value
    };
}

//TODO: FIX THE O(N^2) problem here, use hashes or octree?
template <int dim>
std::vector<Constraint> mesh_continuity(const Mesh<dim>& m, double eps) {

    std::vector<Constraint> constraints;
    for (unsigned int i = 0; i < m.facets.size(); i++) {
        for (unsigned int d1 = 0; d1 < dim; d1++) {
            for (unsigned int j = i + 1; j < m.facets.size(); j++) {
                for (unsigned int d2 = 0; d2 < dim; d2++) {
                    if (!all(fabs(m.facets[i].vertices[d1] - 
                                  m.facets[j].vertices[d2]) < eps)) {
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

template 
std::vector<Constraint> mesh_continuity<2>(const Mesh<2>& m, double eps);
template 
std::vector<Constraint> mesh_continuity<3>(const Mesh<3>& m, double eps);

