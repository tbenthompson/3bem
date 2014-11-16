#include "constraint.h"
#include <iostream>
#include <algorithm>

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
    for (std::size_t i = 0; i < total_dofs; i++) {
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

std::vector<double> ConstraintMatrix::get_unconstrained(const std::vector<double>& all) {
    std::vector<double> out;
    for (std::size_t i = 0; i < all.size(); i++) {
        auto dof_and_constraint = c_map.find(i);
        if (dof_and_constraint == c_map.end()) {
            out.push_back(all[i]);
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
