#ifndef __AAABBBEEEEDDD_CONSTRAINT_H
#define __AAABBBEEEEDDD_CONSTRAINT_H
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "numbers.h"

namespace tbem {

/* This module handles imposing general linear constraints on an iterative
 * solution method.
 *
 * Constraints are used to ensure continuity between displacement values on
 * adjacent elements and to enforce boundary conditions. Constraints can also
 * be used to remove a null space from a problem. If the only boundary conditions 
 * are on the forces, then rigid body motions can superimposed resulting in an
 * infinite number of possible solutions. Imposing the constraint that average
 * displacement is zero is one of the simplest ways of solving this problem.
 */

template <typename T>
struct GeneralLinearTerm {
    const int dof;
    const T weight;

    bool operator==(const GeneralLinearTerm<T>& other) const {
        return dof == other.dof && weight == other.weight;
    }

    friend std::ostream& operator<<(std::ostream& os, const GeneralLinearTerm<T>& lt) {
        os << "(" << lt.dof << ", " << lt.weight << ")";
        return os;
    }
};

typedef GeneralLinearTerm<double> LinearTerm;

/* A constraint is composed of LinearTerms and a right hand side constant.
 * This structure represents one constraint equation, with each variable
 * times constant * term included in the "terms" vector and the RHS constant
 * in "rhs". So, for a constraint:
 * 3*x_0 - x_1 + 4*x_2 = 13.7
 * the equivalent ConstraintEQ would be 
 * ConstraintEQ eqtn{{LinearTerm{0,3}, LinearTerm{1,-1}, LinearTerm{2,4}}, 13.7};
 */
struct ConstraintEQ {
    const std::vector<LinearTerm> terms;
    const double rhs;

    friend std::ostream& operator<<(std::ostream& os, const ConstraintEQ& c) {
        os << "ConstraintEQ[[";
        os << "(rhs, " << c.rhs << "), ";
        for (auto t: c.terms) {
            os << t << ", ";
        }
        os << "]]";
        return os;
    }
};

ConstraintEQ boundary_condition(int dof, double value) {
    return {{LinearTerm{dof, 1.0}}, value};
}

ConstraintEQ continuity_constraint(int dof1, int dof2) {
    return {
        {LinearTerm{dof1, 1.0}, LinearTerm{dof2, -1.0}},
        0.0
    };
}

int find_last_dof_index(const ConstraintEQ& c);

ConstraintEQ filter_zero_terms(const ConstraintEQ& c, double eps = 1e-15);

/* This structure represents a constraint that has been rearranged so that
 * the constrained dof is alone on the left hand side.
 * A constraint:
 * 3*x_0 - x_1 + 4*x_2 = 13.7
 * with the constrained_dof = 2, would be rearranged to:
 * x_2 = (1/4)(13.7 + x_1 - 3*x_0)
 */
struct RearrangedConstraintEQ {
    // Q: Why aren't these variables marked const? 
    // A: Constructing a standard library map object 
    // (std::map or std::unordered_map) can be done either by adding in
    // items one by one or building a list of pairs and then inserting those
    // all at once. An immutable type can only be inserted as as list
    // of pairs because the map preallocates (and preconstructs) all of its
    // space and an immutable type has no assignment operator (you can't change
    // it!). But, the algorithm for ensuring the constraints are lower triangular
    // requires building the map item by item, so immutability is difficult for
    // this data type. Since this is only used internally, it's not too big
    // of a deal. Just don't change its member!
    int constrained_dof;
    std::vector<LinearTerm> terms;
    double rhs;

    friend std::ostream& operator<<(std::ostream& os, 
                                    const RearrangedConstraintEQ& c) {
        os << "RearrangedConstraintEQ[[";
        os << "(constrained_dof=" << c.constrained_dof << ", 1), ";
        os << "(rhs, " << c.rhs << "), ";
        for (auto t: c.terms) {
            os << t << ", ";
        }
        os << "]]";
        return os;
    }
};

RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
                                                  int constrained_index);

ConstraintEQ substitute(const ConstraintEQ& c, int constrained_dof_index,
                        const RearrangedConstraintEQ& subs_in);

typedef std::unordered_map<int,RearrangedConstraintEQ> ConstraintMapT;

bool is_constrained(const ConstraintMapT& dof_constraint_map, int dof);

RearrangedConstraintEQ make_lower_triangular(const ConstraintEQ& c,
                                             const ConstraintMapT& map);

struct ConstraintMatrix {
    const ConstraintMapT map;

    static ConstraintMatrix from_constraints(
            const std::vector<ConstraintEQ>& constraints);

    /* Accepts a reduced DOF vector and returns the full DOF vector. */
    template <typename T>
    std::vector<T> get_all(const std::vector<T>& in, int total_dofs) const;

    /* Accepts a full DOF vector and returns the reduced DOF vector.
     */
    template <typename T>
    std::vector<T> get_reduced(const std::vector<T>& all) const;

    template <typename T>
    std::vector<GeneralLinearTerm<T>>
    add_term_with_constraints(const GeneralLinearTerm<T>& entry) const;
};

template <typename T>
std::vector<T> ConstraintMatrix::get_all(const std::vector<T>& in, int total_dofs) const {
    std::vector<T> out(total_dofs); 

    int next_reduced_dof = 0;

    for (int dof_index = 0; dof_index < total_dofs; dof_index++) {
        if(is_constrained(map, dof_index)) {
            continue;
        }
        out[dof_index] = in[next_reduced_dof];
        next_reduced_dof++;
    }

    for (int dof_index = 0; dof_index < total_dofs; dof_index++) {
        if(!is_constrained(map, dof_index)) {
            continue;
        }
        auto constraint = map.find(dof_index)->second;
        auto val = constant<T>::make(constraint.rhs);
        for (size_t j = 0; j < constraint.terms.size(); j++) {
            val += constraint.terms[j].weight * out[constraint.terms[j].dof];
        }
        out[dof_index] = val;
    }
    return out;
}

template <typename T>
std::vector<GeneralLinearTerm<T>>
ConstraintMatrix::add_term_with_constraints(const GeneralLinearTerm<T>& entry) const {
    if (!is_constrained(map, entry.dof)) {
        return {entry};
    }

    const auto& constraint = map.find(entry.dof)->second;
    const auto& terms = constraint.terms;
    std::vector<GeneralLinearTerm<T>> out_terms;
    for (size_t i = 0; i < terms.size(); i++) {
        T recurse_weight = terms[i].weight * entry.weight;
        GeneralLinearTerm<T> new_entry{terms[i].dof, recurse_weight};
        const auto& terms = add_term_with_constraints(new_entry);
        for (const auto& t: terms) {
            out_terms.push_back(std::move(t));
        }
    }
    return out_terms;
}

template <typename T>
std::vector<T> ConstraintMatrix::get_reduced(const std::vector<T>& all) const {
    std::vector<T> condensed_dofs(all.size(), zeros<T>::make());
    for (size_t dof_idx = 0; dof_idx < all.size(); dof_idx++) {
        GeneralLinearTerm<T> term_to_add{(int)dof_idx, all[dof_idx]};
        auto expanded_term = add_term_with_constraints(term_to_add);
        for (const auto& t: expanded_term) {
            condensed_dofs[t.dof] += t.weight;
        }
    }

    std::vector<T> out;
    for (size_t dof_idx = 0; dof_idx < all.size(); dof_idx++) {
        if (is_constrained(map, dof_idx)) {
            continue;
        }
        out.push_back(condensed_dofs[dof_idx]);
    }

    return out;
}

} // END namespace tbem

#endif
