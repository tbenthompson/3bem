#ifndef __AAABBBEEEEDDD_CONSTRAINT_H
#define __AAABBBEEEEDDD_CONSTRAINT_H
#include <vector>
#include <map>
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

struct LinearTerm {
    const int dof;
    const double weight;

    bool operator==(const LinearTerm& other) const {
        return dof == other.dof && weight == other.weight;
    }

    friend std::ostream& operator<<(std::ostream& os, const LinearTerm& lt) {
        os << "(" << lt.dof << ", " << lt.weight << ")";
        return os;
    }
};

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

ConstraintEQ continuity_constraint(int dof1, int dof2) {
    return {
        {LinearTerm{dof1, 1.0}, LinearTerm{dof2, -1.0}},
        0.0
    };
}

/* This structure represents a constraint that has been rearranged so that
 * the constrained dof is alone on the left hand side.
 * A constraint:
 * 3*x_0 - x_1 + 4*x_2 = 13.7
 * with the constrained_dof = 2, would be rearranged to:
 * x_2 = (1/4)(13.7 + x_1 - 3*x_0)
 */
struct RearrangedConstraintEQ {
    // Q: Why aren't these variables marked const? 
    // A: It isn't possible to include this class as a member of a standard
    // library map (std::map or std::unordered_map) if there isn't a valid
    // assignment operator. A valid assignment operator requires mutable state.
    // Try to avoid taking advantage of this mutable state.
    // Since RearrangedConstraintEQ is only used internally, this isn't a big
    // deal.
    int constrained_dof;
    std::vector<LinearTerm> terms;
    double rhs;

    friend std::ostream& operator<<(std::ostream& os, const RearrangedConstraintEQ& c) {
        os << "RearrangedConstraintEQ[[";
        os << "(constrained_dof, " << c.constrained_dof << "), ";
        os << "(rhs, " << c.rhs << "), ";
        for (auto t: c.terms) {
            os << t << ", ";
        }
        os << "]]";
        return os;
    }
};

struct ConstraintMatrix {
    typedef std::map<int,RearrangedConstraintEQ> MapT;
    const MapT map;

    ConstraintMatrix add_constraints(const std::vector<ConstraintEQ>& constraints);
    static 
    ConstraintMatrix from_constraints(const std::vector<ConstraintEQ>& constraints);

    /* Accepts a reduced DOF vector and returns the full DOF vector. */
    template <typename T>
    std::vector<T> get_all(const std::vector<T>& in, int total_dofs) const;

    /* Accepts a full DOF vector and returns the reduced DOF vector.
     */
    template <typename T>
    std::vector<T> get_reduced(const std::vector<T>& all) const;
};

inline RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
                                                  int constrained_index) {
    assert(constrained_index < (int)c.terms.size());
    const auto& constrained_term = c.terms[constrained_index]; 

    std::vector<LinearTerm> divided_negated_terms;
    for (const auto& t: c.terms) {
        if (t == constrained_term) {
            continue;
        }
        LinearTerm divided_negated{t.dof, -t.weight / constrained_term.weight};
        divided_negated_terms.push_back(divided_negated);
    }

    double divided_rhs = c.rhs / constrained_term.weight;

    return RearrangedConstraintEQ{
        constrained_term.dof,
        divided_negated_terms,
        divided_rhs
    };
}

inline bool is_constrained(const ConstraintMatrix::MapT& dof_constraint_map,
                           int dof) {
    auto it = dof_constraint_map.find(dof);
    if (it == dof_constraint_map.end()) {
        return false;
    }
    return true;
}

inline int find_last_dof_index(const ConstraintEQ& c) {
    auto highest_dof_term = std::max_element(c.terms.begin(), c.terms.end(),
        [&] (const LinearTerm& a, const LinearTerm& b) {
            return a.dof < b.dof;
        });
    return std::distance(c.terms.begin(), highest_dof_term);
}

typedef std::vector<LinearTerm>::const_iterator LinearTermIterator;
LinearTermIterator find_term_with_dof(const std::vector<LinearTerm>& terms, int dof) {
    return std::find_if(terms.begin(), terms.end(),
        [&] (const LinearTerm& sub_lt) {
            return sub_lt.dof == dof;
        }
    );
}

template <typename T>
bool none_found(const typename T::const_iterator& it, const T& vec) {
    return it == std::end(vec);
}

template <typename T>
void remove(std::vector<T>& vec, size_t pos)
{
    typename std::vector<T>::iterator it = vec.begin();
    std::advance(it, pos);
    vec.erase(it);
}

ConstraintEQ substitute(const ConstraintEQ& c,
                        int constrained_dof_index,
                        const RearrangedConstraintEQ& subs_in) {
    const auto& constrained_term = c.terms[constrained_dof_index];
    assert(constrained_term.dof == subs_in.constrained_dof);
    double multiplicative_factor = constrained_term.weight;

    std::vector<LinearTerm> out_terms;
    auto which_subs_terms_unused = integers(subs_in.terms.size());
    for (size_t i = 0; i < c.terms.size(); i++) {
        if ((int)i == constrained_dof_index) {
            continue;
        }
        const auto& term = c.terms[i];
        const auto& subs_term = find_term_with_dof(subs_in.terms, term.dof);

        if (none_found(subs_term, subs_in.terms)) {
            out_terms.push_back(term); 
        } else {
            double additive_weight = subs_term->weight * multiplicative_factor;
            out_terms.push_back(LinearTerm{term.dof, term.weight + additive_weight});
            int subs_term_index = std::distance(subs_in.terms.begin(), subs_term);
            assert(subs_term_index < (int)subs_in.terms.size());
            remove(which_subs_terms_unused, subs_term_index);
        } 
    }

    for (int subs_term_index: which_subs_terms_unused) {
        const auto& subs_term = subs_in.terms[subs_term_index];
        double additive_weight = subs_term.weight * multiplicative_factor;
        out_terms.push_back(LinearTerm{subs_term.dof, additive_weight});
    }

    double out_rhs = c.rhs - multiplicative_factor * subs_in.rhs;
    return {out_terms, out_rhs};
}

RearrangedConstraintEQ make_lower_triangular(const ConstraintEQ& c,
                                             const ConstraintMatrix::MapT map) {
    int last_dof_index = find_last_dof_index(c);
    int last_dof = c.terms[last_dof_index].dof;
    if (is_constrained(map, last_dof)) {
        auto last_dof_constraint = map.find(last_dof)->second;
        ConstraintEQ c_subs = substitute(c, last_dof_index, last_dof_constraint);
        return make_lower_triangular(c_subs, map);
    }
    return isolate_term_on_lhs(c, last_dof_index);
}

ConstraintMatrix ConstraintMatrix::from_constraints(
        const std::vector<ConstraintEQ>& constraints) {
    ConstraintMatrix c;
    return c.add_constraints(constraints);
};

ConstraintMatrix ConstraintMatrix::add_constraints(
        const std::vector<ConstraintEQ>& constraints) {
    MapT new_map = map;

    for (size_t i = 0; i < constraints.size(); i++) {
        const auto& c = constraints[i];
        auto lower_tri_constraint = make_lower_triangular(c, new_map);
        new_map[lower_tri_constraint.constrained_dof] =
            std::move(lower_tri_constraint);
    }

    return ConstraintMatrix{new_map};
}

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
    

} // END namespace tbem

#endif
