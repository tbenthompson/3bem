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

struct LinearTerm {
    const size_t dof;
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

inline ConstraintEQ boundary_condition(size_t dof, double value) {
    return {{LinearTerm{dof, 1.0}}, value};
}

inline ConstraintEQ continuity_constraint(size_t dof1, size_t dof2) {
    return {
        {LinearTerm{dof1, 1.0}, LinearTerm{dof2, -1.0}},
        0.0
    };
}

size_t find_last_dof_index(const ConstraintEQ& c);

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
    size_t constrained_dof;
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

typedef std::unordered_map<int,RearrangedConstraintEQ> ConstraintMatrix;

RearrangedConstraintEQ isolate_term_on_lhs(const ConstraintEQ& c, 
    size_t constrained_index);

ConstraintEQ substitute(const ConstraintEQ& c, size_t constrained_dof_index,
    const RearrangedConstraintEQ& subs_in);

bool is_constrained(const ConstraintMatrix& dof_constraint_map, size_t dof);

RearrangedConstraintEQ make_lower_triangular(const ConstraintEQ& c,
    const ConstraintMatrix& map);

ConstraintMatrix from_constraints(const std::vector<ConstraintEQ>& constraints);

/* Accepts a reduced DOF vector and returns the full DOF vector. */
std::vector<double> distribute_vector(const ConstraintMatrix& matrix, 
    const std::vector<double>& in, size_t total_dofs);

/* Accepts a full DOF vector and returns the reduced DOF vector.  */
std::vector<double> condense_vector(const ConstraintMatrix& matrix,
    const std::vector<double>& all);

struct Matrix {
    const size_t n_rows;
    const size_t n_cols;
    const std::vector<double> data;
};

struct MatrixEntry 
{
    const size_t loc[2];
    const double value;
};

Matrix
condense_matrix(const ConstraintMatrix& constraint_matrix, const Matrix& matrix);

} // END namespace tbem

#endif
