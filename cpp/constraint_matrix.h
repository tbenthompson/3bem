#ifndef TBEMASDASDQ12123_CONSTRAINT_MATRIX_H
#define TBEMASDASDQ12123_CONSTRAINT_MATRIX_H
#include <vector>
#include <unordered_map>
#include "constraint.h"
#include "dense_operator.h"
#include "sparse_operator.h"

namespace tbem {

typedef std::unordered_map<int,RearrangedConstraintEQ> ConstraintMatrix;

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

DenseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const DenseOperator& matrix);

SparseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const SparseOperator& matrix);

} // end namespace tbem

#endif

