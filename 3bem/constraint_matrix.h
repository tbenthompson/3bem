#ifndef __ASDASDQ12123_CONSTRAINT_MATRIX_H
#define __ASDASDQ12123_CONSTRAINT_MATRIX_H
#include <vector>
#include <unordered_map>
#include "constraint.h"
#include "operator.h"
#include "function.h"

namespace tbem {

typedef std::unordered_map<int,RearrangedConstraintEQ> ConstraintMatrix;

bool is_constrained(const ConstraintMatrix& dof_constraint_map, size_t dof);

RearrangedConstraintEQ make_lower_triangular(const ConstraintEQ& c,
    const ConstraintMatrix& map);

ConstraintMatrix from_constraints(const std::vector<ConstraintEQ>& constraints);

/* Accepts a reduced DOF vector and returns the full DOF vector. */
Function distribute_vector(const ConstraintMatrix& matrix, 
    const Function& in, size_t total_dofs);

/* Accepts a full DOF vector and returns the reduced DOF vector.  */
Function condense_vector(const ConstraintMatrix& matrix, const Function& all);

struct MatrixEntry 
{
    const size_t loc[2];
    const double value;
};

Operator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const Operator& matrix);

BlockOperator condense_block_operator(const std::vector<ConstraintMatrix>& row_cms,
    const std::vector<ConstraintMatrix>& col_cms, const BlockOperator& op);

} // end namespace tbem

#endif

