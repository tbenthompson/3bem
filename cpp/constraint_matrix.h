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

ConstraintEQ make_lower_triangular(const ConstraintEQ& c, const ConstraintMatrix& map);

ConstraintMatrix from_constraints(const std::vector<ConstraintEQ>& constraints);

/* Accepts a reduced DOF vector and returns the full DOF vector. 
 */
std::vector<double> distribute_vector(const ConstraintMatrix& matrix, 
    const std::vector<double>& in, size_t total_dofs);

/* Accepts a full DOF vector and returns the reduced DOF vector.  
 */
std::vector<double> condense_vector(const ConstraintMatrix& matrix,
    const std::vector<double>& all);

/* Condenses a matrix. Note that this function and the equivalent that operates
 * on sparse matrices assume that all constraints are homogeneous -- that the
 * rhs of the constraint is zero. To deal with inhomogeneous, a preprocessing 
 * step can be performed by multiplying the uncondensed matrix by a vector 
 * resulting from applying distribute_vector to the all-zeros vector. This
 * should represent the influence of the inhomogeneous component of the 
 * constraints. As a result, adding this term to the rhs will allow treating
 * all constraints as homogeneous.
 */
DenseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const DenseOperator& matrix);

/* See note for the dense version 
 */
SparseOperator condense_matrix(const ConstraintMatrix& row_cm,
    const ConstraintMatrix& col_cm, const SparseOperator& matrix);

/* Set the rhs of each constraint equation to 0. */
ConstraintMatrix homogenize_constraints(const ConstraintMatrix& cm);

} // end namespace tbem

#endif

