#include "matrix_free_operator.h"
#include "petsc_facade.h"
#include "vectorx.h"
#include <iostream>
    
namespace tbem {

MatrixFreeOperator::MatrixFreeOperator(NearfieldPtr nearfield):
    nearfield(std::move(nearfield)),
    shape{this->nearfield->n_rows(), this->nearfield->n_cols()}
{}

MatrixFreeOperator::MatrixFreeOperator(NearfieldPtr nearfield,
        OperatorI::Ptr farfield_op):
    nearfield(std::move(nearfield)),
    shape{this->nearfield->n_rows(), this->nearfield->n_cols()},
    farfield_op(std::move(farfield_op))
{}

size_t MatrixFreeOperator::n_rows() const
{
    return shape.n_rows;
}

size_t MatrixFreeOperator::n_cols() const
{
    return shape.n_cols;
}

VectorX MatrixFreeOperator::apply(const VectorX& x) const 
{
    auto res = nearfield->mat_vec_prod(x);
    if (farfield_op) {
        return res + farfield_op->apply(x);
    } else {
        return res;
    }
}

} //end namespace tbem
