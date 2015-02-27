#include "matrix_free_operator.h"
#include "petsc_facade.h"
#include "vectorx.h"
    
namespace tbem {

MatrixFreeOperator::MatrixFreeOperator(size_t n_rows, size_t n_cols):
    shape{n_rows, n_cols}
{
    return;
}

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
    return x;
}

} //end namespace tbem
