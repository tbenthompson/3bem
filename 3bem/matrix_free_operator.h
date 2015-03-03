#ifndef __ALJLQH1123_MATRIX_FREE_OPERATOR_H
#define __ALJLQH1123_MATRIX_FREE_OPERATOR_H

#include <memory>
#include "operator.h"
#include "fwd_vectorx.h"

namespace tbem {

template <typename T>
struct BlockOperator;

struct PETScSparseMatWrapper;

struct MatrixFreeOperator: public OperatorI {
    typedef std::unique_ptr<const PETScSparseMatWrapper> NearfieldPtr;
    typedef std::function<VectorX(const VectorX&)> FarfieldFunc;

    NearfieldPtr nearfield;
    const OperatorShape shape;
    OperatorI::Ptr farfield_op;

    MatrixFreeOperator(NearfieldPtr nearfield);
    MatrixFreeOperator(NearfieldPtr nearfield, OperatorI::Ptr farfield_op);

    virtual size_t n_rows() const override;
    virtual size_t n_cols() const override;
    virtual VectorX apply(const VectorX& x) const override;
};

typedef BlockOperator<MatrixFreeOperator> BlockMatrixFreeOperator;

} // end namespace tbem

#endif
