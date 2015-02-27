#ifndef __ALJLQH1123_MATRIX_FREE_OPERATOR_H
#define __ALJLQH1123_MATRIX_FREE_OPERATOR_H

#include <memory>
#include "operator.h"
#include "fwd_vectorx.h"

namespace tbem {

struct PETScMatWrapper;

struct MatrixFreeOperator: public OperatorI {
    const OperatorShape shape;
    std::unique_ptr<PETScMatWrapper> nearfield;

    MatrixFreeOperator(size_t n_rows, size_t n_cols);

    virtual size_t n_rows() const override;
    virtual size_t n_cols() const override;
    virtual VectorX apply(const VectorX& x) const override;
};

} // end namespace tbem

#endif
