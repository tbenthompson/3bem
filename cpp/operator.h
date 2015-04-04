#ifndef ADKJASdLL_OPERATOR_H
#define ADKJASdLL_OPERATOR_H

#include "vectorx.h"

namespace tbem {

struct OperatorShape {
    const size_t n_rows;
    const size_t n_cols;
};

struct OperatorI {
    virtual size_t n_rows() const = 0;
    virtual size_t n_cols() const = 0;
    virtual VectorX apply(const VectorX& x) const = 0;
};

} // end namespace tbem

#endif
