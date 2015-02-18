#ifndef __NMMMNBNBBHSKSKS_KERNEL_H
#define __NMMMNBNBBHSKSKS_KERNEL_H

#include "vec.h"

namespace tbem {

template <size_t dim, size_t n_rows, size_t n_cols>
struct Kernel {
    typedef Vec<double,n_rows> OutType;
    typedef Vec<double,n_cols> InType;
    typedef Vec<Vec<double,n_cols>,n_rows> OperatorType;

    virtual OperatorType operator()(double r2, const Vec<double,dim>& delta, 
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const = 0;
};

} //End namespace tbem

#endif
