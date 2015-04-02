#ifndef __NMMMNBNBBHSKSKS_KERNEL_H
#define __NMMMNBNBBHSKSKS_KERNEL_H

#include "vec.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct Kernel {
    const static size_t n_rows = R;
    const static size_t n_cols = C;
    typedef Vec<double,R> OutType;
    typedef Vec<double,C> InType;
    typedef Vec<Vec<double,C>,R> OperatorType;

    virtual OperatorType operator()(double r2, const Vec<double,dim>& delta, 
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const = 0;
};

} //End namespace tbem

#endif
