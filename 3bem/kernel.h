#ifndef __NMMMNBNBBHSKSKS_KERNEL_H
#define __NMMMNBNBBHSKSKS_KERNEL_H

#include "vec.h"

namespace tbem {

template <size_t dim, typename O, typename I, typename Op>
struct Kernel {
    typedef O OutType;
    typedef I InType;
    typedef Op OperatorType;

    virtual OperatorType operator()(double r2, const Vec<double,dim>& delta, 
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const = 0;
};

template <size_t dim>
using ScalarKernel = Kernel<dim,Vec1<double>,Vec1<double>,Vec1<Vec1<double>>>;

} //End namespace tbem

#endif
