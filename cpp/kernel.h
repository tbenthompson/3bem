#ifndef __NMMMNBNBBHSKSKS_KERNEL_H
#define __NMMMNBNBBHSKSKS_KERNEL_H

#include "vec_ops.h"
#include "geometry.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct Kernel {
    const static size_t n_rows = R;
    const static size_t n_cols = C;
    typedef Vec<double,R> OutType;
    typedef Vec<double,C> InType;
    typedef Vec<Vec<double,C>,R> OperatorType;

    virtual OperatorType operator()(
        const Vec<double,dim>& obs_pt, 
        const Vec<double,dim>& src_pt, 
        const Vec<double,dim>& obs_normal, 
        const Vec<double,dim>& src_normal) const 
    {
        const auto d = src_pt - obs_pt;
        const auto r2 = dot_product(d, d);
        return call(r2, d, obs_normal, src_normal);
    }


    virtual OperatorType call(double r2, const Vec<double,dim>& delta, 
        const Vec<double,dim>& nobs, const Vec<double,dim>& nsrc) const = 0;
};

} //End namespace tbem

#endif
