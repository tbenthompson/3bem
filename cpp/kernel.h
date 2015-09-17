#ifndef TBEMNMMMNBNBBHSKSKS_KERNEL_H
#define TBEMNMMMNBNBBHSKSKS_KERNEL_H

#include "vec_ops.h"
#include "geometry.h"
#include <memory>

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct Kernel {
    typedef Vec<double,R> OutType;
    typedef Vec<double,C> InType;
    typedef Vec<Vec<double,C>,R> OperatorType;

    const static size_t n_rows = R;
    const static size_t n_cols = C;

    OperatorType operator()(
        const Vec<double,dim>& obs_pt, 
        const Vec<double,dim>& src_pt, 
        const Vec<double,dim>& obs_normal, 
        const Vec<double,dim>& src_normal) const 
    {
        const auto d = src_pt - obs_pt;
        const auto r2 = dot_product(d, d);
        //TODO: 1e-12 should be dependent on length scale... 
        if (r2 < 1e-12) {
            //TODO: log something in this circumstance?
            return zeros<OperatorType>::make();
        }
        return call(r2, d, obs_normal, src_normal);
    }

    virtual OperatorType call(double r2, const Vec<double,dim>& delta, 
        const Vec<double,dim>& nobs, const Vec<double,dim>& nsrc) const = 0;

    virtual std::unique_ptr<Kernel<dim,R,C>> clone() const = 0;
};


} //End namespace tbem

#endif
