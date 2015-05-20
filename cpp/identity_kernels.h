#ifndef __HDHDHDHHDHDHDH_IDENTITY_KERNELS_H
#define __HDHDHDHHDHDHDH_IDENTITY_KERNELS_H

#include "kernel.h"
#include "vec.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct IdentityTensor: public Kernel<dim,R,C>
{
    typedef Vec<Vec<double,C>,R> OperatorType;
    OperatorType call_with_no_params() const {
        auto out = zeros<OperatorType>::make();
        for (size_t i = 0; i < R; i++) {
            for (size_t j = 0; j < C; j++) {
                if (i == j) {
                    out[i][j] = 1.0;
                }
            }
        }
        return out;
    }

    OperatorType call(double r2, const Vec<double,dim>& delta,
        const Vec<double,dim>& nobs, const Vec<double,dim>& nsrc) const 
    {
        (void)r2; (void)delta; (void)nobs; (void)nsrc;
        return call_with_no_params();
    }

    virtual std::unique_ptr<Kernel<dim,R,C>> clone() const
    {
        return std::unique_ptr<Kernel<dim,R,C>>(new IdentityTensor<dim,R,C>());
    }
};

template <size_t dim>
using IdentityScalar = IdentityTensor<dim,1,1>;

}//end namespace tbem

#endif
