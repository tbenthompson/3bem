#ifndef __HDHDHDHHDHDHDH_IDENTITY_KERNELS_H
#define __HDHDHDHHDHDHDH_IDENTITY_KERNELS_H

#include "vec.h"

namespace tbem {

template <size_t dim>
struct IdentityScalar {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    OperatorType call_with_no_params() const {
        return 1.0;
    }

    OperatorType operator()(const double& r2, const Vec<double,dim>& delta,
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const 
    {
        return call_with_no_params();
    }
};

template <size_t dim, size_t n_rows, size_t n_cols>
struct IdentityTensor {
    typedef Vec<double,n_rows> OutType;
    typedef Vec<double,n_cols> InType;
    typedef Vec<Vec<double,n_cols>,n_rows> OperatorType;

    OperatorType call_with_no_params() const {
        auto out = zeros<OperatorType>::make();
        for (size_t i = 0; i < n_rows; i++) {
            for (size_t j = 0; j < n_cols; j++) {
                if (i == j) {
                    out[i][j] = 1.0;
                }
            }
        }
        return out;
    }

    OperatorType operator()(const double& r2, const Vec<double,dim>& delta,
        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const 
    {
        return call_with_no_params();
    }
};

}//end namespace tbem

#endif
