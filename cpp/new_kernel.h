#ifndef TBEMNMMMNBNBBHSKSKS_NEW_KERNEL_H
#define TBEMNMMMNBNBBHSKSKS_NEW_KERNEL_H

#include "vec.h"
#include <memory>

namespace tbem {

template <size_t dim>
struct NEWKernel {

    std::vector<double> operator()(
        const std::vector<Vec<double,dim>>& obs_pts, 
        const std::vector<Vec<double,dim>>& src_pts, 
        const std::vector<Vec<double,dim>>& obs_normals, 
        const std::vector<Vec<double,dim>>& src_normals) const;

    virtual size_t n_component_rows() const = 0;
    virtual size_t n_component_cols() const = 0;

    virtual std::unique_ptr<NEWKernel<dim>> clone() const = 0;
};

} //End namespace tbem

#endif
