#ifndef TBEM_NEW_LAPLACE_KERNELS_H
#define TBEM_NEW_LAPLACE_KERNELS_H

#include "new_kernel.h"
#include "geometry.h"

namespace tbem {

<%
    name = 'NEWLaplaceSingle'
    component_rows = 1
    component_cols = 1
%>
template <size_t dim>
struct ${name};

<%
    dim = 2
    import math
    import tbempy.templating
    denom = 1.0 / (2 * math.pi)
    expr = 'std::log(r) * ' + str(denom)
%>

template <>
struct ${name}<${dim}>: public NEWKernel<${dim}> 
{
    std::vector<double> operator()(
        const std::vector<Vec<double,${dim}>>& obs_pts, 
        const std::vector<Vec<double,${dim}>>& src_pts, 
        const std::vector<Vec<double,${dim}>>& obs_normals, 
        const std::vector<Vec<double,${dim}>>& src_normals) const
    {
        (void)obs_normals;
        (void)src_normals;
        size_t n_obs = obs_pts.size();
        size_t n_src = src_pts.size();
        std::vector<double> out_matrix(
            n_obs * n_src * ${component_rows} * ${component_cols}
        );
        auto component_entries = ${component_rows} * ${component_cols};
        for (size_t i = 0; i < n_obs; i++) {
            auto row_start_index = i * n_src * component_entries;
            for (size_t j = 0; j < n_src; j++) {
                auto col_offset = j * component_entries;
                auto r = dist(obs_pts[i], src_pts[j]);
                for (size_t d1 = 0; d1 < ${component_rows}; d1++) {
                    for (size_t d2 = 0; d2 < ${component_cols}; d2++) {
                        auto component_offset = d1 * ${component_cols} + d2;
                        auto index = row_start_index + col_offset + component_offset;
                        out_matrix[index] = ${expr};
                    }
                }
            }
        }
        return out_matrix;
    }

    virtual size_t n_component_rows() const {return 1;}
    virtual size_t n_component_cols() const {return 1;}

    virtual std::unique_ptr<NEWKernel<${dim}>> clone() const
    {
        return std::unique_ptr<NEWKernel<${dim}>>(new ${name}<2>());
    }
};

} //end namespace tbem

#endif
