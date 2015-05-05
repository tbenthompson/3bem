#include "nbody_operator.h"
#include "vec_ops.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
std::vector<double> BlockDirectNBodyOperator<dim,R,C>::apply(
    const std::vector<double>& x) const
{
    std::vector<double> out(R * data.obs_locs.size(), 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto kernel_val = K(data.obs_locs[i], data.src_locs[j],
                data.obs_normals[i], data.src_normals[j]);
            auto entry = data.src_weights[j] * kernel_val;
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    auto& out_val = out[d1 * data.obs_locs.size() + i];
                    auto in = x[d2 * data.src_locs.size() + j];
                    out_val += entry[d1][d2] * in;
                }
            }
        }
    }
    return out;
}

template struct BlockDirectNBodyOperator<2,1,1>;
template struct BlockDirectNBodyOperator<2,2,2>;
template struct BlockDirectNBodyOperator<3,1,1>;
template struct BlockDirectNBodyOperator<3,3,3>;

} // end namespace tbem
