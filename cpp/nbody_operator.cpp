#include "nbody_operator.h"
#include "vectorx.h"
#include "vec_ops.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
BlockVectorX BlockDirectNBodyOperator<dim,R,C>::apply(const BlockVectorX& x) const
{
    BlockVectorX out(R, VectorX(data.obs_locs.size(), 0.0));
#pragma omp parallel for
    for (size_t i = 0; i < data.obs_locs.size(); i++) {
        for (size_t j = 0; j < data.src_locs.size(); j++) {
            auto d = data.src_locs[j] - data.obs_locs[i];
            auto r2 = dot_product(d, d);
            auto kernel_val = K(r2, d, data.src_normals[j], data.obs_normals[i]);
            auto entry = data.src_weights[j] * kernel_val;
            for (size_t d1 = 0; d1 < R; d1++) {
                for (size_t d2 = 0; d2 < C; d2++) {
                    out[d1][i] += entry[d1][d2] * x[d2][j];
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
