#ifndef __JLQJWE67155151_NBODY_OPERATOR_H
#define __JLQJWE67155151_NBODY_OPERATOR_H

#include <cstdlib>
#include <vector>
#include "vec.h"
#include "kernel.h"
#include "operator.h"
#include "fwd_vectorx.h"
#include "nbody_data.h"

namespace tbem {

template <size_t dim, size_t R, size_t C>
struct BlockDirectNBodyOperator {
    const NBodyData<dim> data;
    const Kernel<dim,R,C>& K;

    BlockVectorX apply(const BlockVectorX& x) const;
};

} // end namespace tbem

#endif