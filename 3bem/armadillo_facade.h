#ifndef __kljewtlkha_ARMADILLO_H
#define __kljewtlkha_ARMADILLO_H

#include <vector>
#include <cstdlib>

namespace tbem {

struct DenseOperator;
DenseOperator arma_invert(const DenseOperator& op);
double arma_cond(const DenseOperator& op);

} // end namespace tbem

#endif
