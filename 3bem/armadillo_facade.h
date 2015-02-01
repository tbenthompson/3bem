#ifndef __kljewtlkha_ARMADILLO_H
#define __kljewtlkha_ARMADILLO_H

#include <vector>
#include <cstdlib>

namespace tbem {

struct Operator;
Operator arma_invert(const Operator& op);
double arma_cond(const Operator& op);

} // end namespace tbem

#endif
