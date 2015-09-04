#ifndef TBEM_CONTAINING_TRIANGLE_H
#define TBEM_CONTAINING_TRIANGLE_H

#include <vector>
#include "vec.h"

namespace tbem {

std::pair<bool,size_t> find_containing_tri_idx(const Vec<double,2>& query_pt,
    const std::vector<Vec<size_t,3>>& tris,
    const std::vector<Vec<double,2>>& pts);

} //end namespace tbem

#endif
