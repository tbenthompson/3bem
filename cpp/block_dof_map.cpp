#include "block_dof_map.h"

namespace tbem {

size_t BlockDOFMap::get_past_end_dof(size_t component) const {
    assert(component < n_components);
    if (component + 1 == n_components) {
        return n_dofs;
    } else {
        return start_positions[component + 1];
    }
}

BlockDOFMap build_block_dof_map(std::vector<size_t> component_sizes) {
    size_t n_cumulative_dofs = 0;
    std::vector<size_t> start_positions;
    for (size_t i = 0; i < component_sizes.size(); i++) {
        start_positions.push_back(n_cumulative_dofs);
        n_cumulative_dofs += component_sizes[i];
    }
    return {
        component_sizes.size(),
        n_cumulative_dofs,
        start_positions
    };
}

} //end namespace tbem
