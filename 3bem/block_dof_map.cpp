#include "block_dof_map.h"
#include "vectorx.h"

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

BlockDOFMap block_dof_map_from_functions(const BlockVectorX& fncs) {
    std::vector<size_t> sizes;
    for (const auto& f: fncs) {
        sizes.push_back(f.size());
    }
    return build_block_dof_map(sizes);
}

VectorX 
concatenate(const BlockDOFMap& dof_map, const BlockVectorX& fncs) 
{
    assert(fncs.size() == dof_map.n_components);

    VectorX out(dof_map.n_dofs);
    for (size_t f_idx = 0; f_idx < dof_map.n_components; f_idx++) {
        auto start = dof_map.start_positions[f_idx];
        for (size_t i = 0; i < fncs[f_idx].size(); i++) {
            out[start + i] = fncs[f_idx][i];
        }
    }
    
    return out;
}

BlockVectorX expand(const BlockDOFMap& dof_map, const VectorX& data)
{
    BlockVectorX out(dof_map.n_components);

    for (size_t i = 0; i < dof_map.n_components; i++) {
        auto start_dof = dof_map.start_positions[i];
        auto past_end_dof = dof_map.get_past_end_dof(i);
        out[i].resize(past_end_dof - start_dof);
        std::copy(data.begin() + start_dof, data.begin() + past_end_dof, out[i].begin());
    }

    return out;
}

} //end namespace tbem
