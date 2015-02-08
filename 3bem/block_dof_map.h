#ifndef __12345667890123_BLOCK_DOF_MAP_H
#define __12345667890123_BLOCK_DOF_MAP_H

#include <vector>
#include <cstdlib>
#include <cassert>
typedef std::vector<double> Function;

namespace tbem {

struct BlockDOFMap {
    const size_t n_components;
    const size_t n_dofs;
    const std::vector<size_t> start_positions;

    size_t get_past_end_dof(size_t component) const;
};

BlockDOFMap build_block_dof_map(std::vector<size_t> component_sizes);

BlockDOFMap block_dof_map_from_functions(const std::vector<Function>& fncs);

Function 
concatenate(const BlockDOFMap& dof_map, const std::vector<Function>& fncs);

std::vector<Function>
expand(const BlockDOFMap& dof_map, const std::vector<double>& data);

} // end namespace tbem

#endif
