#include <cassert>
#include "function.h"

namespace tbem {

ConcatenatedFunction concatenate(const std::vector<std::vector<double>>& fncs) 
{
    std::vector<size_t> component_lengths;
    size_t total_length = 0;
    for (const auto& f: fncs) {
        total_length += f.size();
        component_lengths.push_back(f.size());
    }

    std::vector<double> data(total_length);
    size_t fnc_begin_idx = 0;
    for (const auto& f: fncs) {
        for (size_t i = 0; i < f.size(); i++) {
            data[fnc_begin_idx + i] = f[i];
        }
        fnc_begin_idx += f.size(); 
    }
    return {fncs.size(), data, component_lengths}; 
}

std::vector<std::vector<double>> expand(const ConcatenatedFunction& block_fnc,
    const std::vector<double>& replacement_data) 
{
    assert(replacement_data.size() == block_fnc.data.size());

    std::vector<std::vector<double>> out(block_fnc.components);
    size_t start_pos = 0;
    for (size_t i = 0; i < block_fnc.components; i++) {
        auto this_comp_len = block_fnc.component_lengths[i];
        auto end_pos = start_pos + this_comp_len;
        out[i].resize(this_comp_len); 
        std::copy(
            replacement_data.begin() + start_pos,
            replacement_data.begin() + end_pos,
            out[i].begin()
        );
        start_pos += this_comp_len;
    }
    return out;
}

std::vector<std::vector<double>> expand(const ConcatenatedFunction& block_fnc) 
{
    return expand(block_fnc, block_fnc.data);
}

} // end namespace tbem
