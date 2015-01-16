#ifndef __rH1JH11LHLHAHAh_FUNCTION_H
#define __rH1JH11LHLHAHAh_FUNCTION_H

#include <vector>
#include <cstdlib>

namespace tbem {

struct ConcatenatedFunction 
{
    const size_t components;
    const std::vector<double> data;
    const std::vector<size_t> component_lengths;
};

ConcatenatedFunction concatenate(const std::vector<std::vector<double>>& fncs); 

std::vector<std::vector<double>> expand(const ConcatenatedFunction& block_fnc,
    const std::vector<double>& replacement_data);

std::vector<std::vector<double>> expand(const ConcatenatedFunction& block_fnc); 

} // end namespace tbem

#endif
