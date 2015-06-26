#ifndef TBEMDSHLHJASDHJ_RICHARDSON_H
#define TBEMDSHLHJASDHJ_RICHARDSON_H

#include <vector>
#include <cmath>
#include <cassert>

namespace tbem {

template <typename T>
T richardson_limit(double step_ratio, const std::vector<T>& values) 
{
    assert(values.size() > 1);

    auto n_steps = values.size();
    auto last_level = values;
    decltype(last_level) this_level;

    for (size_t m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);

        for (size_t i = 0; i < n_steps - m; i++) {
            auto mult = std::pow(step_ratio, m);
            auto factor = 1.0 / (mult - 1.0);
            auto low = last_level[i];
            auto high = last_level[i + 1];
            auto moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        last_level = this_level;
    }
    return this_level[0];
}

} //end namespace

#endif
