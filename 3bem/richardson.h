#ifndef __DSHLHJASDHJ_RICHARDSON_H
#define __DSHLHJASDHJ_RICHARDSON_H

#include <vector>
#include <cassert>

namespace tbem {

template <typename T>
T richardson_limit(const std::vector<T>& values) {
    assert(values.size() > 1);

    int n_steps = values.size();
    std::vector<T> last_level = values;
    std::vector<T> this_level;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, m);
            const double factor = 1.0 / (mult - 1.0);
            T low = last_level[i];
            T high = last_level[i + 1];
            T moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        last_level = this_level;
    }
    return this_level[0];
}

} //end namespace

#endif
