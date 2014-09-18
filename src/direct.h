#ifndef __DIRECT_H
#define __DIRECT_H
#include <array>
#include <vector>
#include <functional>

std::vector<double> direct_n_body(std::array<std::vector<double>,3>& src_locs,
                                  std::array<std::vector<double>,3>& obs_locs,
                                  std::function<double (double, double, double, 
                                                        double, double, double)> kernel,
                                  std::vector<double>& values) 
{
    std::vector<double> out_vals(obs_locs[0].size());
    for (unsigned int i = 0; i < obs_locs[0].size(); ++i) {
        out_vals[i] = 0.0;
        for (unsigned int j = 0; j < src_locs[0].size(); ++j) {
            out_vals[i] += kernel(obs_locs[0][i], obs_locs[1][i], obs_locs[2][i],
                                  src_locs[0][j], src_locs[1][j], src_locs[2][j]);
        }
    }
    return out_vals;
}

#endif
