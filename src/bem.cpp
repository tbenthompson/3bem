#include <cassert>
#include "bem.h"


/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
template <typename T>
T richardson_step(const std::vector<T>& values) {
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
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
}

template double richardson_step(const std::vector<double>&);
template Vec2<double> richardson_step(const std::vector<Vec2<double>>&);
template Vec3<double> richardson_step(const std::vector<Vec3<double>>&);

std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x) {
    assert(n_obs_dofs * x.size() == A.size());
    std::vector<double> res(n_obs_dofs, 0.0);
#pragma omp parallel for
    for (int i = 0; i < n_obs_dofs; i++) {
        for (std::size_t j = 0; j < x.size(); j++) {
            res[i] += A[i * x.size() + j] * x[j]; 
        }
    }
    return res;
}
