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
template Vec3<double> richardson_step(const std::vector<Vec3<double>>&);
template Vec2<double> richardson_step(const std::vector<Vec2<double>>&);

std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x) {
    std::vector<double> res(n_obs_dofs, 0.0);
#pragma omp parallel for
    for (int i = 0; i < n_obs_dofs; i++) {
        for (std::size_t j = 0; j < x.size(); j++) {
            res[i] += A[i * x.size() + j] * x[j]; 
        }
    }
    return res;
}

std::vector<double> mass_term(const Problem<3>& p,
                              const QuadStrategy<3>& qs) {
    int n_obs_dofs = 3 * p.obs_mesh.facets.size();
    std::vector<double> integrals(n_obs_dofs, 0.0);
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        FaceInfo<3> obs_face(p.obs_mesh.facets[obs_idx]);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];
            int dof = 3 * obs_idx;
            Vec3<double> face_vals = {
                p.src_strength[dof], p.src_strength[dof + 1], p.src_strength[dof + 2]
            };
            double interp_val = linear_interp<3>(qpt.x_hat, face_vals);

            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp<3>(qpt.x_hat, unit<double,3>(v)); 
                integrals[dof + v] += obs_face.jacobian * obs_basis_eval * 
                                      interp_val * qpt.w;
            }
        }
    }
    return integrals;
}

double get_len_scale(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face].vertices)) / q;
}
