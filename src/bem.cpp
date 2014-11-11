#include <cassert>
#include "bem.h"
#include "mesh.h"
#include "adaptive_quad.h"

FaceInfo::FaceInfo(const Mesh3D& mesh, int face_index):
    face_index(face_index),
    face(mesh.faces[face_index]),
    corners(index3(mesh.vertices, face)),
    unscaled_normal(tri_unscaled_normal(corners)),
    area(tri_area(unscaled_normal)),
    jacobian(area * 2.0),
    normal(unscaled_normal / jacobian)
{}

ObsPt ObsPt::from_face(const QuadRule2d& obs_quad,
                              const FaceInfo& obs_face, int idx) {
    return {
        //TODO: need to divide by q or something like that
        std::sqrt(obs_face.area),
        ref_to_real(obs_quad[idx].x_hat, obs_face.corners),
        obs_face.normal 
    };
}

Vec3<double> eval_quad_pt(const std::array<double,2>& x_hat,
                          const Kernel& kernel,
                          const FaceInfo& face,
                          const Vec3<double>& obs_loc,
                          const Vec3<double>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.corners);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return (kernel_val * face.jacobian) * linear_basis(x_hat);
}

double appx_face_dist2(const Vec3<double>& pt,
                       const std::array<Vec3<double>,3>& vs) {
    double d0 = dist2(pt, vs[0]);
    double d1 = dist2(pt, vs[1]);
    double d2 = dist2(pt, vs[2]);
    return std::min(d0, std::min(d1, d2));
}

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

static unsigned int interacts = 0;
static unsigned far_field_pairs = 0;
static unsigned near_field_pairs = 0;
static unsigned adjacent_pairs = 0;

Vec3<double> near_field(const Problem& p, const QuadStrategy& qs,
                        const ObsPt& obs, const FaceInfo& src_face,
                        const double dist2) {
    std::vector<Vec3<double>> near_steps(qs.n_singular_steps, {0,0,0});
    for (int nf = 0; nf < qs.n_singular_steps; nf++) {
        double nfdn = 5 * obs.len_scale * qs.singular_steps[nf];
        auto nf_obs_pt = obs.loc + nfdn * obs.normal;
         
        if (dist2 > 0.5 * src_face.area) {
            for (unsigned int i = 0; i < qs.src_near_quad.size(); i++) {
                near_steps[nf] += qs.src_near_quad[i].w *
                    eval_quad_pt(qs.src_near_quad[i].x_hat, p.K, src_face,
                                 nf_obs_pt, obs.normal);
            }
            interacts += qs.src_near_quad.size();
            near_field_pairs++;
        } else {
            Vec3<double> ns = adaptive_integrate<Vec3<double>>(
                [&] (double x_hat) {
                    if (x_hat == 1.0) {
                        return Vec3<double>{0.0,0.0,0.0};
                    }
                    return adaptive_integrate<Vec3<double>>([&] (double y_hat) {
                            Vec3<double> val  = eval_quad_pt({x_hat, y_hat}, p.K, src_face,
                                                nf_obs_pt, obs.normal);
                            interacts++;
                            return val;
                        }, 0.0, 1 - x_hat, qs.singular_tol);
                }, 0.0, 1.0, qs.singular_tol);
            near_steps[nf] += ns;
            adjacent_pairs++;
        }
    }
    return richardson_step(near_steps);
}

std::vector<double> integral_equation_vector(const Problem& p,
                                             const QuadStrategy& qs,
                                             const ObsPt& obs) {
    std::vector<double> result(p.src_mesh.vertices.size(), 0.0);
    for (unsigned int i = 0; i < p.src_mesh.faces.size(); i++) {
        FaceInfo src_face(p.src_mesh, i);
        const double dist2 = appx_face_dist2(obs.loc, src_face.corners);

        Vec3<double> integrals;
        if (dist2 < pow(qs.far_threshold, 2) * src_face.area) {
            integrals = near_field(p, qs, obs, src_face, dist2);
        } else {
            // farfield
            integrals = {0.0, 0.0, 0.0};
            for (unsigned int i = 0; i < qs.src_far_quad.size(); i++) {
                integrals += qs.src_far_quad[i].w *
                    eval_quad_pt(qs.src_far_quad[i].x_hat, p.K, src_face,
                                 obs.loc, obs.normal);
            }
            interacts += qs.src_far_quad.size();
            far_field_pairs++;
        }
        for (int b = 0; b < 3; b++) {
            result[src_face.face[b]] += integrals[b];
        }
    }
    return result;
}

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(const Problem& p, const QuadStrategy& qs,
                              const ObsPt& obs) {
    double result = 0.0;
    auto row = integral_equation_vector(p, qs, obs);
    for (unsigned int i = 0; i < row.size(); i++) {
        result += row[i] * p.src_strength[i];
    }
    return result;
}

//TODO: Use a sparse matrix storage format here.
std::vector<std::vector<double>> interact_matrix(const Problem& p,
                                                 const QuadStrategy& qs) {
    int n_obs_basis = p.obs_mesh.vertices.size();
    int n_src_basis = p.src_mesh.vertices.size();
    std::vector<std::vector<double>> matrix(n_obs_basis,
            std::vector<double>(n_src_basis, 0.0));
// #pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt::from_face(qs.obs_quad, obs_face, obs_q);

            const auto row = integral_equation_vector(p, qs, pt);

            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp(qs.obs_quad[obs_q].x_hat,
                                                      unit<double>(v)); 
                int b = obs_face.face[v];
#pragma omp critical
                for (int i = 0; i < n_src_basis; i++) {
                    matrix[b][i] += obs_basis_eval * row[i] * qs.obs_quad[obs_q].w * obs_face.jacobian;
                }
            }
        }
    }
    std::cout << interacts << " " << far_field_pairs << " " << near_field_pairs << " " << adjacent_pairs << std::endl;
    return matrix;
}

std::vector<double> bem_mat_mult(const Mat& A, const std::vector<double>& x) {
    std::vector<double> res(A.size(), 0.0);
    for (unsigned int i = 0; i < A.size(); i++) {
        for (unsigned int j = 0; j < A[i].size(); j++) {
            res[i] += A[i][j] * x[j]; 
        }
    }
    return res;
}

std::vector<double> direct_interact(const Problem& p,
                                    const QuadStrategy& qs) {

    return bem_mat_mult(interact_matrix(p, qs), p.src_strength);
}

std::vector<double> mass_term(const Problem& p,
                              const QuadStrategy& qs) {
    std::vector<double> integrals(p.obs_mesh.vertices.size(), 0.0);
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];
            auto face_vals = index3(p.src_strength, obs_face.face);
            double interp_val = linear_interp(qpt.x_hat, face_vals);

            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp(qpt.x_hat, unit<double>(v)); 
                integrals[obs_face.face[v]] += obs_face.jacobian * obs_basis_eval * 
                                          interp_val * qpt.w;
            }
        }
    }
    return integrals;
}

double get_len_scale(Mesh3D& mesh, int which_face, int q) {
    auto face = mesh.faces[which_face];
    return std::sqrt(tri_area(index3(mesh.vertices, face))) / q;
}
