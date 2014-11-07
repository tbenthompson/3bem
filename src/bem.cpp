#include <cassert>
#include "bem.h"
#include "mesh.h"

FaceInfo::FaceInfo(const Mesh& mesh, int face_index):
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
        ref_to_real(obs_quad[idx].x_hat[0], obs_quad[idx].x_hat[1], obs_face.corners),
        obs_face.normal 
    };
}

SrcPointInfo::SrcPointInfo(const QuadRule2d& quad_rule,
              const Kernel& kernel,
              const FaceInfo& face,
              const Vec3<double>& obs_loc,
              const Vec3<double>& obs_n,
              int q_index) {
    auto qpt = quad_rule[q_index];
    basis = linear_basis(qpt.x_hat[0], qpt.x_hat[1]);

    const auto src_pt = ref_to_real(qpt.x_hat[0], qpt.x_hat[1], face.corners);
    const auto d = src_pt - obs_loc;
    const auto r2 = hypot2(d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);

    weighted_kernel = kernel_val * qpt.w * face.jacobian;
}


Vec3<double> basis_integrals(const QuadRule2d& quad_rule,
                        const Kernel& kernel,
                        const FaceInfo& face,
                        const Vec3<double>& obs_loc,
                        const Vec3<double>& obs_n) {

    Vec3<double> result = {0,0,0};
    for (unsigned int src_q = 0; src_q < quad_rule.size(); src_q++) {
        SrcPointInfo pt(quad_rule, kernel, face, obs_loc, obs_n, src_q);
        result += pt.weighted_kernel * pt.basis;
    }
    return result;
}

double integral(const QuadRule2d& quad_rule,
           const Kernel& kernel,
           const FaceInfo& face,
           const Vec3<double>& src_vals,
           const Vec3<double>& obs_loc,
           const Vec3<double>& obs_n) {
    auto basis = basis_integrals(quad_rule, kernel, face, obs_loc, obs_n);
    return dot(basis, src_vals);
}


//TODO: test this
double appx_face_dist2(const Vec3<double>& pt,
                       const std::array<Vec3<double>,3>& vs) {
    double d0 = dist2(pt, vs[0]);
    double d1 = dist2(pt, vs[1]);
    double d2 = dist2(pt, vs[2]);
    return std::min(d0, std::min(d1, d2));
}


/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
//TODO: Use the diligenti mapping quadrature for the nearly 
//singular quadratures required.
template <typename T>
T richardson_step(const std::vector<T>& values) {
    assert(values.size() > 1);
    int n_steps = values.size();
    std::vector<T> last_level = values;
    std::vector<T> this_level;
    int error_order = 1;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, error_order);
            const double factor = 1.0 / (mult - 1.0);
            T low = last_level[i];
            T high = last_level[i + 1];
            T moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        error_order++;
        last_level = this_level;
    }
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
}

template double richardson_step(const std::vector<double>&);
template Vec3<double> richardson_step(const std::vector<Vec3<double>>&);

std::vector<double> integral_equation_vector(const Problem& p,
                                             const QuadStrategy& qs,
                                             const ObsPt& obs) {
    std::vector<double> result(p.src_mesh.vertices.size(), 0.0);
    for (unsigned int i = 0; i < p.src_mesh.faces.size(); i++) {
        FaceInfo src_face(p.src_mesh, i);
        const double dist2 = appx_face_dist2(obs.loc, src_face.corners);

        std::array<double,3> integrals;
        if (dist2 < pow(qs.far_threshold,2) * src_face.area) {
            std::vector<Vec3<double>> near_steps(qs.n_singular_steps, {0,0,0});
            auto near_quad = qs.get_near_quad(dist2 < src_face.area);
            for (int nf = 0; nf < qs.n_singular_steps; nf++) {

                double nfdn = 5 * obs.len_scale * qs.singular_steps[nf];
                auto nf_obs_pt = obs.loc + nfdn * obs.normal;

                near_steps[nf] += basis_integrals(*near_quad[nf], p.K,
                                           src_face, nf_obs_pt, obs.normal);
            }
            integrals = richardson_step(near_steps);
        } else {
            //farfield
            integrals = basis_integrals(qs.src_far_quad, p.K, src_face,
                                              obs.loc, obs.normal);
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

void outer_integral(std::vector<double>& integrals,
                    const std::array<int,3>& obs_face, double x_hat, double y_hat,
                    double jacobian, double eval, double q_wt) {
    for(int v = 0; v < 3; v++) {
        double obs_basis_eval = linear_interp(x_hat, y_hat, unit<double>(v)); 
#pragma omp critical
        {
            integrals[obs_face[v]] += jacobian * obs_basis_eval * eval * q_wt;
        }
    }
}

//TODO: Use a sparse matrix storage format here.
std::vector<std::vector<double>> interact_matrix(const Problem& p,
                                                 const QuadStrategy& qs) {
    int n_obs_basis = p.obs_mesh.vertices.size();
    int n_src_basis = p.src_mesh.vertices.size();
    std::vector<std::vector<double>> matrix(n_obs_basis,
            std::vector<double>(n_src_basis, 0.0));
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt::from_face(qs.obs_quad, obs_face, obs_q);

            const auto row = integral_equation_vector(p, qs, pt);


            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp(qs.obs_quad[obs_q].x_hat[0],
                                                      qs.obs_quad[obs_q].x_hat[1],
                                                      unit<double>(v)); 
                int b = obs_face.face[v];
#pragma omp critical
                for (int i = 0; i < n_src_basis; i++) {
                    matrix[b][i] += obs_basis_eval * row[i] * qs.obs_quad[obs_q].w * obs_face.jacobian;
                }
            }
        }
    }
    return matrix;
}



std::vector<double> direct_interact(const Problem& p,
                                    const QuadStrategy& qs) {

    std::vector<double> integrals(p.obs_mesh.vertices.size(), 0.0);
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt::from_face(qs.obs_quad, obs_face, obs_q);

            const double inner_integral = eval_integral_equation(p, qs, pt);

            outer_integral(integrals, obs_face.face, qs.obs_quad[obs_q].x_hat[0],
                           qs.obs_quad[obs_q].x_hat[1], obs_face.jacobian,
                           inner_integral, qs.obs_quad[obs_q].w);
        }
    }

    return integrals;
}

std::vector<double> mass_term(const Problem& p,
                              const QuadStrategy& qs) {
    std::vector<double> integrals(p.obs_mesh.vertices.size(), 0.0);
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            const double strength_interp = 
                linear_interp(qs.obs_quad[obs_q].x_hat[0], qs.obs_quad[obs_q].x_hat[1], index3(p.src_strength, obs_face.face));

            outer_integral(integrals, obs_face.face, qs.obs_quad[obs_q].x_hat[0], qs.obs_quad[obs_q].x_hat[1], 
                           obs_face.jacobian, strength_interp, qs.obs_quad[obs_q].w);
        }
    }
    return integrals;
}

double get_len_scale(Mesh& mesh, int which_face, int q) {
    auto face = mesh.faces[which_face];
    return std::sqrt(tri_area(index3(mesh.vertices, face))) / q;
}
