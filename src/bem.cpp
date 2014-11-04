#include <cassert>
#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include "vec.h"
#include "quadrature.h"
#include "taylor.h"

NearEval::NearEval(int n_steps):
    n_steps(n_steps),
    dist(n_steps)
{
    for (int nf = 0; nf < n_steps; nf++) {
        //TODO: need a much better distribution of points in the nearfield.
        //TODO: use the diligenti and aimi distribution per the 
        //nearly_singular_quad_test example
        quad.push_back(tri_gauss((int)pow(2, nf + 2)));
        dist[nf] = initial_dist / (pow(2, nf + 1));
    }
}

inline Vec3<double> ref_to_real(double x_hat, double y_hat, 
                                const std::array<Vec3<double>,3>& locs) {
    return {
        linear_interp(x_hat, y_hat, {locs[0][0], locs[1][0], locs[2][0]}),
        linear_interp(x_hat, y_hat, {locs[0][1], locs[1][1], locs[2][1]}),
        linear_interp(x_hat, y_hat, {locs[0][2], locs[1][2], locs[2][2]})
    };
}

FaceInfo::FaceInfo(const Mesh& mesh, int face_index):
    face(mesh.faces[face_index]),
    corners(index3(mesh.vertices, face))
{
    const auto unscaled_normal = tri_unscaled_normal(corners);
    area = tri_area(unscaled_normal);
    jacobian = area * 2.0;
    normal = unscaled_normal / jacobian;
}

struct SrcPointInfo {
    SrcPointInfo(const QuadratureRule2D& quad_rule,
              const Kernel& kernel,
              const FaceInfo& face,
              const Vec3<double>& obs_loc,
              const Vec3<double>& obs_n,
              int q_index) {
        const double x_hat = quad_rule.x_hat[q_index];
        const double y_hat = quad_rule.y_hat[q_index];

        const auto src_pt = ref_to_real(x_hat, y_hat, face.corners);
        const auto d = src_pt - obs_loc;
        const double r2 = hypot2(d);
        const double kernel_val = kernel(r2, d, face.normal, obs_n);

        const double q_wt = quad_rule.weights[q_index];
        weighted_kernel = kernel_val * q_wt * face.jacobian;

        basis = linear_basis(x_hat, y_hat);
    }
    
    Vec3<double> basis;
    double weighted_kernel;
};


std::array<double,3> basis_integrals(const QuadratureRule2D& quad_rule,
                                     const Kernel& kernel,
                                     const FaceInfo& face,
                                     const Vec3<double>& obs_loc,
                                     const Vec3<double>& obs_n) {

    Vec3<double> result = {0,0,0};
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        SrcPointInfo pt(quad_rule, kernel, face, obs_loc, obs_n, src_q);
        result += pt.weighted_kernel * pt.basis;
    }
    return result;
}


double integral(const QuadratureRule2D& quad_rule,
                const Kernel& kernel,
                const FaceInfo& face,
                const Vec3<double>& src_vals,
                const Vec3<double>& obs_loc,
                const Vec3<double>& obs_n) {

    double result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        SrcPointInfo pt(quad_rule, kernel, face, obs_loc, obs_n, src_q);
        result += pt.weighted_kernel * dot(pt.basis, src_vals);
    }
    return result;
}

//TODO: test this
double appx_face_dist2(const Vec3<double>& pt,
                       const std::array<Vec3<double>,3>& vs) {
    double d0 = dist2(pt, vs[0]);
    double d1 = dist2(pt, vs[1]);
    double d2 = dist2(pt, vs[2]);
    return std::min(d0, std::min(d1, d2));
}

Vec3<double> near_field_point(double ref_dist,
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double len_scale, 
                              double length_factor = 5.0) {
    double nfdn = length_factor * len_scale * ref_dist;

    // The new observation point moved a little bit off the
    // source surface.
    return obs_pt + nfdn * obs_normal;
}


std::vector<double> integral_equation_vector(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const Kernel& kernel,
                              const TaylorKernel& t_kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const double far_threshold) {
    std::vector<double> result(src_mesh.vertices.size(), 0.0);
    for (unsigned int i = 0; i < src_mesh.faces.size(); i++) {
        FaceInfo src_face(src_mesh, i);
        const double dist2 = appx_face_dist2(obs_pt, src_face.corners);

        //TODO: Better way of distinguishing nearfield and far-field.
        //a further hierarchy -- identical, adjacent, near, far
        std::array<double,3> integrals;
        if (dist2 < pow(far_threshold,2) * src_face.area) {
            std::vector<Vec3<double>> near_steps(near_eval.n_steps, {0,0,0});
            for (int nf = 0; nf < near_eval.n_steps; nf++) {
                auto nf_obs_pt = near_field_point(near_eval.dist[nf], obs_pt,
                                                  obs_normal, obs_len_scale);
                near_steps[nf] += basis_integrals(near_eval.quad[nf], kernel,
                                               src_face, nf_obs_pt,
                                               obs_normal);
            }
            integrals = richardson_step(near_steps);
        } else {
            //farfield
            integrals = basis_integrals(src_quad, kernel, src_face,
                                              obs_pt, obs_normal);
        }
        for (int b = 0; b < 3; b++) {
            result[src_face.face[b]] += integrals[b];
        }
    }
    return result;
}

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const Kernel& kernel,
                              const TaylorKernel& t_kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              const std::vector<double>& src_strength,
                              const double far_threshold) {
    double result = 0.0;
    std::vector<double> near_steps(near_eval.n_steps, 0.0);
    for (unsigned int i = 0; i < src_mesh.faces.size(); i++) {
        FaceInfo src_face(src_mesh, i);
        Vec3<double> src_vals = index3(src_strength, src_face.face);

        // Square of the approximate distance from the source to observation.
        const double dist2 = appx_face_dist2(obs_pt, src_face.corners);

        //TODO: Better way of distinguishing nearfield and far-field.
        //a further hierarchy -- identical, adjacent, near, far
        if (dist2 < pow(far_threshold,2) * src_face.area) {
            //nearfield
            
            // for (int nf = 0; nf < near_eval.n_steps; nf++) {
            //     auto nf_obs_pt = near_field_point(near_eval.dist[nf], obs_pt,
            //                                       obs_normal, obs_len_scale);
            //     near_steps[nf] += integral(near_eval.quad[nf], kernel,
            //                                src_face, src_vals, nf_obs_pt,
            //                                obs_normal);
            // }
            auto nf_obs_pt = near_field_point(near_eval.dist[nf], obs_pt,
                                              obs_normal, obs_len_scale);
            near_steps[nf] += integral(near_eval.quad[nf], kernel,
                                       src_face, src_vals, nf_obs_pt,
                                       obs_normal);
        } else {
            //farfield
            double farfield_effect = integral(src_quad, kernel, src_face,
                                              src_vals, obs_pt, obs_normal);
            result += farfield_effect;
        }
    }
    double nearfield_effect = richardson_step(near_steps);
    result += nearfield_effect;
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


struct ObsPointInfo {
    ObsPointInfo(const QuadratureRule2D& obs_quad,
                 const FaceInfo& obs_face,
                 int idx) {
        // Grab the quadrature point.
        x_hat = obs_quad.x_hat[idx];
        y_hat = obs_quad.y_hat[idx];
        q_wt = obs_quad.weights[idx];

        //TODO: What to use? Need to divide by q?
        len_scale = std::sqrt(obs_face.area);

        // Observation point in real space
        obs_pt = ref_to_real(x_hat, y_hat, obs_face.corners);
    }

    double x_hat;
    double y_hat;
    double q_wt;
    double len_scale;
    Vec3<double> obs_pt;
};

//TODO: Use a sparse matrix storage format here.
std::vector<std::vector<double>> interact_matrix(const Mesh& src_mesh,
                                    const Mesh& obs_mesh,
                                    const QuadratureRule2D& src_quad,
                                    const QuadratureRule2D& obs_quad,
                                    const Kernel& kernel,
                                    const TaylorKernel& t_kernel,
                                    int n_steps,
                                    double far_threshold) {
    NearEval near_eval(n_steps);
    int n_obs_basis = obs_mesh.vertices.size();
    int n_src_basis = src_mesh.vertices.size();
    std::vector<std::vector<double>> matrix(n_obs_basis,
            std::vector<double>(n_src_basis, 0.0));
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < obs_quad.x_hat.size(); obs_q++) {
            ObsPointInfo pt(obs_quad, obs_face, obs_q);

            const auto row =
                integral_equation_vector(src_mesh, src_quad, kernel, t_kernel, near_eval,
                                         pt.obs_pt, obs_face.normal, pt.len_scale,
                                         far_threshold);


            for(int v = 0; v < 3; v++) {
                double obs_basis_eval =
                        linear_interp(pt.x_hat, pt.y_hat, unit<double>(v)); 
                int b = obs_face.face[v];
#pragma omp critical
                for (int i = 0; i < n_src_basis; i++) {
                    {
                        matrix[b][i] += obs_basis_eval * row[i] * pt.q_wt * obs_face.jacobian;
                    }
                }
            }
        }
    }
    return matrix;
}



std::vector<double> direct_interact(const Mesh& src_mesh,
                                    const Mesh& obs_mesh,
                                    const QuadratureRule2D& src_quad,
                                    const QuadratureRule2D& obs_quad,
                                    const Kernel& kernel,
                                    const TaylorKernel& t_kernel,
                                    const std::vector<double>& src_strength,
                                    int n_steps,
                                    double far_threshold) {
    NearEval near_eval(n_steps);
    std::vector<double> integrals(obs_mesh.vertices.size(), 0.0);
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < obs_quad.x_hat.size(); obs_q++) {
            ObsPointInfo pt(obs_quad, obs_face, obs_q);

            const double inner_integral =
                eval_integral_equation(src_mesh, src_quad, kernel, t_kernel, near_eval,
                                       pt.obs_pt, obs_face.normal, pt.len_scale,
                                       src_strength, far_threshold);

            outer_integral(integrals, obs_face.face, pt.x_hat, pt.y_hat, 
                           obs_face.jacobian, inner_integral, pt.q_wt);
        }
    }

    return integrals;
}

std::vector<double> mass_term(const Mesh& obs_mesh,
                              const QuadratureRule2D& obs_quad,
                              const std::vector<double>& strengths) {

    std::vector<double> integrals(obs_mesh.vertices.size(), 0.0);
    for (std::size_t obs_idx = 0; obs_idx < obs_mesh.faces.size(); obs_idx++) {
        FaceInfo obs_face(obs_mesh, obs_idx);
        for (std::size_t obs_q = 0; obs_q < obs_quad.x_hat.size(); obs_q++) {
            ObsPointInfo pt(obs_quad, obs_face, obs_q);

            const double strength_interp = 
                linear_interp(pt.x_hat, pt.y_hat, index3(strengths, obs_face.face));

            outer_integral(integrals, obs_face.face, pt.x_hat, pt.y_hat, 
                           obs_face.jacobian, strength_interp, pt.q_wt);
        }
    }
    return integrals;
}

double get_len_scale(Mesh& mesh, int which_face, int q) {
    auto face = mesh.faces[which_face];
    return std::sqrt(tri_area(index3(mesh.vertices, face))) / q;
}
