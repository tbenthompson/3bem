#include <cassert>
#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include "vec.h"
#include "quadrature.h"

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

double integral(const QuadratureRule2D& quad_rule,
                const KernelFnc& kernel,
                const std::array<Vec3<double>,3>& src_locs,
                const Vec3<double>& src_vals,
                const Vec3<double>& obs_loc,
                const Vec3<double>& obs_n) {

    //Compute normal and triangle area 
    //TODO: For linear elements, these could be done as a preprocessing step.
    //How to abstract this so that it work for both linear and high order basis?
    const auto unscaled_normal = tri_unscaled_normal(src_locs);
    const double src_area = tri_area(unscaled_normal);
    const double jacobian = src_area * 2.0;
    const auto scaled_normal = unscaled_normal / jacobian;

    double result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        const double x_hat = quad_rule.x_hat[src_q];
        const double y_hat = quad_rule.y_hat[src_q];
        const auto src_pt = ref_to_real(x_hat, y_hat, src_locs);
        const double interp_value = linear_interp(x_hat, y_hat, src_vals);

        const auto d = src_pt - obs_loc;

        const double r2 = hypot2(d);

        const double kernel_val = kernel(r2, d, scaled_normal, obs_n);
        const double q_wt = quad_rule.weights[src_q];

        result += q_wt * interp_value * kernel_val * jacobian;
    }
    return result;
}

/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
double richardson_step(const std::vector<double>& values) {
    assert(values.size() > 1);
    int n_steps = values.size();
    std::vector<double> last_level = values;
    std::vector<double> this_level;
    int error_order = 1;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, error_order);
            const double factor = 1.0 / (mult - 1.0);
            double low = last_level[i];
            double high = last_level[i + 1];
            double moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        error_order++;
        //TODO: Consider steps of two in error_order as an optional feature
        //TODO: Consider allowing setting the maximum error and then 
        // adaptively building
        //TODO: Integrate this into the main loop.
        //TODO: Use the diligenti mapping quadrature for the nearly 
        //singular quadratures required.
        last_level = this_level;
    }
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
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
    // double nfdn = 5 * (src_len_scale / src_quad.x_hat.size()) * 
    //               near_eval.dist[nf];
    double nfdn = length_factor * len_scale * ref_dist;

    // The new observation point moved a little bit off the
    // source surface.
    return obs_pt + nfdn * obs_normal;
}

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(const Mesh& src_mesh,
                              const QuadratureRule2D& src_quad,
                              const KernelFnc& kernel,
                              const NearEval& near_eval, 
                              const Vec3<double>& obs_pt,
                              const Vec3<double>& obs_normal,
                              double obs_len_scale,
                              std::vector<double>& src_strength,
                              const double far_threshold) {
    double result = 0.0;
    std::vector<double> near_steps(near_eval.n_steps, 0.0);
    for (unsigned int i = 0; i < src_mesh.faces.size(); i++) {
        auto src_face = src_mesh.faces[i];
        std::array<Vec3<double>,3> src_locs = {
            src_mesh.vertices[src_face[0]],
            src_mesh.vertices[src_face[1]],
            src_mesh.vertices[src_face[2]]
        };
        Vec3<double> src_vals = {
            src_strength[src_face[0]],
            src_strength[src_face[1]],
            src_strength[src_face[2]]
        };

        // Square of the threshold 
        // TODO: Make threshold a parameter
        const double threshold2 = far_threshold * far_threshold;

        // Square of the approximate source "length scale"
        const double src_area = tri_area(src_locs);

        // Square of the approximate distance from the source to observation.
        const double dist2 = appx_face_dist2(obs_pt, src_locs);

        //TODO: Better way of distinguishing nearfield and far-field.
        //a further hierarchy -- identical, adjacent, near, far
        if (dist2 < threshold2 * src_area) {
            //nearfield
            //cout
            
            for (int nf = 0; nf < near_eval.n_steps; nf++) {
                auto nf_obs_pt = near_field_point(near_eval.dist[nf], obs_pt,
                                                  obs_normal, obs_len_scale);
                near_steps[nf] += integral(near_eval.quad[nf], kernel,
                                           src_locs, src_vals, nf_obs_pt,
                                           obs_normal);
            }
        } else {
            //farfield
            double farfield_effect = integral(src_quad, kernel, src_locs,
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

std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule2D& src_quad,
                                    QuadratureRule2D& obs_quad,
                                    const KernelFnc& kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps,
                                    double far_threshold) {
    NearEval near_eval(n_steps);
    std::vector<double> integrals(obs_mesh.vertices.size(), 0.0);
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < obs_mesh.faces.size(); obs_idx++) {
        auto obs_face = obs_mesh.faces[obs_idx];
        auto obs_verts = index3(obs_mesh.vertices, obs_face);
        for (std::size_t obs_q = 0; obs_q < obs_quad.x_hat.size(); obs_q++) {
            double x_hat = obs_quad.x_hat[obs_q];
            double y_hat = obs_quad.y_hat[obs_q];
            double q_wt = obs_quad.weights[obs_q];

            const auto unscaled_normal = tri_unscaled_normal(obs_verts);
            const double obs_area = tri_area(unscaled_normal);
            double jacobian = obs_area * 2.0;
            const auto obs_n = unscaled_normal / jacobian;

            //TODO: What to use?
            const double obs_len_scale = std::sqrt(obs_area);

            // Observation point in real space
            const auto obs_pt = ref_to_real(x_hat, y_hat, obs_verts);

            const double inner_integral =
                eval_integral_equation(src_mesh, src_quad, kernel,
                                       near_eval, obs_pt, obs_n, obs_len_scale,
                                       src_strength, far_threshold);
            assert(!std::isnan(inner_integral));
            outer_integral(integrals, obs_face, x_hat, y_hat, 
                           jacobian, inner_integral, q_wt);
        }
    }

    return integrals;
}

std::vector<double> mass_term(const Mesh& obs_mesh,
                              const QuadratureRule2D& obs_quad,
                              const std::vector<double>& strengths) {

    std::vector<double> integrals(obs_mesh.vertices.size(), 0.0);

    for (std::size_t obs_idx = 0; obs_idx < obs_mesh.faces.size(); obs_idx++) {
        auto obs_face = obs_mesh.faces[obs_idx];
        auto obs_verts = index3(obs_mesh.vertices, obs_face);
        for (std::size_t obs_q = 0; obs_q < obs_quad.x_hat.size(); obs_q++) {
            double x_hat = obs_quad.x_hat[obs_q];
            double y_hat = obs_quad.y_hat[obs_q];
            double q_wt = obs_quad.weights[obs_q];

            const double strength_interp = linear_interp(x_hat, y_hat,
                                                         index3(strengths, obs_face));

            const double obs_area = tri_area(obs_verts);
            double jacobian = obs_area * 2.0;

            outer_integral(integrals, obs_face, x_hat, y_hat, 
                           jacobian, strength_interp, q_wt);
        }
    }
    return integrals;
}

double get_len_scale(Mesh& mesh, int which_face, int q) {
    auto face = mesh.faces[which_face];
    return std::sqrt(tri_area({
        mesh.vertices[face[0]],
        mesh.vertices[face[1]],
        mesh.vertices[face[2]]
    })) / q;
}
