#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include <iostream>
#include <cassert>

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

//TODO: Pass in observation normal.
double integral(QuadratureRule2D& quad_rule,
                KernelFnc& kernel,
                const std::array<std::array<double,3>,3>& src_locs,
                const std::array<double,3>& src_vals,
                const std::array<double,3>& obs_loc) {

    //Compute normal and triangle area 
    //TODO: For linear elements, these could be done as a preprocessing step.
    //How to abstract this so that it work for both linear and high order basis?
    const std::array<double,3> unscaled_normal = tri_unscaled_normal(src_locs);
    const double src_area = tri_area(unscaled_normal);
    const std::array<double,3> scaled_normal = {
        0.5 * unscaled_normal[0] / src_area,
        0.5 * unscaled_normal[1] / src_area,
        0.5 * unscaled_normal[2] / src_area
    };

    double result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.x_hat.size(); src_q++) {
        const double x_hat = quad_rule.x_hat[src_q];
        const double y_hat = quad_rule.y_hat[src_q];
        const double src_x = linear_interp(x_hat, y_hat,
                                    {src_locs[0][0], src_locs[1][0], src_locs[2][0]});
        const double src_y = linear_interp(x_hat, y_hat,
                                    {src_locs[0][1], src_locs[1][1], src_locs[2][1]});
        const double src_z = linear_interp(x_hat, y_hat,
                                    {src_locs[0][2], src_locs[1][2], src_locs[2][2]});
        const double interp_value = linear_interp(x_hat, y_hat, src_vals);

        const std::array<double,3> d = {
            obs_loc[0] - src_x,
            obs_loc[1] - src_y,
            obs_loc[2] - src_z
        };

        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

        const double kernel_val = kernel(r2, d, scaled_normal);
        const double q_wt = quad_rule.weights[src_q];

        const double jacobian = src_area * 2.0;
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
double appx_face_dist2(std::array<double,3> pt,
                       const std::array<std::array<double,3>,3> vs) {
    double d0 = dist2<3>(pt, vs[0]);
    double d1 = dist2<3>(pt, vs[1]);
    double d2 = dist2<3>(pt, vs[2]);
    return std::min(d0, std::min(d1, d2));
}

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(Mesh& src_mesh,
                              QuadratureRule2D& src_quad,
                              KernelFnc& kernel,
                              NearEval& near_eval, 
                              std::array<double,3> obs_pt,
                              std::array<double,3> obs_normal,
                              std::vector<double>& src_strength) {
    double result = 0.0;
    std::vector<double> near_steps(near_eval.n_steps, 0.0);
    for (unsigned int i = 0; i < src_mesh.faces.size(); i++) {
        auto src_face = src_mesh.faces[i];
        std::array<std::array<double,3>,3> src_locs = {
            src_mesh.vertices[src_face[0]],
            src_mesh.vertices[src_face[1]],
            src_mesh.vertices[src_face[2]]
        };
        std::array<double,3> src_vals = {
            src_strength[src_face[0]],
            src_strength[src_face[1]],
            src_strength[src_face[2]]
        };

        // Square of the approximate source "length scale"
        const double src_area = tri_area(src_locs);
        const double src_len_scale = std::sqrt(src_area);
        const double threshold = 36;

        // Square of the approximate distance from the source to observation.
        const double dist2 = appx_face_dist2(obs_pt, src_locs);

        //TODO: Better way of distinguishing nearfield and far-field.
        //a further hierarchy -- identical, adjacent, near, far
        if (dist2 < threshold * src_area) {
            //nearfield
            for (int nf = 0; nf < near_eval.n_steps; nf++) {
                double nfdn = 5 * (src_len_scale / src_quad.x_hat.size()) * 
                              near_eval.dist[nf];

                // The new observation point moved a little bit off the
                // source surface.
                std::array<double,3> nf_obs_pt = {
                    obs_pt[0] + nfdn * obs_normal[0],
                    obs_pt[1] + nfdn * obs_normal[1],
                    obs_pt[2] + nfdn * obs_normal[2]
                };

                near_steps[nf] += integral(near_eval.quad[nf], kernel,
                                           src_locs, src_vals, nf_obs_pt);
            }
        } else {
            //farfield
            result += integral(src_quad, kernel, src_locs, src_vals, obs_pt);
        }
    }
    result += richardson_step(near_steps);
    return result;
}
