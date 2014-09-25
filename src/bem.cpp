#include "bem.h"
#include "numerics.h"
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <iomanip>

void refine_edge(Mesh& new_mesh, std::array<int, 2> seg) {
    // Find the new vertex and segments.
    const auto v0 = new_mesh.vertices[seg[0]];
    const auto v1 = new_mesh.vertices[seg[1]];
    const std::array<double,2> midpt = {(v0[0] + v1[0]) / 2, (v0[1] + v1[1]) / 2};

    int next_vert = new_mesh.vertices.size();
    const std::array<int,2> new_seg0 = {seg[0], next_vert};
    const std::array<int,2> new_seg1 = {next_vert, seg[1]};

    new_mesh.vertices.push_back(midpt);
    new_mesh.segments.push_back(new_seg0);
    new_mesh.segments.push_back(new_seg1);
}

Mesh refine_mesh(const Mesh& m, std::vector<int> refine_these) {
    Mesh new_mesh;
    //Keep all the vertices (refinement should never remove a vertex).
    new_mesh.vertices = m.vertices;
    
    // Sort the refined edges so that we only have to check the
    // next one at any point in the loop.
    std::sort(refine_these.begin(), refine_these.end());

    // The next refined edge.
    int current = 0;

    for (int i = 0; i < (int)m.segments.size(); i++) {
        if (i == refine_these[current]) {
            current += 1;
            refine_edge(new_mesh, m.segments[i]);
        } else {
            new_mesh.segments.push_back(m.segments[i]);
        }
    }

    return new_mesh;
}

NearEval::NearEval(int n_steps):
    n_steps(n_steps),
    quad(n_steps),
    dist(n_steps)
{
    for (int nf = 0; nf < n_steps; nf++) {
        //TODO: need a much better distribution of points in the nearfield.
        //Gauss seems unlikely to be optimal.
        quad[nf] = gauss((int)pow(2, nf + 2));
        dist[nf] = initial_dist / (pow(2, nf + 1));
    }
}

inline double linear_interp(double x_hat, double v0_val, double v1_val) {
    return 0.5 * ((1 + x_hat) * v0_val + (1 - x_hat) * v1_val);
}



double dist2(std::array<double, 2> v0, std::array<double, 2> v1) {
    double dx = (v0[0] - v1[0]);
    double dy = (v0[1] - v1[1]);
    return dx * dx + dy * dy;
}

//TODO: test this
double appx_segment_distance(std::array<double, 2> pt,
                             std::array<double, 2> v0, 
                             std::array<double, 2> v1) {
    double d0 = dist2(pt, v0);
    double d1 = dist2(pt, v1);
    return std::sqrt(std::min(d0, d1));
}

double integral(QuadratureRule& quad_rule,
                KernelFnc& kernel,
                const std::array<double, 2>& src_v0,
                const std::array<double, 2>& src_v1,
                const double src_length,
                const double v0_val, 
                const double v1_val,
                const double obs_x, 
                const double obs_y) {
    double result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.size(); src_q++) {
        const double x_hat = quad_rule[src_q].first;
        //TODO: Generalize to use an arbitrary polynomial shape function.
        const double src_x = ref_to_real(x_hat, src_v0[0], src_v1[0]);
        const double src_y = ref_to_real(x_hat, src_v0[1], src_v1[1]);
        const double map_deriv_x = src_v1[0] - src_v0[0];
        const double map_deriv_y = src_v1[1] - src_v0[1];
        const double map_deriv_length = hypot(map_deriv_x, map_deriv_y);
        const double n_x = -map_deriv_y / map_deriv_length;
        const double n_y = map_deriv_x / map_deriv_length;
        const double quad_weight = quad_rule[src_q].second; 

        //TODO: Generalize to use an arbitrary polynomial shape function.
        const double interp_value = linear_interp(x_hat, v0_val, v1_val);

        //TODO: Is it src - obs or obs - src?
        const double dx = obs_x - src_x;
        const double dy = obs_y - src_y;
        const double r = hypot(dx, dy);

        //TODO: Need a way to specify which kernel parameters are necessary
        const double kernel_val = kernel(r, dx, dy, n_x, n_y);
        const double jacobian = src_length / 2.0;
        result += quad_weight * interp_value * kernel_val * jacobian;
    }
    return result;
}

/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
double richardson_step(std::vector<double>& values) {
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

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(Mesh& src_mesh,
                              QuadratureRule& src_quad,
                              NearEval& near_eval, 
                              std::array<double, 2> obs_pt,
                              std::array<double, 2> obs_normal,
                              KernelFnc& kernel,
                              std::vector<double>& src_strength) {
    double result = 0.0;
    std::vector<double> near_steps(near_eval.n_steps, 0.0);
    for (int i = 0; i < src_mesh.segments.size(); i++) {
        auto src_seg = src_mesh.segments[i];
        // std::cout << src_seg[0] << " " << src_seg[1] << std::endl;
        auto src_v0 = src_mesh.vertices[src_seg[0]];
        auto src_v1 = src_mesh.vertices[src_seg[1]];
        double src_length = sqrt(dist2(src_v0, src_v1));
        double dist = appx_segment_distance(obs_pt, src_v0, src_v1);

        double v0_val = src_strength[src_seg[0]];
        double v1_val = src_strength[src_seg[1]];
        if (dist < 6 * src_length) {
            //nearfield
            for (int nf = 0; nf < near_eval.n_steps; nf++) {
                double nfdz = 5 * (src_length / src_quad.size()) * near_eval.dist[nf];
                double obs_x = obs_pt[0] + nfdz * obs_normal[0]; 
                double obs_y = obs_pt[1] + nfdz * obs_normal[1]; 
                near_steps[nf] += 
                    integral(near_eval.quad[nf], kernel, src_v0, src_v1, 
                             src_length, v0_val, v1_val, obs_x,
                             obs_y);
            }
        } else {
            //farfield
            result += integral(src_quad, kernel, src_v0, src_v1, 
                          src_length, v0_val, v1_val, obs_pt[0],
                          obs_pt[1]);
        }
    }
    result += richardson_step(near_steps);
    return result;
}
                                           

std::vector<double> direct_interact(Mesh& src_mesh,
                                    Mesh& obs_mesh,
                                    QuadratureRule src_quad,
                                    QuadratureRule obs_quad,
                                    KernelFnc kernel,
                                    std::vector<double>& src_strength,
                                    int n_steps) {
    int n_obs_segs = obs_mesh.segments.size();
    int nq_obs = obs_quad.size();
    int n_obs = n_obs_segs * nq_obs;
    std::vector<double> obs_value(n_obs);
    NearEval near_eval(n_steps);

// #pragma omp parallel for
    for (int obs_idx = 0; obs_idx < n_obs_segs; obs_idx++) {
        auto obs_seg = obs_mesh.segments[obs_idx];
        auto obs_v0 = obs_mesh.vertices[obs_seg[0]];
        auto obs_v1 = obs_mesh.vertices[obs_seg[1]];
        for (int obs_q = 0; obs_q < nq_obs; obs_q++) {
            double obs_x_hat = obs_quad[obs_q].first;
            double obs_x = ref_to_real(obs_x_hat, obs_v0[0], obs_v1[0]);
            double obs_y = ref_to_real(obs_x_hat, obs_v0[1], obs_v1[1]);
            double map_deriv_x = obs_v1[0] - obs_v0[0];
            double map_deriv_y = obs_v1[1] - obs_v0[1];
            double map_deriv_length = hypot(map_deriv_x, map_deriv_y);
            double n_x = -map_deriv_y / map_deriv_length;
            double n_y = map_deriv_x / map_deriv_length;
            obs_value[obs_idx * nq_obs + obs_q] =
                eval_integral_equation(src_mesh, src_quad, near_eval,
                                       {obs_x, obs_y}, {n_x, n_y},
                                       kernel, src_strength);
        }
    }
    return obs_value;
}

