#include "bem.h"
#include "numerics.h"
#include <iostream>
#include <algorithm>
#include <unordered_map>

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

NearEval::NearEval(int n_steps, int n_obs):
    n_steps(n_steps),
    steps(n_steps, std::vector<double>(n_obs)),
    quad(n_steps),
    dist(n_steps)
{
    for (int nf = 0; nf < n_steps; nf++) {
        quad[nf] = gauss((int)pow(2, nf + 1));
        dist[nf] = initial_dist / (pow(2, nf));
    }
}

void NearEval::zero_nears(int i) {
    for (int nf = 0; nf < n_steps; nf++) {
        steps[nf][i] = 0.0;
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

double appx_segment_distance(std::array<double, 2> v00,
                             std::array<double, 2> v01, 
                             std::array<double, 2> v10, 
                             std::array<double, 2> v11) {
    double d0 = dist2(v00, v10);
    double d1 = dist2(v00, v11);
    double d2 = dist2(v01, v10);
    double d3 = dist2(v01, v11);
    return std::min(d0, std::min(d1, std::min(d2, d3)));
}

// double 

double integral(QuadratureRule& quad_rule,
                KernelFnc& kernel,
                const std::array<double, 2>& src_v0,
                const std::array<double, 2>& src_v1,
                const double src_length,
                const double v0_val, 
                const double v1_val,
                double obs_x, 
                double obs_y,
                double dz) {
    double result = 0.0;
    for (unsigned int src_q = 0; src_q < quad_rule.size(); src_q++) {
        double x_hat = quad_rule[src_q].first;
        double src_x = ref_to_real(x_hat, src_v0[0], src_v1[0]);
        double src_y = ref_to_real(x_hat, src_v0[1], src_v1[1]);
        double quad_weight = quad_rule[src_q].second; 

        double interp_value = linear_interp(x_hat, v0_val, v1_val);

        double dx = obs_x - src_x;
        double dy = obs_y - src_y;
        double r = sqrt(dx * dx + dy * dy + dz * dz);

        double kernel_val = kernel(r, dx, dy);
        double jacobian = src_length / 2.0;
        result += quad_weight * interp_value * kernel_val * jacobian;
    }
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
    NearEval near_eval(n_steps, n_obs);

    for (int obs_idx = 0; obs_idx < n_obs_segs; obs_idx++) {
        auto obs_seg = obs_mesh.segments[obs_idx];
        auto obs_v0 = obs_mesh.vertices[obs_seg[0]];
        auto obs_v1 = obs_mesh.vertices[obs_seg[1]];
        double obs_length2 = dist2(obs_v0, obs_v1);
        double five_avg_quad_spacing = 5 * sqrt(obs_length2) / nq_obs;
        for (int obs_q = 0; obs_q < nq_obs; obs_q++) {
            double obs_x = ref_to_real(obs_quad[obs_q].first, obs_v0[0], obs_v1[0]);
            double obs_y = ref_to_real(obs_quad[obs_q].first, obs_v0[1], obs_v1[1]);
            int obs_value_idx = obs_idx * nq_obs + obs_q;
            obs_value[obs_value_idx] = 0.0;
            near_eval.zero_nears(obs_value_idx);

            for (auto src_seg: src_mesh.segments) {
                auto src_v0 = src_mesh.vertices[src_seg[0]];
                auto src_v1 = src_mesh.vertices[src_seg[1]];
                auto src_length = sqrt(dist2(src_v0, src_v1));
                double dist2 = appx_segment_distance(obs_v0, obs_v1,
                                                    src_v0, src_v1);

                auto v0_val = src_strength[src_seg[0]];
                auto v1_val = src_strength[src_seg[1]];
                if (dist2 < sqrt(3) * obs_length2) {
                    //nearfield
                    for (int nf = 0; nf < n_steps; nf++) {
                        double nfdz = five_avg_quad_spacing * near_eval.dist[nf];
                        double addme = integral(near_eval.quad[nf],
                                                kernel, src_v0, src_v1, 
                                                src_length, v0_val,
                                                v1_val, obs_x,
                                                obs_y, nfdz);
                        near_eval.steps[nf][obs_value_idx] += addme;
                    }
                } else {
                    //farfield
                     obs_value[obs_value_idx] += 
                         integral(src_quad, kernel, src_v0, src_v1, 
                                  src_length, v0_val, v1_val, obs_x,
                                  obs_y, 0.0);
                }
            }
            obs_value[obs_value_idx] += near_eval.steps[n_steps - 1][obs_value_idx];
        }
    }
    return obs_value;
}

//     int n_obs = obs.center.size();
//     int n_src = src.center.size();
//     std::cout << "Total interactions: " << (n_obs * n_src) << std::endl;
//     std::vector<double> obs_value(n_obs);
// 
// 
//     NearEval near_eval(n_steps, n_obs);
// 
//     int a = 0;
//     for (int i = 0; i < n_obs; i++) {
//         obs_value[i] = 0.0;
//         near_eval.zero_nears(i);
// 
//         for (int j = 0; j < n_src; j++) {
//             int which_ref = j % src.ref_center.size();
//             auto obs_loc = obs.center[i];
//             auto src_loc = src.center[j];
//             auto src_left = src.left[j];
//             auto src_right = src.right[j];
// 
// 
//             double dx = obs_loc[0] - src_loc[0];
//             double dy = obs_loc[1] - src_loc[1];
//             double r = sqrt(dx * dx + dy * dy);
// 
//             // this should be subsegment length, not segment length
//             double subseg_length = sqrt(pow(src_left[0] - src_right[0], 2) +
//                                         pow(src_left[1] - src_right[1], 2));
// 
//             auto seg = src_mesh.segments[src.owner[j]];
//             // auto v0 = src_mesh.vertices[seg[0]];
//             // auto v1 = src_mesh.vertices[seg[1]];
//             // double seg_length = sqrt(pow(v0[0] - v1[0], 2) +
//             //                             pow(v0[1] - v1[1], 2));
// 
//             double v0_val = src_strength[seg[0]];
//             double v1_val = src_strength[seg[1]];
//             double x_hat = src.ref_center[which_ref];
//             double interp_value = linear_interp(x_hat, v0_val, v1_val);
// 
//             // if near-field, do something more complex.
//             // std::cout << r << " " << seg_length << std::endl;
//             if (r < 0 * subseg_length) {
//                 //create a series of new result vectors that hold the steps of
//                 //the sequence converging to the correct result.
//                 //add to these result vectors the integral over this subseg
//                 //the integral over this subseg should be re-quadratured, 
//                 //maybe using a simple trapezoidal rule or something. since,
//                 //the point is nearby gauss does not have a big advantage 
//                 //anymore and trapezoidal is very simple.
//                 //
//                 //then perform the richardson extrapolation after the whole
//                 //summation loop is over
//                 //
//                 // maybe use 2 + nf as the number of quadrature points. ad hoc
//                 // but perhaps effective
//                 a++;
//                 for (int nf = 0; nf < n_steps; nf++) {
//                     double nfdz = r * near_eval.dist[nf];
//                     for (unsigned int qpi = 0;
//                          qpi < near_eval.quad[nf].size();
//                          qpi++) {
//                         auto qp = near_eval.quad[nf][qpi];
//                         double nfx_hat = ref_to_real(qp.first,
//                                                      src.ref_left[which_ref],
//                                                      src.ref_right[which_ref]);
//                         double nfx = ref_to_real(qp.first, src_left[0], src_right[0]);
//                         double nfy = ref_to_real(qp.first, src_left[1], src_right[1]);
//                         double nfdx = obs_loc[0] - nfx;
//                         double nfdy = obs_loc[1] - nfy;
//                         double nfr = sqrt(nfdx * nfdx + nfdy * nfdy + nfdz * nfdz);
//                         double interp_value = linear_interp(nfx_hat, v0_val, v1_val);
//                         double kernel_val = kernel(nfr, nfdx, nfdy);
//                         double jacobian = subseg_length / 2.0;
//                         near_eval.steps[nf][i] += 
//                             qp.second * kernel_val * interp_value * jacobian;
//                         // std::cout << which_ref << 
//                         //     " " << qp.second << 
//                         //     " " << kernel_val << 
//                         //     " " << interp_value << 
//                         //     " " << jacobian << std::endl;
//                     }
//                 }
//             } else {
//                 // if far-field, just compute the interaction.
//                 double kernel_val = kernel(r, dx, dy); 
//                 double jacobian = subseg_length;
//                 double quad_weight = src.ref_weight[which_ref];
//                 obs_value[i] += kernel_val * interp_value * jacobian;
//                 std::cout << which_ref << 
//                     " " << quad_weight << 
//                     " " << kernel_val << 
//                     " " << interp_value << 
//                     " " << jacobian << std::endl;
//             }
//         }
//         obs_value[i] += near_eval.steps[n_steps - 1][i];
//     }
//     std::cout << "Near-field interactions: " << a << std::endl;
//     return obs_value;
// }
