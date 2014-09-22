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

Subsegments build_ref_subsegs(const QuadratureRule& quad) {
    Subsegments ref;
    ref.ref_left.resize(quad.size());
    ref.ref_center.resize(quad.size());
    ref.ref_weight.resize(quad.size());
    ref.ref_right.resize(quad.size());

    ref.ref_left[0] = -1.0;
    for (unsigned int i = 0; i < quad.size(); i++) {
        ref.ref_center[i] = quad[i].first;
        ref.ref_weight[i] = quad[i].second;
        // std::cout << quad[i].second << std::endl;
    }
    for (unsigned int i = 1; i < quad.size(); i++) {
        ref.ref_left[i] = (quad[i - 1].first + quad[i].first) / 2.0;
        ref.ref_right[i - 1] = ref.ref_left[i];
        /* std::cout << ref_subsegs.left[i] << std::endl; */
    }
    ref.ref_right[quad.size() - 1] = 1.0;

    return ref;
}

Subsegments get_src_obs(Mesh& m, const QuadratureRule& quad_rule) {
    auto subsegs = build_ref_subsegs(quad_rule);
    for (unsigned int i = 0; i < m.segments.size(); i++) {
        const auto v0 = m.vertices[m.segments[i][0]];
        const auto v1 = m.vertices[m.segments[i][1]];
        for (unsigned int k = 0; k < quad_rule.size(); k++) {
            const double lx = ref_to_real(subsegs.ref_left[k], v0[0], v1[0]);
            const double ly = ref_to_real(subsegs.ref_left[k], v0[1], v1[1]);
            const double cx = ref_to_real(subsegs.ref_center[k], v0[0], v1[0]);
            const double cy = ref_to_real(subsegs.ref_center[k], v0[1], v1[1]);
            const double rx = ref_to_real(subsegs.ref_right[k], v0[0], v1[0]);
            const double ry = ref_to_real(subsegs.ref_right[k], v0[1], v1[1]);
            subsegs.left.push_back({lx, ly});
            subsegs.center.push_back({cx, cy});
            subsegs.right.push_back({rx, ry});
            subsegs.owner.push_back(i);
        }
    }
    return subsegs;
}


NearEval::NearEval(int n_near_steps, int n_obs):
    n_near_steps(n_near_steps),
    near_steps(n_near_steps, std::vector<double>(n_obs)),
    near_quad(n_near_steps),
    near_dist(n_near_steps)
{
    for (int nf = 0; nf < n_near_steps; nf++) {
        near_quad[nf] = gauss((int)pow(2, nf + 1));
        near_dist[nf] = initial_dist / (pow(2, nf));
    }
}

void NearEval::zero_nears(int i) {
    for (int nf = 0; nf < n_near_steps; nf++) {
        near_steps[nf][i] = 0.0;
    }
}

inline double linear_interp(double x_hat, double v0_val, double v1_val) {
    return 0.5 * ((1 + x_hat) * v0_val + (1 - x_hat) * v1_val);
}

inline double laplace_single(double r, double dx, double dy) {
    return (1.0 / (4 * M_PI * r));
}

inline double one(double, double, double) {
    return 1.0;
}

std::vector<double> direct_interact(Mesh& src_mesh,
                                    Subsegments& src, 
                                    Subsegments& obs,
                                    std::vector<double> src_strength,
                                    int n_near_steps) {
    int n_obs = obs.center.size();
    int n_src = src.center.size();
    std::cout << "Total interactions: " << (n_obs * n_src) << std::endl;
    std::vector<double> obs_value(n_obs);

    // std::function<double (double, double, double)> kernel = laplace_single;
    std::function<double (double, double, double)> kernel = one;

    NearEval near_eval(n_near_steps, n_obs);

    int a = 0;
    for (int i = 0; i < n_obs; i++) {
        obs_value[i] = 0.0;
        near_eval.zero_nears(i);

        for (int j = 0; j < n_src; j++) {
            int which_ref = j % src.ref_center.size();
            auto obs_loc = obs.center[i];
            auto src_loc = src.center[j];
            auto src_left = src.left[j];
            auto src_right = src.right[j];


            double dx = obs_loc[0] - src_loc[0];
            double dy = obs_loc[1] - src_loc[1];
            double r = sqrt(dx * dx + dy * dy);

            // this should be subsegment length, not segment length
            double subseg_length = sqrt(pow(src_left[0] - src_right[0], 2) +
                                        pow(src_left[1] - src_right[1], 2));

            auto seg = src_mesh.segments[src.owner[j]];
            // auto v0 = src_mesh.vertices[seg[0]];
            // auto v1 = src_mesh.vertices[seg[1]];
            // double seg_length = sqrt(pow(v0[0] - v1[0], 2) +
            //                             pow(v0[1] - v1[1], 2));

            double v0_val = src_strength[seg[0]];
            double v1_val = src_strength[seg[1]];
            double x_hat = src.ref_center[which_ref];
            double interp_value = linear_interp(x_hat, v0_val, v1_val);

            // if near-field, do something more complex.
            // std::cout << r << " " << seg_length << std::endl;
            if (r < 0 * subseg_length) {
                //create a series of new result vectors that hold the steps of
                //the sequence converging to the correct result.
                //add to these result vectors the integral over this subseg
                //the integral over this subseg should be re-quadratured, 
                //maybe using a simple trapezoidal rule or something. since,
                //the point is nearby gauss does not have a big advantage 
                //anymore and trapezoidal is very simple.
                //
                //then perform the richardson extrapolation after the whole
                //summation loop is over
                //
                // maybe use 2 + nf as the number of quadrature points. ad hoc
                // but perhaps effective
                a++;
                for (int nf = 0; nf < n_near_steps; nf++) {
                    double nfdz = r * near_eval.near_dist[nf];
                    for (unsigned int qpi = 0;
                         qpi < near_eval.near_quad[nf].size();
                         qpi++) {
                        auto qp = near_eval.near_quad[nf][qpi];
                        double nfx_hat = ref_to_real(qp.first,
                                                     src.ref_left[which_ref],
                                                     src.ref_right[which_ref]);
                        double nfx = ref_to_real(qp.first, src_left[0], src_right[0]);
                        double nfy = ref_to_real(qp.first, src_left[1], src_right[1]);
                        double nfdx = obs_loc[0] - nfx;
                        double nfdy = obs_loc[1] - nfy;
                        double nfr = sqrt(nfdx * nfdx + nfdy * nfdy + nfdz * nfdz);
                        double interp_value = linear_interp(nfx_hat, v0_val, v1_val);
                        double kernel_val = kernel(nfr, nfdx, nfdy);
                        double jacobian = subseg_length / 2.0;
                        near_eval.near_steps[nf][i] += 
                            qp.second * kernel_val * interp_value * jacobian;
                        // std::cout << which_ref << 
                        //     " " << qp.second << 
                        //     " " << kernel_val << 
                        //     " " << interp_value << 
                        //     " " << jacobian << std::endl;
                    }
                }
            } else {
                // if far-field, just compute the interaction.
                double kernel_val = kernel(r, dx, dy); 
                double jacobian = subseg_length;
                double quad_weight = src.ref_weight[which_ref];
                obs_value[i] += kernel_val * interp_value * jacobian;
                std::cout << which_ref << 
                    " " << quad_weight << 
                    " " << kernel_val << 
                    " " << interp_value << 
                    " " << jacobian << std::endl;
            }
        }
        obs_value[i] += near_eval.near_steps[n_near_steps - 1][i];
    }
    std::cout << "Near-field interactions: " << a << std::endl;
    return obs_value;
}
