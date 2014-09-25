#include "UnitTest++.h"
#include "bem.h"
#include "numerics.h"
#include "test_shared.h"
#include <iostream>
#include <cmath>
#include <iomanip>

TEST(RefineMesh) {
    Mesh m = square_mesh();
    Mesh refined = refine_mesh(m, {0, 2});
    double correct[2][2] = {{0.5, 0.0}, {0.5, 1.0}};
    CHECK_ARRAY_CLOSE(refined.vertices[4], correct[0], 2, 1e-14);
    CHECK_ARRAY_CLOSE(refined.vertices[5], correct[1], 2, 1e-14);

    int new_segs[4][2] = {{0, 4}, {4, 1}, {2, 5}, {5, 3}};
    CHECK_ARRAY_EQUAL(refined.segments[0], new_segs[0], 2);
    CHECK_ARRAY_EQUAL(refined.segments[1], new_segs[1], 2);
    CHECK_ARRAY_EQUAL(refined.segments[3], new_segs[2], 2);
    CHECK_ARRAY_EQUAL(refined.segments[4], new_segs[3], 2);
}

TEST(HeavyLoadMeshRefine) {
    auto m = refined_square_mesh(17);
    
    TIC
    m = refine_mesh(m, naturals(m.segments.size()));
    TOC("refine_mesh to " + std::to_string(m.segments.size()) + " segments");
}

TEST(InteractOneKernel) {
    return;
    for (int refinement = 0; refinement < 5; refinement++) {
        for (int q_order = 0; q_order < 3; q_order++) {
            for (int n_near_steps = 1; n_near_steps < 4; n_near_steps++) {
                auto m = refined_square_mesh(refinement); 

                int n_verts = m.vertices.size();
                std::vector<double> src_str(n_verts);
                for (int i = 0; i < n_verts; i++) {
                    src_str[i] = 1.0;
                }

                auto obs_quad = double_exp(q_order, 0.3);
                auto src_quad = gauss(q_order + 1);
                auto obs_values = direct_interact(m, m, src_quad, obs_quad,
                                                  one, src_str, n_near_steps);
                for (auto o: obs_values) {
                    CHECK_CLOSE(o, 4.0, 1e-8);
                }
            }
        }
    }
}

TEST(CircleMesh) {
    int n_verts = 64;
    std::array<double, 2> center = {20.0, 0.0};
    Mesh src_circle = circle_mesh(center, 19.0, n_verts);
    CHECK_CLOSE(src_circle.vertices[16][0], 20.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[32][0], 1.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[48][0], 20.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[0][0], 39.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[16][1], 19.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[32][1], 0.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[48][1], -19.0, 1e-4);
    CHECK_CLOSE(src_circle.vertices[0][1], 0.0, 1e-4);
}

double exact_single(double obs_x, double obs_y) {
    return -0.0795775 * (
        -2 * obs_y * atan((-1 + obs_x) / obs_y) +
        2 * obs_y * atan((1 + obs_x) / obs_y) -
        (-1 + obs_x) * (-2 + log(pow((1 - obs_x), 2) + pow(obs_y, 2))) +
        (1 + obs_x) * (-2 + log(pow((1 + obs_x), 2) + pow(obs_y, 2))));
}

double exact_double(double x, double y) {
    return (atan((1 - x) / y) + atan((1 + x) / y)) / (2 * M_PI);
}

TEST(OneSegment) {
    std::array<double, 2> v0 = {-1.0, 0.0};
    std::array<double, 2> v1 = {1.0, 0.0};

    auto quad = gauss(15);
    std::vector<std::function<double (double, double)>> exact =
        {exact_single, exact_double};
    std::vector<KernelFnc> kernel = {laplace_single, laplace_double};
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                double obs_x = -5.0 + 10 * (i / 19.0);
                double obs_y = -5.0 + 10 * (j / 19.0);

                double result = integral(quad, kernel[k], v0, v1, 2.0, 1.0, 1.0,
                                         obs_x, obs_y);
                double exact_val = exact[k](obs_x, obs_y);
                CHECK_CLOSE(result, exact_val, 1e-4);
            }
        }
    }
}

TEST(DifferentApproaches) {
    //what is the value of the double layer potential approaching from
    //each direction?
    std::array<double, 2> v0 = {-1.0, 0.0};
    std::array<double, 2> v1 = {1.0, 0.0};
    auto quad = gauss(500);
    KernelFnc kernel = laplace_double;
    double up = integral(quad, kernel, v0, v1, 2.0, 1.0, 1.0,
                             0.0, 0.01);
    double down = integral(quad, kernel, v0, v1, 2.0, 1.0, 1.0,
                             0.0, -0.01);
    CHECK_CLOSE(up, 0.5, 1e-2);
    CHECK_CLOSE(down, -0.5, 1e-2);
    // std::cout << "UP: " << up << std::endl;
    // std::cout << "DOWN: " << down << std::endl;
}

TEST(ConstantLaplace) {
    int n_verts = 100;
    std::array<double, 2> center = {20.0, 0.0};
    Mesh src_circle = circle_mesh(center, 19.0, n_verts);
    KernelFnc double_layer = laplace_double;
    std::vector<double> u(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = 1.0;
    }
    for (double i = 1.0; i < 18.0; i++) {
        Mesh obs_circle = circle_mesh(center, i, n_verts);
        auto obs_values = direct_interact(src_circle, obs_circle, 
                                          gauss(2), double_exp(0, 0.3),
                                          laplace_double, u, 5);
        for (auto o: obs_values) {
            CHECK_CLOSE(o, 1.0, 1e-4);
        }
    }
}


double log_loc(std::array<double, 2> vert) {
    double x = vert[0];
    double y = vert[1];
    return log(hypot(x, y));
}

double log_loc_deriv(std::array<double, 2> vert, std::array<double, 2> dir) {
    double x = vert[0];
    double y = vert[1];
    double r = hypot(x, y);
    double dx = x / (r * r);
    double dy = y / (r * r);
    return dir[0] * dx + dir[1] * dy;
}

double theta_loc(std::array<double, 2> vert) {
    return atan(vert[1] / vert[0]);
}

double theta_loc_deriv(std::array<double, 2> vert, std::array<double, 2> dir) {
    double x = vert[0];
    double y = vert[1];
    double dy = 1.0 / (x * (1 + (y * y / (x * x))));
    double dx = (-y / x) * dy;
    return dir[0] * dx + dir[1] * dy;
}

TEST(LaplaceTest) {
    // f(z) = log(r) + i \theta 
    // is a harmonic function and thus both the real and imaginary parts
    // should satisfy the laplace equation.
    // By using the values of these functions on the boundary, I should be
    // able to use the eval_integral_equation function to evaluate the
    // correct values in the interior from the single and double layer 
    // potentials.
    int n_verts = 2000;
    std::array<double, 2> center = {20.0, 0.0};
    Mesh src_circle = circle_mesh(center, 19.0, n_verts);

    KernelFnc double_layer = laplace_double;
    KernelFnc single_layer = laplace_single;
    auto harmonic_fnc = log_loc;
    auto harmonic_fnc_deriv = log_loc_deriv;
    // auto harmonic_fnc = theta_loc;
    // auto harmonic_fnc_deriv = theta_loc_deriv;
    std::vector<double> u(n_verts);
    std::vector<double> du(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = harmonic_fnc(src_circle.vertices[i]);
        double normal_x = -src_circle.vertices[i][0] + center[0];
        double normal_y = -src_circle.vertices[i][1] + center[1];
        double normal_mag = hypot(normal_x, normal_y);
        normal_x /= normal_mag;
        normal_y /= normal_mag;
        du[i] = harmonic_fnc_deriv(src_circle.vertices[i], {normal_x, normal_y});
    }
    
    NearEval near_eval(5);
    auto quad = gauss(4);
    for (int i = 0; i < 20; i++) {
        double x = 1.0 + i * (38.0 / 19.0);
        double y = 10 * sin(2 * M_PI * (i / 19.0));
        double n_x = 1.0;
        if (x > 20) {
            n_x = -1.0;
        }
        double n_y = 0.0;
        double dbl_layer_value = 
            eval_integral_equation(src_circle, quad, near_eval, {x, y},
                                   {n_x, n_y}, double_layer, u);

        double sgl_layer_value = 
            eval_integral_equation(src_circle, quad, near_eval, {x, y},
                                   {n_x, n_y}, single_layer, du);
        double result = dbl_layer_value - sgl_layer_value;
        double correct = harmonic_fnc({x, y});
        //TODO: Think about PDE-wide Richardson extrapolation
        CHECK_CLOSE(result, correct, 2.0e-3);
        // std::cout << std::scientific << std::setprecision(15) << 
        //              std::fabs(correct - result) << " " << result << std::endl;
    }
}

//TODO: Write a test for richardson convergence.
