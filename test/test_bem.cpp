#include "UnitTest++.h"
#include "bem.h"
#include "numerics.h"
#include "test_shared.h"
#include <math.h>
#include <iostream>

Mesh square_mesh() {
    std::vector<std::array<double, 2>> vertices = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0},
    };
    
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };

    Mesh m = {vertices, segs};
    return m;
}

Mesh refined_square_mesh(int levels) {

    Mesh m = square_mesh();
    for (int i = 0; i < levels; i++) {
        m = refine_mesh(m, naturals(m.segments.size()));
    }
    return m;
}

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

TEST(HeavyInteract) {
    auto m = refined_square_mesh(3); 

    int n_verts = m.vertices.size();
    std::vector<double> src_str(n_verts);
    for (int i = 0; i < n_verts; i++) {
        src_str[i] = 1.0;
    }

    TIC
    auto obs_values = direct_interact(m, m, gauss(2), double_exp(2, 0.3),
                                      laplace_single, src_str, 6);
    TOC("direct_interact on " + std::to_string(m.segments.size()) + " segments");
    // for(auto o: obs_values) 
    // {
    //     std::cout << o << std::endl;
    // }
}
