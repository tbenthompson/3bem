#include "UnitTest++.h"
#include "octree.h"
#include "fmm.h"
#include "numerics.h"
#include "test_shared.h"
#include "direct.h"

Octree simple_pts_tree(int n, int cell_max) {
    auto pts = random_pts(n);
    return Octree(pts, cell_max);
}

TEST(SetupFMM) {
    int n = 500;
    auto oct = simple_pts_tree(n, 10);
    auto strength = random_list(n);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 2);
    CHECK_EQUAL(fmm_info.multipole_weights.size(), oct.cells.size() * 8);
    CHECK_EQUAL(fmm_info.local_weights.size(), oct.cells.size() * 8);
    CHECK_EQUAL(fmm_info.n_exp_pts, 2);
    CHECK_CLOSE(fmm_info.nodes[0][0], 1.0 / sqrt(2), 1e-13);
    CHECK_CLOSE(fmm_info.nodes[1][2], -1.0 / sqrt(2), 1e-13);
    
    //Check that the octrees are referenced and not copied
    CHECK_EQUAL(&fmm_info.src_oct, &oct);
    CHECK_EQUAL(&fmm_info.obs_oct, &oct);
}

TEST(P2M_M2P_OneCell) {
    int n = 100;
    auto oct = simple_pts_tree(n, n+1);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 1);
    P2M_pts_cell(fmm_info, 0);
    CHECK_EQUAL(fmm_info.nodes.size(), 1);
    CHECK_CLOSE(fmm_info.multipole_weights[0], n, 1e-14);

    for(int i = 0; i < n; i++) {
        M2P_cell_pt(fmm_info, 0, i);
        CHECK_CLOSE(fmm_info.obs_effect[i], n, 1e-14);
    }
}

TEST(P2M_M2P_OneKernel) {
    int n = 8000;
    auto oct = simple_pts_tree(n, 10);
    std::vector<double> strength(n, 1.0);
    TIC
    FMMInfo fmm_info(laplace_single, oct, strength, oct, 1);
    P2M(fmm_info);
    for(unsigned int i = 0; i < oct.cells.size(); i++) {
        // CHECK_CLOSE(oct.cells[i].end - oct.cells[i].begin, 
        //             fmm_info.multipole_weights[i], 1e-14);
    }

    treecode_eval(fmm_info);
    TOC("Treecode");
    for (unsigned int i = 0; i < oct.elements.size(); i++) {
        // CHECK_CLOSE(fmm_info.obs_effect[i], n, 1e-14);
    }

    TIC2
    auto exact = direct_n_body(oct.elements, oct.elements, laplace_single, strength);
    TOC("Direct");

    std::vector<double> error(exact.size()); 
    for(unsigned int i = 0; i < exact.size(); i++) {
        error[i] = std::fabs((exact[i] - fmm_info.obs_effect[i]) / exact[i]);
    }
    std::vector<double> zeros(n, 0.0);
    CHECK_ARRAY_CLOSE(error, zeros, n, 1e-2);
}
