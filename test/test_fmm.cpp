#include "UnitTest++.h"
#include "octree.h"
#include "fmm.h"
#include "numerics.h"
#include "util.h"
#include "direct.h"

using namespace tbem;

TEST(SetupFMM) {
    int n = 200;
    auto oct = random_pts_tree(n, 10);
    auto strength = random_list(n);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 2, 9.0);
    CHECK_EQUAL(fmm_info.multipole_weights.size(), oct.cells.size() * 8);
    CHECK_EQUAL(fmm_info.local_weights.size(), oct.cells.size() * 8);
    CHECK_EQUAL(fmm_info.obs_effect.size(), n);
    CHECK_EQUAL(fmm_info.n_exp_pts, 2);
    CHECK_CLOSE(fmm_info.nodes[0][0], 1.0 / sqrt(2), 1e-13);
    CHECK_CLOSE(fmm_info.nodes[1][2], -1.0 / sqrt(2), 1e-13);
    
    //Check that the octrees are referenced and not copied
    CHECK_EQUAL(&fmm_info.src_oct, &oct);
    CHECK_EQUAL(&fmm_info.obs_oct, &oct);
}

TEST(P2M_M2P_OneCell) {
    int n = 100;
    auto oct = random_pts_tree(n, n+1);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 1, 9.0);
    fmm_info.P2M_pts_cell(0);
    CHECK_EQUAL(fmm_info.nodes[0].size(), 1);
    CHECK_CLOSE(fmm_info.multipole_weights[0], n, 1e-14);

    for(int i = 0; i < n; i++) {
        auto cell = fmm_info.src_oct.cells[0];
        fmm_info.M2P_cell_pt(cell.bounds, 0, i);
        CHECK_CLOSE(fmm_info.obs_effect[i], n, 1e-14);
    }
}

TEST(TreecodeOneKernel) {
    int n = 6;
    auto oct = random_pts_tree(n, 5);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 1, 9.0);
    fmm_info.P2M();
    fmm_info.treecode();
    std::vector<double> exact(n, n);
    CHECK_ARRAY_CLOSE(fmm_info.obs_effect, exact, n, 1e-14);
}

TEST(TreecodeLaplace) {
    int n = 200;
    auto oct = random_pts_tree(n, 2);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(laplace_single, oct, strength, oct, 2, 20.0);
    fmm_info.P2M();
    fmm_info.treecode();
    auto exact = direct_n_body(oct.elements, oct.elements, laplace_single, strength);
    std::vector<double> error(n);
    for (unsigned int i = 0; i < error.size(); i++) {
        error[i] = std::fabs((fmm_info.obs_effect[i] - exact[i]) / exact[i]);
    }
    std::vector<double> zeros(n, 0.0);
    CHECK_ARRAY_CLOSE(error, zeros, n, 1e-2);
}

TEST(M2L_Kernel) {
    int n = 70;
    auto oct = random_pts_tree(n, 1);
    auto oct2 = random_pts_tree(n, 1);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(one_kernel, oct, strength, oct2, 1, 1.0);
    fmm_info.P2M();
    for(unsigned int i = 0; i < fmm_info.src_oct.cells.size(); i++) {
        auto src_cell =  fmm_info.src_oct.cells[i];
        for(unsigned int j = 0; j < fmm_info.obs_oct.cells.size(); j++) {
            fmm_info.local_weights[j] = 0.0;
            fmm_info.M2L_cell_cell(i, j);
            CHECK_EQUAL(fmm_info.local_weights[j], src_cell.end - src_cell.begin);
        }
    }
}

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
