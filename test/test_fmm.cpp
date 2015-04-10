#include "UnitTest++.h"
#include "fmm.h"
#include "laplace_kernels.h"
#include "identity_kernels.h"
#include "util.h"

using namespace tbem;

TEST(LU) {
    std::vector<double> matrix{
        2, 1, -1, 0.5
    };
    auto lu = LU_decompose(matrix);
    auto soln = LU_solve(lu, {1,1});
    std::vector<double> correct{
        -0.25, 1.5
    };
    CHECK_ARRAY_CLOSE(soln, correct, 2, 1e-14);
}

TEST(SetupTreecode) {
    int n = 200;
    LaplaceSingle<3> K;
    auto src_pts = random_pts<3>(n);
    auto obs_pts = random_pts<3>(n);
    auto normals = random_pts<3>(n);
    std::vector<double> weights(n, 1.0);
    NBodyData<3> data{src_pts, normals, obs_pts, normals, weights};
    auto tree = make_tree_nbody_operator(K, data, 5, 5, 5);
}

TEST(MakeSurroundingSurface) {
    auto surface = make_surrounding_surface<2>(2);
    CHECK_CLOSE(surface.pts[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.pts[1], (Vec<double,2>{-1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.normals[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.normals[1], (Vec<double,2>{-1.0, 0.0}), 1e-12);
} 

TEST(TreecodeIdentityScalar) {
    int n = 6;
    IdentityScalar<3> K;
    auto src_pts = random_pts<3>(n);
    auto obs_pts = random_pts<3>(n);
    auto normals = random_pts<3>(n);
    std::vector<double> weights(n, 1.0);
    NBodyData<3> data{src_pts, normals, obs_pts, normals, weights};
    auto tree = make_tree_nbody_operator(K, data, 5, 5, 5);
    VectorX out(data.obs_locs.size(), 1.0);
    auto applied = tree.apply({out})[0];
    CHECK_ARRAY_CLOSE(applied, (std::vector<double>(n, n)), n, 1e-14);
}

// TEST(P2M_M2P_OneCell) {
//     int n = 100;
//     auto oct = random_pts_tree(n, n+1);
//     std::vector<double> strength(n, 1.0);
//     FMMInfo fmm_info(one_kernel, oct, strength, oct, 1, 9.0);
//     fmm_info.P2M_pts_cell(0);
//     CHECK_EQUAL(fmm_info.nodes[0].size(), 1);
//     CHECK_CLOSE(fmm_info.multipole_weights[0], n, 1e-14);
// 
//     for(int i = 0; i < n; i++) {
//         auto cell = fmm_info.src_oct.cells[0];
//         fmm_info.M2P_cell_pt(cell.bounds, 0, i);
//         CHECK_CLOSE(fmm_info.obs_effect[i], n, 1e-14);
//     }
// }

// 
// TEST(TreecodeLaplace) {
//     int n = 200;
//     auto oct = random_pts_tree(n, 2);
//     std::vector<double> strength(n, 1.0);
//     FMMInfo fmm_info(laplace_single, oct, strength, oct, 2, 20.0);
//     fmm_info.P2M();
//     fmm_info.treecode();
//     auto exact = direct_n_body(oct.elements, oct.elements, laplace_single, strength);
//     std::vector<double> error(n);
//     for (unsigned int i = 0; i < error.size(); i++) {
//         error[i] = std::fabs((fmm_info.obs_effect[i] - exact[i]) / exact[i]);
//     }
//     std::vector<double> zeros(n, 0.0);
//     CHECK_ARRAY_CLOSE(error, zeros, n, 1e-2);
// }
// 
// TEST(M2L_Kernel) {
//     int n = 70;
//     auto oct = random_pts_tree(n, 1);
//     auto oct2 = random_pts_tree(n, 1);
//     std::vector<double> strength(n, 1.0);
//     FMMInfo fmm_info(one_kernel, oct, strength, oct2, 1, 1.0);
//     fmm_info.P2M();
//     for(unsigned int i = 0; i < fmm_info.src_oct.cells.size(); i++) {
//         auto src_cell =  fmm_info.src_oct.cells[i];
//         for(unsigned int j = 0; j < fmm_info.obs_oct.cells.size(); j++) {
//             fmm_info.local_weights[j] = 0.0;
//             fmm_info.M2L_cell_cell(i, j);
//             CHECK_EQUAL(fmm_info.local_weights[j], src_cell.end - src_cell.begin);
//         }
//     }
// }

int main(int, char const *[])
{
    int retval = UnitTest::RunAllTests();
    return retval;
}
