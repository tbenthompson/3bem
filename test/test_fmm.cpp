#include "UnitTest++.h"
#include "fmm.h"
#include "laplace_kernels.h"
#include "identity_kernels.h"
#include "util.h"
#include "dense_operator.h"

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

template <size_t dim>
NBodyData<dim> ones_data(size_t n) {
    auto src_pts = random_pts<dim>(n);
    auto obs_pts = random_pts<dim>(n);
    auto normals = random_pts<dim>(n);
    std::vector<double> weights(n, 1.0);
    return NBodyData<dim>{src_pts, normals, obs_pts, normals, weights};
}

TEST(SetupTreecode) {
    auto n = 200;
    LaplaceSingle<3> K;
    auto data = ones_data<3>(n);
    TreeNBodyOperator<3,1,1> tree(K, data, 5, 5, 5);
}

TEST(MakeSurroundingSurface) {
    auto surface = make_surrounding_surface<2>(4);
    CHECK_CLOSE(surface.pts[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.pts[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
    CHECK_CLOSE(surface.pts[3], (Vec<double,2>{0.0, -1.0}), 1e-12);
    CHECK_CLOSE(surface.normals[0], (Vec<double,2>{1.0, 0.0}), 1e-12);
    CHECK_CLOSE(surface.normals[1], (Vec<double,2>{0.0, 1.0}), 1e-12);
} 

TEST(P2M_2D) {
    auto n = 30000;
    size_t order = 15;
    LaplaceSingle<2> K;
    auto data = ones_data<2>(n);
    TreeNBodyOperator<2,1,1> tree(K, data, 100, order, 1.0);
    TreeNBodyOperator<2,1,1> tree2(K, data, n + 1, order, 2.0);
    BlockVectorX x({VectorX(data.src_weights)});
    TIC
    auto out = tree.apply(x);
    TOC("FMM");
    TIC2 
    auto out2 = tree2.apply(x);
    TOC("DIRECT");
    double sum = 0.0;
    // for (size_t i = 0; i < n; i++) {
    //     auto error = std::fabs((out2[0][i] - out[0][i]) / out2[0][i]);
    //     std::cout << error << std::endl;
    // }
}

// TEST(TreecodeIdentityScalar) {
//     int n = 6;
//     IdentityScalar<3> K;
//     auto src_pts = random_pts<3>(n);
//     auto obs_pts = random_pts<3>(n);
//     auto normals = random_pts<3>(n);
//     std::vector<double> weights(n, 1.0);
//     NBodyData<3> data{src_pts, normals, obs_pts, normals, weights};
//     TreeNBodyOperator<3,1,1> tree(K, data, 5, 5, 5);
//     VectorX out(data.obs_locs.size(), 1.0);
//     auto applied = tree.apply({out})[0];
//     CHECK_ARRAY_CLOSE(applied, (std::vector<double>(n, n)), n, 1e-14);
// }

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
