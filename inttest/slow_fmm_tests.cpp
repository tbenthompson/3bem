#include "UnitTest++.h"
#include "octree.h"
#include "fmm.h"
#include "numerics.h"
#include "util.h"
#include "direct.h"

using namespace tbem;

void test_fmm_one(double mac2) {
    int n = 20;
    auto oct = random_pts_tree(n, 4);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(one_kernel, oct, strength, oct, 1, mac2);
    fmm_info.P2M();
    for (unsigned int i = 0; i < oct.cells.size(); i++) {
        auto c = oct.cells[i];
        CHECK_EQUAL(fmm_info.multipole_weights[i], c.end - c.begin);
    }
    fmm_info.fmm();
    fmm_info.L2P();
    std::vector<double> exact(n, n);
    CHECK_ARRAY_CLOSE(fmm_info.obs_effect, exact, n, 1e-14);
}

TEST(FMMOneOnlyP2P) {
    //Super big MAC forces P2P always
    test_fmm_one(1e30);
}

TEST(FMMOneOnlyM2L) {
    //Negative one MAC forces M2L always.
    test_fmm_one(-1.0);
}


TEST(FMMOneBothP2PM2L) {
    for (double mac2 = 1.0; mac2 < 50.0; mac2 += 5.0) {
        test_fmm_one(mac2);
    }
}

void test_fmm_laplace(double mac2, double error_max, bool non_one_strength) {
    int n = 2000;
    auto oct = random_pts_tree(n, 10);
    std::vector<double> strength(n, 1.0);
    if (non_one_strength) {
        strength = random_list(n);
    }
    FMMInfo fmm_info(laplace_single, oct, strength, oct, 3, mac2);
    fmm_info.P2M();
    fmm_info.fmm();
    fmm_info.L2P();
    auto exact = direct_n_body(oct.elements, oct.elements, laplace_single, strength);
    std::vector<double> error(n);
    for (unsigned int i = 0; i < error.size(); i++) {
        error[i] = std::fabs((fmm_info.obs_effect[i] - exact[i]) / exact[i]);
    }
    std::vector<double> zeros(n, 0.0);
    CHECK_ARRAY_CLOSE(error, zeros, n, error_max);
}

TEST(FMMLaplaceP2POnly) {
    //Super big MAC forces only P2P
    test_fmm_laplace(1e30, 1e-14, false);
}

TEST(FMMLaplace) {
    test_fmm_laplace(10.0, 1e-3, false);
}

TEST(FMMLaplaceNonOneStrength) {
    test_fmm_laplace(5.0, 3e-2, true);
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
