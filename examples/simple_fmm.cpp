#include "octree.h"
#include "fmm.h"
#include "numerics.h"
#include "test_shared.h"
#include "direct.h"

int main() {
    int n = (int)1e3;
    int pts_per_cell = 1000000;
    auto pts = random_pts(n);
    Octree oct(pts, pts_per_cell);
    std::vector<double> strength(n, 1.0);
    FMMInfo fmm_info(laplace_single, oct, strength, oct, 3, 10.0);
    TIC;
    fmm_info.P2M();
    TOC("P2M");
    TIC2;
    fmm_info.fmm();
    TOC("M2L + P2P");
    TIC2;
    fmm_info.L2P();
    TOC("L2P");
}
