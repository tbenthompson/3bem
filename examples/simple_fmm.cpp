#include "octree.h"
#include "fmm.h"
#include "numerics.h"
#include "util.h"
#include "direct.h"

int main() {
    int n = (int)1e6;
    int pts_per_cell = 100;
    int n_exp_pts = 2;
    double mac2 = 4.0;
    auto pts = random_pts(n);
    std::vector<double> strength(n, 1.0);
    TIC;
    Octree oct(pts, pts_per_cell);
    TOC("Octree");
    TIC2;
    FMMInfo fmm_info(laplace_single, oct, strength, oct, n_exp_pts, mac2);
    TOC("Structures");
    TIC2;
    fmm_info.P2M();
    TOC("P2M");
    TIC2;
    fmm_info.fmm();
    TOC("M2L + P2P");
    std::cout << "Total possible billions of interactions: " << ((long)n * (long)n)/1e9 << std::endl;
    TIC2;
    fmm_info.L2P();
    TOC("L2P");

    TIC2;
    direct_n_body(oct.elements, oct.elements, laplace_single, strength);
    TOC("Direct");
}
