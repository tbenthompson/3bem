#include "bem.h"
#include "util.h"

int main() {
    auto m = circle_mesh({0.0, 0.0}, 1.0, 8000); 
    /* m = refined_square_mesh(11); */

    int n_verts = m.vertices.size();
    std::vector<double> src_str(n_verts, 1.0);

    TIC
    auto obs_values = direct_interact(m, m, gauss(2), double_exp(1, 0.3),
                                      laplace_single, src_str, 6);
    TOC("direct_interact on " + std::to_string(m.segments.size()) + " segments");
}

