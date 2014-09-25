#include "bem.h"
#include "test_shared.h"

int main() {
    auto m = refined_square_mesh(2); 

    int n_verts = m.vertices.size();
    std::vector<double> src_str(n_verts);
    for (int i = 0; i < n_verts; i++) {
        src_str[i] = 1.0;
    }

    TIC
    auto obs_values = direct_interact(m, m, gauss(2), double_exp(1, 0.3),
                                      laplace_single, src_str, 6);
    TOC("direct_interact on " + std::to_string(m.segments.size()) + " segments");
    // for(auto o: obs_values) 
    // {
    //     std::cout << o << std::endl;
    // }
}

