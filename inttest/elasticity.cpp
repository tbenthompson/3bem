#include "bem.h"
#include "mesh_gen.h"
#include "quadrature.h"
#include "kernels.h"
#include "petsc_interface.h"
#include "util.h"

int main() {
    //TODO: Super broken! Problems:
    //-- doing ~9x as much work because the problem is not vectored
    //-- constraints are becoming a problem -- need to loosen constraints
    //   across the fault
    Mesh fault = rect_mesh(
        {-1, 0, -2}, {-1, 0, 0},
        {1, 0, 0}, {1, 0, -2}
    );
    fault = refine_clean(fault, 1);

    double width = 4;
    Mesh surface = rect_mesh(
        {-width, -width, 0}, {-width, width, 0},
        {width, width, 0}, {width, -width, 0}
    );
    surface = refine_clean(surface, 5);

    double far_threshold = 2.0;
    int near_field = 4;
    auto q_src = tri_gauss(2);
    auto q_obs = tri_gauss(3);
    NearEval ne(near_field);

    ElasticKernels ek(30e9, 0.25);
    
    std::size_t n_fault_verts = fault.vertices.size();
    std::size_t n_surface_verts = surface.vertices.size();
    //TODO: Clearly need a better way of handling vector problems!
    std::array<std::vector<double>,3> du = {
        std::vector<double>(n_fault_verts, 1.0),   
        std::vector<double>(n_fault_verts, 0.0),   
        std::vector<double>(n_fault_verts, 0.0)
    };

    std::array<std::vector<double>,3> rhs = {
        std::vector<double>(n_surface_verts, 0.0),   
        std::vector<double>(n_surface_verts, 0.0),   
        std::vector<double>(n_surface_verts, 0.0)
    };

    using namespace std::placeholders;

    std::array<std::array<Kernel,3>,3> h = {{ {
            std::bind(&ElasticKernels::hypersingular<0,0>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<0,1>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<0,2>, ek, _1, _2, _3, _4)
        }, {
            std::bind(&ElasticKernels::hypersingular<1,0>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<1,1>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<1,2>, ek, _1, _2, _3, _4)
        }, {
            std::bind(&ElasticKernels::hypersingular<2,0>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<2,1>, ek, _1, _2, _3, _4),
            std::bind(&ElasticKernels::hypersingular<2,2>, ek, _1, _2, _3, _4)
        }
    }};

    // for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 3; j++) {
    //         auto res = direct_interact(fault, surface, q_src, q_obs,
    //                                    h[k][j], du[j], near_field, far_threshold);
    //         for (unsigned int i = 0; i < res.size(); i++) {
    //             rhs[k][i] += res[i];
    //         }
    //     }
    // }

    // std::vector<double> full_rhs(3 * surface.vertices.size());
    // for (int d = 0; d < 3; d++) {
    //     std::copy(rhs[d].begin(), rhs[d].end(), full_rhs.begin() + d * n_surface_verts);
    // }

    // int count = 0;
    // auto surface_disp = solve_system(full_rhs, 1e-5,
    //     [&] (std::vector<double>& x, std::vector<double>& y) {
    //         std::cout << "iteration " << count << std::endl;
    //         count++;
    //         std::array<std::vector<double>,3> x_temp;
    //         for (int i = 0; i < 3; i++) {
    //             x_temp[i] = std::vector<double>(n_surface_verts);
    //             std::copy(x.begin() + i * n_surface_verts,
    //                       x.begin() + (i + 1) * n_surface_verts,
    //                       x_temp[i].begin());
    //         }

    //         std::vector<double> y_temp(y.size(), 0.0);
    //         for (int k = 0; k < 3; k++) {
    //             for (int j = 0; j < 3; j++) {
    //                 auto res = direct_interact(fault, surface, q_src, q_obs,
    //                                            h[k][j], x_temp[j], near_field, 
    //                                            far_threshold);
    //                 for (unsigned int i = 0; i < res.size(); i++) {
    //                     y_temp[k * n_surface_verts + i] -= res[i];
    //                 }
    //             }
    //         }
    //         std::copy(y_temp.begin(), y_temp.end(), y.begin());
    //     });

    // std::array<std::vector<double>,3> soln;
    // for (int i = 0; i < 3; i++) {
    //     soln[i] = std::vector<double>(n_surface_verts);
    //     std::copy(surface_disp.begin() + i * n_surface_verts,
    //               surface_disp.begin() + (i + 1) * n_surface_verts,
    //               soln[i].begin());
    // }
    // hdf_out("strike_slip0.hdf5", surface, soln[0]); 
    // hdf_out("strike_slip1.hdf5", surface, soln[1]); 
    // hdf_out("strike_slip2.hdf5", surface, soln[2]); 
}
