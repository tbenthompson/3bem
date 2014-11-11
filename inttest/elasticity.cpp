#include "bem.h"
#include "mesh_gen.h"
#include "quadrature.h"
#include "kernels.h"
#include "petsc_interface.h"
#include "util.h"

int main() {
    //TODO: Problems:
    //-- doing ~9x as much work because the problem is not vectored
    //-- constraints are becoming a problem -- need to loosen constraints
    //   across the fault

    double surf_width = 12;
    int refine_surf = 6;
    double far_threshold = 3.0;
    int near_quad_pts = 3;
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;
    double singular_tolerance = 1e-2;

    auto fault = rect_mesh(
        {-1, 0, -3.0}, {-1, 0, -1.0},
        {1, 0, -1.0}, {1, 0, -3.0}
    );
    fault = refine_clean(fault, 2);

    auto surface = rect_mesh(
        {-surf_width, -surf_width, 0}, {-surf_width, surf_width, 0},
        {surf_width, surf_width, 0}, {surf_width, -surf_width, 0}
    );
    surface = refine_clean(surface, refine_surf);
    std::cout << surface.faces.size() << std::endl;

    QuadStrategy qs(obs_quad_pts, src_quad_pts, near_quad_pts,
                    near_steps, far_threshold, singular_tolerance);

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

    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            Problem p = {fault, surface, ek.hypersingular_mat[k][j], du[j]};
            auto res = direct_interact(p, qs);
            for (unsigned int i = 0; i < res.size(); i++) {
                rhs[k][i] += res[i];
            }
        }
    }

    std::vector<double> full_rhs(3 * surface.vertices.size());
    for (int d = 0; d < 3; d++) {
        std::copy(rhs[d].begin(), rhs[d].end(), full_rhs.begin() + d * n_surface_verts);
    }

    std::array<std::array<std::vector<std::vector<double>>,3>,3> mats;
    TIC
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            Problem p = {surface, surface, ek.hypersingular_mat[k][j], {}};
            mats[k][j] = interact_matrix(p, qs);
        }
    }
    TOC("Building matrices")

    int count = 0;
    auto surface_disp = solve_system(full_rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            std::array<std::vector<double>,3> x_temp;
            for (int i = 0; i < 3; i++) {
                x_temp[i] = std::vector<double>(n_surface_verts);
                std::copy(x.begin() + i * n_surface_verts,
                          x.begin() + (i + 1) * n_surface_verts,
                          x_temp[i].begin());
            }

            std::vector<double> y_temp(y.size(), 0.0);
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (unsigned int mi = 0; mi < mats[k][j].size(); mi++) {
                        for (unsigned int ni = 0; ni < mats[k][j].size(); ni++) {
                            y_temp[k * n_surface_verts + mi] += 
                                -mats[k][j][mi][ni] * x_temp[j][ni];
                        }
                    }
                }
            }
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });

    std::array<std::vector<double>,3> soln;
    for (int i = 0; i < 3; i++) {
        soln[i] = std::vector<double>(n_surface_verts);
        std::copy(surface_disp.begin() + i * n_surface_verts,
                  surface_disp.begin() + (i + 1) * n_surface_verts,
                  soln[i].begin());
    }
    hdf_out("strike_slip0.hdf5", surface, soln[0]); 
    hdf_out("strike_slip1.hdf5", surface, soln[1]); 
    hdf_out("strike_slip2.hdf5", surface, soln[2]); 
}
