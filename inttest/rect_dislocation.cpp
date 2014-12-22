#include "bem.h"
#include "constraint.h"
#include "mesh_gen.h"
#include "quadrature.h"
#include "kernels.h"
#include "petsc_interface.h"
#include "util.h"

using namespace tbem;

int main() {
    //TODO: Clearly need a better way of handling vector problems!
    //-- doing ~9x as much work because the problem is not vectored/tensored

    double surf_width = 4;
    int refine_surf = 4;
    double far_threshold = 3.0;
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;
    double singular_tolerance = 1e-3;

    auto fault = rect_mesh(
        {-1, 0, -3.0}, {-1, 0, -0.0},
        {1, 0, -0.0}, {1, 0, -3.0}
    ).refine_repeatedly(refine_surf - 2);

    auto surface = rect_mesh(
        {-surf_width, -surf_width, 0}, {-surf_width, surf_width, 0},
        {surf_width, surf_width, 0}, {surf_width, -surf_width, 0}
    ).refine_repeatedly(refine_surf);

    auto raw_constraints =
        ConstraintMatrix::from_constraints(mesh_continuity<3>(surface));
    auto constraints = apply_discontinuities<3>(surface, fault, raw_constraints);
    std::cout << surface.facets.size() << std::endl;

    QuadStrategy<3> qs(obs_quad_pts, src_quad_pts,
                    near_steps, far_threshold, singular_tolerance);

    ElasticKernels<3> ek(30e9, 0.25);
    
    std::size_t n_fault_dofs = 3 * fault.facets.size();
    std::size_t n_surface_dofs = 3 * surface.facets.size();

    std::array<std::vector<double>,3> du = {
        std::vector<double>(n_fault_dofs, 1.0),   
        std::vector<double>(n_fault_dofs, 0.0),   
        std::vector<double>(n_fault_dofs, 0.0)
    };

    std::array<std::vector<double>,3> all_dofs_rhs = {
        std::vector<double>(n_surface_dofs, 0.0),   
        std::vector<double>(n_surface_dofs, 0.0),   
        std::vector<double>(n_surface_dofs, 0.0)
    };

    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            Problem<3> p = {fault, surface, ek.hypersingular_mat[k][j], du[j]};
            auto res = direct_interact(p, qs);
            for (unsigned int i = 0; i < res.size(); i++) {
                all_dofs_rhs[k][i] += res[i];
            }
        }
    }

    std::array<std::vector<double>,3> rhs = {
        constraints.get_reduced(all_dofs_rhs[0]),
        constraints.get_reduced(all_dofs_rhs[1]),
        constraints.get_reduced(all_dofs_rhs[2])
    };

    int n_reduced_surface_dofs = rhs[0].size();

    std::vector<double> vector_rhs(3 * n_reduced_surface_dofs);
    for (int d = 0; d < 3; d++) {
        std::copy(rhs[d].begin(), rhs[d].end(), vector_rhs.begin() + 
                                                d * n_reduced_surface_dofs);
    }

    std::array<std::array<std::vector<double>,3>,3> mats;
    TIC
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            Problem<3> p = {surface, surface, ek.hypersingular_mat[k][j], {}};
            mats[k][j] = interact_matrix(p, qs);
        }
    }
    TOC("Building matrices")

    int count = 0;
    auto surface_disp = solve_system(vector_rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            std::array<std::vector<double>,3> x_temp;
            for (int i = 0; i < 3; i++) {
                auto x_temp_reduced = std::vector<double>(n_surface_dofs);
                std::copy(x.begin() + i * n_reduced_surface_dofs,
                          x.begin() + (i + 1) * n_reduced_surface_dofs,
                          x_temp_reduced.begin());
                x_temp[i] = constraints.get_all(x_temp_reduced, n_surface_dofs);
            }

            std::array<std::vector<double>,3> y_temp;
            for (int k = 0; k < 3; k++) {
                auto y_temp_full = std::vector<double>(n_surface_dofs, 0.0);
                for (int j = 0; j < 3; j++) {
                    for (unsigned int mi = 0; mi < n_surface_dofs; mi++) {
                        for (unsigned int ni = 0; ni < n_surface_dofs; ni++) {
                            y_temp_full[mi] += 
                                -mats[k][j][mi * n_surface_dofs + ni] * x_temp[j][ni];
                        }
                    }
                }
                auto y_temp_reduced = constraints.get_reduced(y_temp_full);
                std::copy(y_temp_reduced.begin(), y_temp_reduced.end(),
                          y.begin() + k * n_reduced_surface_dofs);
            }
        });

    std::array<std::vector<double>,3> soln;
    for (int i = 0; i < 3; i++) {
        auto reduced_soln = std::vector<double>(n_reduced_surface_dofs);
        std::copy(surface_disp.begin() + i * n_reduced_surface_dofs,
                  surface_disp.begin() + (i + 1) * n_reduced_surface_dofs,
                  reduced_soln.begin());
        soln[i] = constraints.get_all(reduced_soln, n_surface_dofs);
    }
    hdf_out_surface<3>(
        "rect_dislocation.hdf5", surface, {soln[0], soln[1], soln[2]}
    ); 
}
