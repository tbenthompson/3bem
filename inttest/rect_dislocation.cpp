#include "3bem.h"
#include "elastic_kernels.h"
#include "disloc_shared.h"

using namespace tbem;

int main() {
    double surf_width = 4;
    int refine_surf = 3;
    double far_threshold = 3.0;
    int near_steps = 5;
    int src_quad_pts = 2;
    int obs_quad_pts = 2;
    double near_tol = 1e-3;

    QuadStrategy<3> qs(obs_quad_pts, src_quad_pts,
                    near_steps, far_threshold, near_tol);

    auto fault = rect_mesh(
        {-1, 0, -3.0}, {-1, 0, -0.0},
        {1, 0, -0.0}, {1, 0, -3.0}
    ).refine_repeatedly(refine_surf - 1);

    auto surface = rect_mesh(
        {-surf_width, -surf_width, 0}, {-surf_width, surf_width, 0},
        {surf_width, surf_width, 0}, {surf_width, -surf_width, 0}
    ).refine_repeatedly(refine_surf);

    ElasticHypersingular<3> hyp(30e9, 0.25);
    
    std::cout << "Number of surface DOFs: " << surface.n_dofs() << std::endl;

    VectorX dux(fault.n_dofs(), 1.0);
    VectorX duyz(fault.n_dofs(), 0.0);
    BlockVectorX du{dux, duyz, duyz};

    DislocationProblem<3> abc(hyp, surface, fault, qs, du);
    abc.solve();
}
