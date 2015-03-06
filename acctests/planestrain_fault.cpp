#include "disloc_shared.h"

using namespace tbem;

int main() {
    // HALF SPACE THRUST FAULT IN PLANE STRAIN

    // Fault mesh.
    auto fault = line_mesh({-1, -1}, {0, 0}).refine_repeatedly(0);

    // Simple low accuracy quadrature strategy.
    QuadStrategy<2> qs(3);

    // Earth's surface
    auto surface = line_mesh({-100, 0.0}, {100, 0.0}).refine_repeatedly(9);

    auto constraint_matrix = surf_fault_constraints(surface.begin(), fault.begin());

    double shear_mod = 30e9;
    double poisson = 0.25;
    ElasticHypersingular<2> hyp(shear_mod, poisson);

    double slip = 1;
    VectorX duxy(fault.n_dofs(), slip);
    BlockVectorX du{duxy, duxy};

    DislocationProblem<2> abc(hyp, surface, fault, qs, du);
    abc.solve();
}
