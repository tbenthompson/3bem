#include "UnitTest++.h"
#include "basis.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"

using namespace tbem;

TEST(Interpolate) {
    auto m = sphere_mesh({0,0,0}, 1).refine_repeatedly(2);
    auto res = interpolate<3>(m, [](const Vec<double,3>& x) {return x[0];});
    for (unsigned int i = 0; i < m.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], m.facets[i].vertices[d][0], 1e-14);
        }
    }
}

TEST(ConstrainedInterpolate) {
    auto mesh = sphere_mesh({0,0,0}, 1).refine_repeatedly(2);
    auto raw_constraints = mesh_continuity<3>(mesh);
    auto constraint_matrix = ConstraintMatrix::from_constraints(raw_constraints);
    auto res = constrained_interpolate<3>(
        mesh, 
        [](const Vec<double,3>& x) {return x[0];},
        constraint_matrix
    );
    for (unsigned int i = 0; i < mesh.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], mesh.facets[i].vertices[d][0], 1e-14);
        }
    }
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
