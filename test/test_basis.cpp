#include "UnitTest++.h"
#include "basis.h"
#include "mesh.h"
#include "mesh_gen.h"

using namespace tbem;

TEST(Interpolate) {
    auto m = sphere_mesh({0,0,0}, 1).refine_repeatedly(4);
    auto res = interpolate<3>(m, [](const Vec<double,3>& x) {return x[0];});
    for (unsigned int i = 0; i < m.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], m.facets[i].vertices[d][0], 1e-14);
        }
    }
}

TEST(ConstrainedInterpolate) {
    auto m = sphere_mesh({0,0,0}, 1).refine_repeatedly(4);
    auto constraints = ConstraintMatrix::from_constraints(mesh_continuity<3>(m));
    auto res = constrained_interpolate<3>(m, 
                [](const Vec<double,3>& x) {return x[0];}, constraints);
    for (unsigned int i = 0; i < m.facets.size(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], m.facets[i].vertices[d][0], 1e-14);
        }
    }
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
