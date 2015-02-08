#include "UnitTest++.h"
#include "basis.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "constraint_builder.h"
#include "util.h"

using namespace tbem;

TEST(Interpolate) {
    auto m = sphere_mesh({0,0,0}, 1, 2);
    auto res = interpolate<3>(m, [](const Vec<double,3>& x) {return x[0];});
    for (unsigned int i = 0; i < m.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], m.facets[i][d][0], 1e-14);
        }
    }
}

TEST(ConstrainedInterpolate) {
    auto mesh = sphere_mesh({0,0,0}, 1, 2);

    auto continuity = mesh_continuity(mesh.begin());
    auto constraints = convert_to_constraints(continuity);
    auto constraint_matrix = from_constraints(constraints);

    auto res = constrained_interpolate<3>(
        mesh, 
        [](const Vec<double,3>& x) {return x[0];},
        constraint_matrix
    );
    for (unsigned int i = 0; i < mesh.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            CHECK_CLOSE(res[3 * i + d], mesh.facets[i][d][0], 1e-14);
        }
    }
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
