#include "UnitTest++.h"
#include "basis.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "continuity_builder.h"
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

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
