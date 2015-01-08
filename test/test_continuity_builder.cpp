#include "UnitTest++.h"
#include "continuity_builder.h"
#include "mesh_gen.h"

using namespace tbem;

TEST(ConstraintMesh) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(2);
    auto constraints = mesh_continuity<3>(sphere);
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    CHECK_EQUAL(3 * sphere.facets.size(), 384);
    CHECK_EQUAL(matrix.map.size(), 318);
    auto my_c = matrix.map.begin()->second.terms;
    CHECK_EQUAL(my_c[0].weight, 1);
}

TEST(GetReducedToCountTheNumberOfVerticesOnASphereApproximation) {
    auto sphere = sphere_mesh({0, 0, 0}, 1).refine_repeatedly(0);
    auto constraints = mesh_continuity<3>(sphere);
    auto matrix = ConstraintMatrix::from_constraints(constraints);
    std::vector<double> all(3 * sphere.facets.size(), 1.0);
    auto reduced = matrix.get_reduced(all);
    CHECK_EQUAL(reduced.size(), 6);
    CHECK_ARRAY_EQUAL(reduced, (std::vector<double>(6, 4.0)), 6);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
    // return RunOneTest("CircularConstraints");
}
