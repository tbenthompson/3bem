#include "UnitTest++.h"
#include "function.h"
#include "vec_ops.h"

using namespace tbem;

TEST(Function) {
    Function<double,3> mf{{{0.0, 1.0, 2.0}, {2.0, 3.0, 4.0}, {4.0, 5.0, 6.0}}};
    auto mf_refined = mf.refine();
    CHECK_EQUAL(mf_refined.facets[0][1], 0.5);
}

TEST(FunctionUnion) {
    Function<double,3> mf1{{{0.0, 1.0, 2.0}}};
    Function<double,3> mf2{{{2.0, 3.0, 4.0}, {4.0, 5.0, 6.0}}};
    auto mf_combined = Function<double,3>::create_union({mf1, mf2});
    CHECK_EQUAL(mf_combined.n_facets(), 3);
    CHECK_EQUAL(mf_combined.facets[1][2], 4.0);
}

TEST(CallRefineWithNoFacetsToRefine) {
    Function<double,3> mf1{{{0.0, 1.0, 2.0}}};
    mf1.refine({});
}


int main()
{
    return UnitTest::RunAllTests();
}
