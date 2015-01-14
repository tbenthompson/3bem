#include "UnitTest++.h"
#include "function.h"

using namespace tbem;

TEST(FacetFunction) {
    FacetFunction<double,3> a = {1,2,3};
    CHECK_EQUAL(a.vertices, (Vec<double,3>{1,2,3}));
}

TEST(Function) {
    Function<double,3> mf{{{0.0, 1.0, 2.0}, {2.0, 3.0, 4.0}, {4.0, 5.0, 6.0}}};
    auto mf_refined = mf.refine();
    CHECK_EQUAL(mf_refined.facets[0].vertices[1], 0.5);
}

TEST(FunctionUnion) {
    Function<double,3> mf1{{{0.0, 1.0, 2.0}}};
    Function<double,3> mf2{{{2.0, 3.0, 4.0}, {4.0, 5.0, 6.0}}};
    auto mf_combined = Function<double,3>::form_union({mf1, mf2});
    CHECK_EQUAL(mf_combined.n_facets(), 3);
    CHECK_EQUAL(mf_combined.facets[1].vertices[2], 4.0);
}

TEST(CallRefineWithNoFacetsToRefine) {
    Function<double,3> mf1{{{0.0, 1.0, 2.0}}};
    mf1.refine({});
}


int main()
{
    return UnitTest::RunAllTests();
}
