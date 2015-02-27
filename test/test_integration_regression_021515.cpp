#include "UnitTest++.h"
#include "integral_term.h"
#include "elastic_kernels.h"
#include "bem.h"

using namespace tbem;

TEST(IntegralDoesntHang) {
    UNITTEST_TIME_CONSTRAINT(5000);
    ElasticHypersingular<2> kernel(1.0, 0.25);
    QuadStrategy<2> qs(2, 2, 8, 3.0, 1e1);
    ObsPt<2> pt{
        2.2097086912079610e-03,
        {5.8620689655199998e-01, -5.8620689655199998e-01},
        {1.0, 0.0},
        {-7.0710678118654746e-01, -7.0710678118654746e-01}
    };
    FacetInfo<2> facet{
        {
          {{6.0937500000000000e-01, -6.0937500000000000e-01}, 
          {6.0156250000000000e-01, -6.0156250000000000e-01}}
        },
        1.2207031250000000e-04,
        1.1048543456039806e-02,
        5.5242717280199029e-03,
        {-7.0710678118654746e-01, -7.0710678118654746e-01}
    };
    auto term = make_integral_term(qs, kernel, pt, facet);
    compute_term(term);
}

int main() {
    UnitTest::RunAllTests();
}
