#include "UnitTest++.h"
#include "armadillo_facade.h"
#include "dense_operator.h"

using namespace tbem;

TEST(ArmadilloInvert) {
    std::vector<double> orig_mat{
        {-1, 1,
         2, 1}
    };

    auto arma_result = arma_invert(Operator(2, 2, orig_mat));

    std::vector<double> correct{
        {-1 / 3.0, 1 / 3.0,
         2 / 3.0, 1 / 3.0}
    };

    CHECK_ARRAY_EQUAL(&arma_result[0], correct, 4);
}

int main() {
    return UnitTest::RunAllTests();
}
