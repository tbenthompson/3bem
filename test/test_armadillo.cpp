#include "UnitTest++.h"
#include "armadillo_interface.h"
#include "operator.h"

using namespace tbem;

TEST(ArmadilloInvert) {
    std::vector<double> orig_mat{
        {-1, 1,
         2, 1}
    };

    auto arma_result = arma_invert({2, 2, orig_mat}).data;

    std::vector<double> correct{
        {-1 / 3.0, 1 / 3.0,
         2 / 3.0, 1 / 3.0}
    };

    CHECK_ARRAY_EQUAL(arma_result, correct, 4);
}

int main() {
    return UnitTest::RunAllTests();
}
