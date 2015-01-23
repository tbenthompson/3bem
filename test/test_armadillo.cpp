#include "UnitTest++.h"
#include "armadillo_interface.h"

TEST(Blah) {
    std::vector<double> orig_mat{
        {-1, 1,
         2, 1}
    };

    auto arma_result = armadillo_invert(orig_mat);

    std::vector<double> correct{
        {-1 / 3.0, 1 / 3.0,
         2 / 3.0, 1 / 3.0}
    };

    CHECK_ARRAY_EQUAL(arma_result, correct, 4);
}

int main() {
    return UnitTest::RunAllTests();
}
