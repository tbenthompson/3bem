#include "UnitTest++.h"
#include "matrix_free_operator.h"
#include "petsc_facade.h"

using namespace tbem;

TEST(Create) {
    MatrixFreeOperator op(100, 100);
}

int main() {
    return UnitTest::RunAllTests();
}
