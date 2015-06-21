#include "catch.hpp"
#include "dense_operator.h"

using namespace tbem;

TEST_CASE("add", "[dense]")
{
    DenseOperator A(2, 2, {1, 2, 3, 4});
    DenseOperator B(2, 2, {-1, -2, -3, -4});
    auto C = A.add(B);
    REQUIRE_ARRAY_EQUAL(C.data(), std::vector<double>(4, 0.0), 4);
}
