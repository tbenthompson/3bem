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

TEST_CASE("compose", "[dense]")
{
    DenseOperator A(2, 2, {1, 2, 3, 4});
    DenseOperator B(2, 2, {-1, -2, -3, -4});
    auto result = compose_dense_ops({A, B}, {0, 2}, {0, 2}, {1, 1}, 4, 4);
    REQUIRE(result.n_rows() == 4);
    REQUIRE(result.n_cols() == 4);
    REQUIRE(result[0] == 1);
    REQUIRE(result[1] == 2);
    REQUIRE(result[15] == -4);
}

TEST_CASE("compose overlap", "[dense]")
{
    DenseOperator A(2, 2, {1, 2, 3, 4});
    DenseOperator B(2, 2, {-1, -2, -3, -4});
    auto result = compose_dense_ops({A, B}, {0, 0}, {0, 0}, {1, -2}, 4, 4);
    REQUIRE(result.n_rows() == 4);
    REQUIRE(result.n_cols() == 4);
    REQUIRE(result[0] == 3);
    REQUIRE(result[1] == 6);
    REQUIRE(result[15] == 0);
}
