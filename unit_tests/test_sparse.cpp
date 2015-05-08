#include "catch.hpp"
#include "sparse_operator.h"

using namespace tbem;

TEST_CASE("Empty sparse apply", "[sparse]") {
    BlockSparseOperator op(1, 1, 2, 2, {});
    auto out = op.apply({0.5, 4.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.0, 0.0}, 2);
}

TEST_CASE("Sparse apply", "[sparse]") {
    BlockSparseOperator op(1, 1, 2, 2, {
        {0, 0, 1.0}, {1, 0, 2.0}
    });
    auto out = op.apply({0.5, 0.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 1.0}, 2);
}

TEST_CASE("Block sparse apply", "[sparse]") {
    BlockSparseOperator op(2, 2, 1, 1, {
        {0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 4.0}
    });
    auto out = op.apply({0.5, 7.0});
    REQUIRE_ARRAY_EQUAL(out, std::vector<double>{0.5, 29.0}, 2);
}

