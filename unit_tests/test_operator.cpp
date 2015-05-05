#include "catch.hpp"
#include "dense_operator.h"
#include "block_operator.h"

using namespace tbem;

TEST_CASE("SimpleMatrixMultiply", "[operator]") {
    BlockDenseOperator matrix{
        1, 1, 
        {DenseOperator(2, 2, {0,1,1,0})}
    };
    auto res = matrix.apply({3, 4});
    REQUIRE(res[0] == 4);
    REQUIRE(res[1] == 3);
}

TEST_CASE("SimpleMatrixMultiplyWithComponents", "[operator]") {
    BlockDenseOperator matrix{
        2,2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {0})
        }
    };
    std::vector<double> x{3, 4};
    auto res = matrix.apply(x);
    REQUIRE(res[0] == 4);
    REQUIRE(res[1] == 3);
}

TEST_CASE("CombineComponents", "[operator]") {
    BlockDenseOperator matrix{
        2,2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {2}),
            DenseOperator(1, 1, {3})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[4] = {0,1,2,3};
    REQUIRE(combined_op.n_rows() == 2);
    REQUIRE(combined_op.n_cols() == 2);
    REQUIRE_ARRAY_EQUAL(combined_op.data(), correct, 4);
}

TEST_CASE("CombineComponents2By2", "[operator]") {
    BlockDenseOperator matrix{
        2, 2,
        {
            DenseOperator(2, 2, {0,1,2,3}),
            DenseOperator(2, 2, {4,5,6,7}),
            DenseOperator(2, 2, {8,9,10,11}),
            DenseOperator(2, 2, {12,13,14,15})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[16] = {
        0,1,4,5,
        2,3,6,7,
        8,9,12,13,
        10,11,14,15
    };
    REQUIRE(combined_op.n_rows() == 4);
    REQUIRE(combined_op.n_cols() == 4);
    REQUIRE_ARRAY_EQUAL(combined_op.data(), correct, 16);
}

TEST_CASE("CombineComponentsNonSquareBlocks", "[operator]") {
    BlockDenseOperator matrix{
        2, 2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 3, {1,2,3}),
            DenseOperator(3, 1, {4,8,12}),
            DenseOperator(3, 3, {5,6,7,9,10,11,13,14,15})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[16] = {
        0,1,2,3,
        4,5,6,7,
        8,9,10,11,
        12,13,14,15
    };
    REQUIRE(combined_op.n_rows() == 4);
    REQUIRE(combined_op.n_cols() == 4);
    REQUIRE_ARRAY_EQUAL(combined_op.data(), correct, 16);
}
