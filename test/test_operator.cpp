#include "UnitTest++.h"
#include "dense_operator.h"
#include "vectorx.h"
#include "test_shared.h"

using namespace tbem;

TEST(SimpleMatrixMultiply) {
    BlockDenseOperator matrix{
        1, 1, 
        {DenseOperator(2, 2, {0,1,1,0})}
    };
    auto res = matrix.apply({{3,4}})[0];
    CHECK_EQUAL(res[0], 4);
    CHECK_EQUAL(res[1], 3);
}

TEST(SimpleMatrixMultiplyWithComponents) {
    BlockDenseOperator matrix{
        2,2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {0})
        }
    };
    BlockVectorX x({VectorX({3}), VectorX({4})});
    auto res = matrix.apply(x);
    CHECK_EQUAL(res[0][0], 4);
    CHECK_EQUAL(res[1][0], 3);
}

TEST(CombineComponents) {
    BlockDenseOperator matrix{
        2,2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 1, {1}),
            DenseOperator(1, 1, {2}),
            DenseOperator(1, 1, {3})
        }
    };
    auto combined_op = matrix.combine_components();
    double correct[4] = {0,1,2,3};
    CHECK_EQUAL(combined_op.n_rows(), 2);
    CHECK_EQUAL(combined_op.n_cols(), 2);
    CHECK_ARRAY_EQUAL(combined_op.data(), correct, 4);
}

TEST(CombineComponents2By2) {
    BlockDenseOperator matrix{
        2, 2,
        {
            DenseOperator(2, 2, {0,1,2,3}),
            DenseOperator(2, 2, {4,5,6,7}),
            DenseOperator(2, 2, {8,9,10,11}),
            DenseOperator(2, 2, {12,13,14,15})
        }
    };
    auto combined_op = matrix.combine_components();
    double correct[16] = {
        0,1,4,5,
        2,3,6,7,
        8,9,12,13,
        10,11,14,15
    };
    CHECK_EQUAL(combined_op.n_rows(), 4);
    CHECK_EQUAL(combined_op.n_cols(), 4);
    CHECK_ARRAY_EQUAL(combined_op.data(), correct, 16);
}

TEST(CombineComponentsNonSquareBlocks) {
    BlockDenseOperator matrix{
        2, 2,
        {
            DenseOperator(1, 1, {0}),
            DenseOperator(1, 3, {1,2,3}),
            DenseOperator(3, 1, {4,8,12}),
            DenseOperator(3, 3, {5,6,7,9,10,11,13,14,15})
        }
    };
    auto combined_op = matrix.combine_components();
    double correct[16] = {
        0,1,2,3,
        4,5,6,7,
        8,9,10,11,
        12,13,14,15
    };
    CHECK_EQUAL(combined_op.n_rows(), 4);
    CHECK_EQUAL(combined_op.n_cols(), 4);
    CHECK_ARRAY_EQUAL(combined_op.data(), correct, 16);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
