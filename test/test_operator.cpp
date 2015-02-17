#include "UnitTest++.h"
#include "operator.h"
#include "test_shared.h"


using namespace tbem;

TEST(ReshapeToOperator) {
    Vec2<Vec2<double>> entry{{{0,1},{2,3}}};
    std::vector<Vec2<Vec2<double>>> A{entry};
    auto result = reshape_to_operator(1, 1, A);
    CHECK_EQUAL(result.ops.size(), 4);
    for (size_t i = 0; i < A.size(); i++) {
        for (int d1 = 0; d1 < 2; d1++) {
            for (int d2 = 0; d2 < 2; d2++) {
                CHECK_EQUAL(A[i][d1][d2], (*result.ops[d1 * 2 + d2].data)[i]);
            }
        }
    }
}


TEST(SimpleMatrixMultiply) {
    BlockOperator matrix{
        1, 1, 
        {make_operator(2, 2, {0,1,1,0})}
    };
    auto res = apply_operator(matrix, std::vector<double>{3,4});
    CHECK_EQUAL(res[0], 4);
    CHECK_EQUAL(res[1], 3);
}

TEST(SimpleMatrixMultiplyWithComponents) {
    BlockOperator matrix{
        2,2,
        {
            make_operator(1, 1, {0}),
            make_operator(1, 1, {1}),
            make_operator(1, 1, {1}),
            make_operator(1, 1, {0})
        }
    };
    auto res = apply_operator(matrix, std::vector<std::vector<double>>{{3},{4}});
    CHECK_EQUAL(res[0][0], 4);
    CHECK_EQUAL(res[1][0], 3);
}

TEST(CombineComponents) {
    BlockOperator matrix{
        2,2,
        {
            make_operator(1, 1, {0}),
            make_operator(1, 1, {1}),
            make_operator(1, 1, {2}),
            make_operator(1, 1, {3})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[4] = {0,1,2,3};
    CHECK_EQUAL(combined_op.ops[0].n_rows, 2);
    CHECK_EQUAL(combined_op.ops[0].n_cols, 2);
    CHECK_ARRAY_EQUAL(*combined_op.ops[0].data, correct, 4);
}

TEST(CombineComponents2By2) {
    BlockOperator matrix{
        2, 2,
        {
            make_operator(2, 2, {0,1,2,3}),
            make_operator(2, 2, {4,5,6,7}),
            make_operator(2, 2, {8,9,10,11}),
            make_operator(2, 2, {12,13,14,15})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[16] = {
        0,1,4,5,
        2,3,6,7,
        8,9,12,13,
        10,11,14,15
    };
    CHECK_EQUAL(combined_op.ops[0].n_rows, 4);
    CHECK_EQUAL(combined_op.ops[0].n_cols, 4);
    CHECK_ARRAY_EQUAL(*combined_op.ops[0].data, correct, 16);
}

TEST(CombineComponentsNonSquareBlocks) {
    BlockOperator matrix{
        2, 2,
        {
            make_operator(1, 1, {0}),
            make_operator(1, 3, {1,2,3}),
            make_operator(3, 1, {4,8,12}),
            make_operator(3, 3, {5,6,7,9,10,11,13,14,15})
        }
    };
    auto combined_op = combine_components(matrix);
    double correct[16] = {
        0,1,2,3,
        4,5,6,7,
        8,9,10,11,
        12,13,14,15
    };
    CHECK_EQUAL(combined_op.ops[0].n_rows, 4);
    CHECK_EQUAL(combined_op.ops[0].n_cols, 4);
    CHECK_ARRAY_EQUAL(*combined_op.ops[0].data, correct, 16);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
