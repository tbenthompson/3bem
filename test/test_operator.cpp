#include "UnitTest++.h"
#include "operator.h"


using namespace tbem;

TEST(SimpleMatrixMultiply) {
    MatrixOperator<double,double,double> matrix{
        2,2,
        {0,1,1,0}
    };
    auto res = apply_operator(matrix, {3,4});
    CHECK_ARRAY_EQUAL(res, (std::vector<double>{4,3}), 2);
}

// TODO: Someday I want this test to compile. It'd be a really cool way to do
// block operators.
// TEST(BlockMatrixMultiplyCanWeRecursivelyUseMatrixOperator) {
//     MatrixOperator<
//         std::vector<double>,
//         std::vector<double>,
//         MatrixOperator<double,double,double>
//     > matrix{
//         2,2,
//         {
//             {1,1,{0}},
//             {1,1,{1}},
//             {1,1,{1}},
//             {1,1,{0}}
//         }
//     };
//     auto res = apply_operator(matrix, {{3},{4}});
//     CHECK_EQUAL(res[0][0], 4);
//     CHECK_EQUAL(res[1][0], 3);
// }

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
