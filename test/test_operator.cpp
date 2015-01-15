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
// 
// TEST(InnerProductTensorVector) {
//     Vec2<Vec2<double>> left{{{3.0, 0.0}, {0.0,4.0}}};
//     Vec2<double> right{1.0, 1.0};
//     Vec2<double> correct{{3.0, 4.0}};
//     auto result = apply_operator(left, right);
//     CHECK_EQUAL(result, correct);
// }
// 
// TEST(InnerProduct3TensorVector) {
//     Vec2<Vec2<Vec2<double>>> left{{
//         {{{3.0, 0.0}, {0.0, 4.0}}},
//         {{{3.0, 0.0}, {0.0, 4.0}}}
//     }};
//     Vec2<double> right = {1.0, 1.0};
//     Vec2<Vec2<double>> correct{{{6.0, 0.0}, {0.0,8.0}}};
//     auto result = apply_operator(left, right);
//     CHECK_EQUAL(result, correct);
// }

// TEST(InnerProductTensorTensor) {
//     Vec2<Vec2<double>> left{{{-1.0, 0.0}, {0.0,-2.0}}};
//     Vec2<Vec2<double>> right{{{3.0, 0.0}, {0.0,4.0}}};
// }

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
