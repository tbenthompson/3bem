#include "UnitTest++.h"
#include "operator.h"


using namespace tbem;

TEST(ReshapeToOperator) {
    Vec2<Vec2<double>> entry{{{0,1},{2,3}}};
    std::vector<Vec2<Vec2<double>>> A{entry};
    auto result = reshape_to_operator(1, 1, A);
    CHECK_EQUAL(result.data.size(), 4);
    for (size_t i = 0; i < A.size(); i++) {
        for (int d1 = 0; d1 < 2; d1++) {
            for (int d2 = 0; d2 < 2; d2++) {
                CHECK_EQUAL(A[i][d1][d2], result.data[d1 * 2 + d2][i]);
            }
        }
    }
}


TEST(SimpleMatrixMultiply) {
    MatrixOperator matrix{
        2,2,1,1,
        {{0,1,1,0}}
    };
    auto res = apply_operator(matrix, std::vector<double>{3,4});
    CHECK_EQUAL(res[0], 4);
    CHECK_EQUAL(res[1], 3);
}

TEST(SimpleMatrixMultiplyWithComponents) {
    MatrixOperator matrix{
        1,1,2,2,
        {{0},{1},{1},{0}}
    };
    auto res = apply_operator(matrix, std::vector<std::vector<double>>{{3},{4}});
    CHECK_EQUAL(res[0][0], 4);
    CHECK_EQUAL(res[1][0], 3);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
