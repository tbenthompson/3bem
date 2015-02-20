#include "UnitTest++.h"
#include "vectorx.h"

using namespace tbem;


struct Data {
    BlockVectorX a{{1,2}, {3,4}};
    BlockVectorX b{{-1,-2}, {-3,-4}};
    BlockVectorX c{{0,0}, {0,0}};
};

TEST_FIXTURE(Data, VectorXSize) {
    CHECK_EQUAL(a.size(), 2);
    CHECK_EQUAL(a[0].size(), 2);
}

TEST_FIXTURE(Data, VectorXIdx) {
    CHECK_EQUAL(a[0][0], 1.0);
}

TEST_FIXTURE(Data, BlockVectorXAdd) {
    a += b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, BlockVectorXSub) {
    a -= -b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, BlockVectorXMul) {
    a *= 2;
    a += b * 2;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, MoreMath) {
    auto d = a * b + 2 - (-c);
    BlockVectorX correct{{1, -2}, {-7, -14}};
    CHECK_EQUAL(d, correct);
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
