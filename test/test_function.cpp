#include "UnitTest++.h"
#include "function.h"

using namespace tbem;


struct Data {
    BlockFunction a{{1,2}, {3,4}};
    BlockFunction b{{-1,-2}, {-3,-4}};
    BlockFunction c{{0,0}, {0,0}};
};

TEST_FIXTURE(Data, FunctionSize) {
    CHECK_EQUAL(a.size(), 2);
    CHECK_EQUAL(a[0].size(), 2);
}

TEST_FIXTURE(Data, FunctionIdx) {
    CHECK_EQUAL(a[0][0], 1.0);
}

TEST_FIXTURE(Data, BlockFunctionAdd) {
    a += b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, BlockFunctionSub) {
    a -= -b;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, BlockFunctionMul) {
    a *= 2;
    a += b * 2;
    CHECK_EQUAL(a, c);
}

TEST_FIXTURE(Data, MoreMath) {
    auto d = a * b + 2 - (-c);
    BlockFunction correct{{1, -2}, {-7, -14}};
    CHECK_EQUAL(d, correct);
}


int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
