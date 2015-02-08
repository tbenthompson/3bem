#include "UnitTest++.h"
#include "function.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

using namespace tbem;

struct Data {
    BlockFunction a{ {1,2}, {3,4} };
    BlockFunction b{ {-1,-2}, {-3,-4} };
    BlockFunction c{ {0,0}, {0,0} };
};

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

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
