#include "UnitTest++.h"
#include "numbers.h"

using namespace tbem;

TEST(Range) {
    auto nats5 = range(-1, 8);
    double correct[] = {-1, 0, 1, 2, 3, 4, 5, 6, 7};
    CHECK_ARRAY_EQUAL(nats5, correct, 9);
}

TEST(Linspace) {
    auto vals = linspace(0.0, 1.0, 10);
    CHECK_EQUAL(vals.size(), 10);
    CHECK_EQUAL(vals[0], 0.0);
    CHECK_EQUAL(vals[vals.size() - 1], 1.0);
}

int main() {
    return UnitTest::RunAllTests();
}
