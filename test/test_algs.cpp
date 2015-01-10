#include "UnitTest++.h"
#include "algs.h"
#include "util.h"
#include "shared.h"

using namespace tbem;

TEST(MortonEncode) {
    int split = 2;
    uint64_t center = morton_encode(0, 0, split);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                uint64_t val = morton_encode(k, j, i); 
                if (i < split) {
                    CHECK(val < center);
                } else if (i == split && j == 0 && k == 0) {
                    CHECK(val == center); 
                } else {
                    CHECK(val > center);
                }
            }
        }
    }
}

TEST(MortonOrder) {
    uint64_t counter = 0;
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            for(int k = 0; k < 2; k++) {
                auto code = morton_encode(k, j, i);
                CHECK(code == counter);
                counter++;
            }
        }
    }
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
