#include "UnitTest++.h"
#include "bem.h"

TEST(MeshToSubsegs) {
    std::vector<std::array<double, 2>> vertices = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
    };
    
    std::vector<std::array<int, 2>> segs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
}
