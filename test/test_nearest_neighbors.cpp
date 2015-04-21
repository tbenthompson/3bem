#include "UnitTest++.h"
#include "nearest_neighbors.h"

using namespace tbem;

TEST(SimpleIdenticalPoints) 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = build_octree(es, 1);
    auto pts = identical_points({1.0, 2.0, 0.0}, es, oct);
    CHECK_EQUAL(pts.size(), 1);
    CHECK_EQUAL(pts[0], 0);
}

TEST(TwoIdenticalPoints) 
{
    std::vector<Vec<double,3>> es{
        {1.0, 2.0, 0.0}, {1.0, 2.0, 0.0}, {-1.0, 0.0, -3.0}, {0.0, -2.0, 3.0}
    };
    auto oct = build_octree(es, 1);
    auto pts = identical_points({1.0, 2.0, 0.0}, es, oct);
    CHECK_EQUAL(pts.size(), 2);
}

int main() {
    return UnitTest::RunAllTests();
}
