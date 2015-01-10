#include "UnitTest++.h"
#include "output.h"
#include "mesh_gen.h"
#include "util.h"

using namespace tbem;

TEST(CreateHDFFile) {
    std::string filename = "build/test_hdf.hdf5";
    std::remove(filename.c_str());
    CHECK(!does_file_exist(filename));
    {
        int n_locs = 10;
        auto pts = reinterpret_vector<Vec2<double>>(random_list(2 * n_locs));
        auto data = random_list(n_locs);
        auto o = HDFOutputter(filename);
        out_volume<2>(o, pts, data, 1);
    }
    CHECK(does_file_exist(filename));
}

int main(int, char const* args[])
{
    return UnitTest::RunAllTests();
}
