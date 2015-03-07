#include <unistd.h>
#include "UnitTest++.h"
#include "output.h"
#include "mesh_gen.h"
#include "util.h"

using namespace tbem;

TEST(CreateHDFFile) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {

    }
    std::string filepath = cwd;
    std::string filename = filepath + "/test_out/test_hdf.hdf5";
    std::remove(filename.c_str());
    CHECK(!does_file_exist(filename));
    {
        size_t n_locs = 10;
        std::vector<Vec2<double>> pts;
        for (size_t i = 0; i < n_locs; i++) {
            pts.push_back(random_pt2d());
        }
        auto data = random_list(n_locs);
        auto o = HDFOutputter(filename);
        out_volume<2>(o, pts, {data});
    }
    CHECK(does_file_exist(filename));
}

int main(int, char const* args[])
{
    return UnitTest::RunAllTests();
}
