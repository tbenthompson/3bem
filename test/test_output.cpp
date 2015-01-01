#include "UnitTest++.h"
#include "output.h"
#include "mesh_gen.h"
#include "util.h"

using namespace tbem;

class MockOutputter: public Outputter {
public:
    virtual void write_locations(int dim1, int dim2, const void* data_ptr) {
        loc_dim1 = dim1; 
        loc_dim2 = dim2; 
        loc_data_ptr = data_ptr;
    }

    virtual void write_values(int n_vars, const std::vector<double>& data) {
        vals_n_vars = n_vars;
        vals_data = data;
    }

    int loc_dim1;
    int loc_dim2;
    const void* loc_data_ptr;
    int vals_n_vars;
    std::vector<double> vals_data;
};

TEST(SurfaceOut) {
    auto m = circle_mesh({0,0}, 1).refine_repeatedly(0);
    int n_dofs = 2 * m.facets.size();
    auto data = random_list(n_dofs);
    MockOutputter mo;
    out_surface<2>(mo, m, data, 1);
    CHECK_EQUAL(mo.loc_dim1, m.facets.size());
    CHECK_EQUAL(mo.loc_dim2, 4);
    CHECK_ARRAY_EQUAL((double*)mo.loc_data_ptr, (double*)m.facets.data(),
                      4 * m.facets.size());
    CHECK_EQUAL(mo.vals_n_vars, 1);
    CHECK_ARRAY_EQUAL(mo.vals_data, data, n_dofs);
}

TEST(VolumeOut) {
    int n_locs = 10;
    auto pts = reinterpret_vector<Vec2<double>>(random_list(2 * n_locs));
    auto data = random_list(n_locs);
    MockOutputter mo;
    out_volume<2>(mo, pts, data, 1);
    CHECK_EQUAL(mo.loc_dim1, n_locs);
    CHECK_EQUAL(mo.loc_dim2, 2);
    CHECK_ARRAY_EQUAL((double*)mo.loc_data_ptr, (double*)pts.data(),
                      2 * n_locs);
    CHECK_EQUAL(mo.vals_n_vars, 1);
    CHECK_ARRAY_EQUAL(mo.vals_data, data, n_locs);
}

TEST(VolumeOutVec) {
    int n_locs = 10;
    auto pts = reinterpret_vector<Vec2<double>>(random_list(2 * n_locs));
    auto data = reinterpret_vector<Vec2<double>>(random_list(2 * n_locs));
    MockOutputter mo;
    out_volume<2>(mo, pts, data, 2);
    CHECK_EQUAL(mo.vals_n_vars, 2);
    CHECK_ARRAY_EQUAL(mo.vals_data, (double*)data.data(), 2 * n_locs);
}

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
