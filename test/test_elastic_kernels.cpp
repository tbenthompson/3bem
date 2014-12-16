#include "UnitTest++.h"
#include "kernels.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

template <int k, int j>
void test_elastic_kernel(std::string name) {
    std::ifstream test_data;
    std::string line;
    std::string filename = "test/test_data";
    test_data.open(filename, std::ios::in);
    CHECK(test_data.is_open()); 
    
    while (std::getline(test_data, line)) {
        std::vector<std::string> es;
        std::string entry;
        std::istringstream iss(line);
        while (std::getline(iss, entry,',')) {
            es.push_back(entry);
        }
        if (es[1] != name || std::stoi(es[2]) != k || std::stoi(es[3]) != j) {
            continue;
        }
        double exact = std::stod(es[0]);
        double shear_mod = std::stod(es[4]);
        double poisson_ratio = std::stod(es[5]);

        Vec3<double> src_loc = {std::stod(es[6]), std::stod(es[7]), std::stod(es[8])};
        Vec3<double> obs_loc = {std::stod(es[9]), std::stod(es[10]), std::stod(es[11])};
        Vec3<double> src_n = {std::stod(es[12]), std::stod(es[13]), std::stod(es[14])};
        Vec3<double> obs_n = {std::stod(es[15]), std::stod(es[16]), std::stod(es[17])};
        auto delta = (src_loc - obs_loc);
        double r2 = hypot2(delta);
        ElasticKernels<3> ek(shear_mod, poisson_ratio);
        double attempt = 0;
        if (name == "disp") {
             attempt = ek.displacement<k,j>(r2, delta, src_n, obs_n);
        } else if (name == "trac") {
             attempt = ek.traction<k,j>(r2, delta, src_n, obs_n);
        } else if (name == "adj_trac") {
             attempt = ek.adjoint_traction<k,j>(r2, delta, src_n, obs_n);
        } else if (name == "hyp") {
             attempt = ek.hypersingular<k,j>(r2, delta, src_n, obs_n);
        } 
        double error = std::fabs(exact - attempt) /
                       ((std::fabs(exact) + std::fabs(attempt)) / 2);
        CHECK_CLOSE(error, 0, 1e-8);
    }

    test_data.close();
}

template <int k, int j> 
void test_all() {
    test_elastic_kernel<k,j>("disp");
    test_elastic_kernel<k,j>("trac");
    test_elastic_kernel<k,j>("adj_trac");
    test_elastic_kernel<k,j>("hyp");
}

TEST(ABC) {
    test_all<0,0>(); test_all<0,1>(); test_all<0,2>();
    test_all<1,0>(); test_all<1,1>(); test_all<1,2>();
    test_all<2,0>(); test_all<2,1>(); test_all<2,2>();
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
