#include "UnitTest++.h"
#include "elastic_kernels.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace tbem;


template <int dim>
Vec<double,dim> vec_from_indices(std::vector<std::string> row_entries,
                                 const int indices[]) {
    Vec<double,dim> out;
    for (int d = 0; d < dim; d++) {
        out[d] = std::stod(row_entries[indices[d]]);
    }
    return out;
}

template <int dim>
Vec<Vec<double,dim>,dim> call_kernel(const double r2,
                   const Vec<double,dim>& delta, 
                   const Vec<double,dim>& src_n,
                   const Vec<double,dim>& obs_n,
                   const std::string& name,
                   double shear_modulus,
                   double poisson_ratio) {
    if (name == "disp") {
         return ElasticDisplacement<dim>(shear_modulus, poisson_ratio)
                                        (r2, delta, src_n, obs_n);
    } else if (name == "trac") {
         return ElasticTraction<dim>(shear_modulus, poisson_ratio)
                                    (r2, delta, src_n, obs_n);
    } else if (name == "adj_trac") {
         return ElasticAdjointTraction<dim>(shear_modulus, poisson_ratio)
                                           (r2, delta, src_n, obs_n);
    } else if (name == "hyp") {
         return ElasticHypersingular<dim>(shear_modulus, poisson_ratio)
                                           (r2, delta, src_n, obs_n);
    } 
    throw std::invalid_argument("name must be one of \
                                ['disp', 'trac', 'adj_trac', 'hyp']");
}

template <int dim>
void test_elastic_kernel(std::string name) {
    const int src_loc_indices[3] = {6, 7, 8};
    const int obs_loc_indices[3] = {9, 10, 11};
    const int src_n_indices[3] = {12, 13, 14};
    const int obs_n_indices[3] = {15, 16, 17};

    std::ifstream test_data;
    std::string line;
    std::string filename = "test/test_data_elastic" + std::to_string(dim);
    test_data.open(filename, std::ios::in);
    CHECK(test_data.is_open()); 
    
    while (std::getline(test_data, line)) {
        std::vector<std::string> es;
        std::string entry;
        std::istringstream iss(line);
        while (std::getline(iss, entry,',')) {
            es.push_back(entry);
        }

        int k = std::stoi(es[2]);
        int j = std::stoi(es[3]);
        if (es[1] != name) {
            continue;
        }
        double shear_mod = std::stod(es[4]);
        double poisson_ratio = std::stod(es[5]);

        auto src_loc = vec_from_indices<dim>(es, src_loc_indices);
        auto obs_loc = vec_from_indices<dim>(es, obs_loc_indices);
        auto src_n = vec_from_indices<dim>(es, src_n_indices);
        auto obs_n = vec_from_indices<dim>(es, obs_n_indices);
        auto delta = (src_loc - obs_loc);
        double r2 = hypot2(delta);

        double exact = std::stod(es[0]);
        double attempt = call_kernel<dim>(r2, delta, src_n, obs_n,
                                          name, shear_mod, poisson_ratio)[k][j];
        double error = std::fabs(exact - attempt) /
                       ((std::fabs(exact) + std::fabs(attempt)) / 2);
        // std::cout << "Name: " << name << " k,j,dim:" << k << j << dim << std::endl;
        // std::cout << attempt << " " << exact << std::endl;
        CHECK_CLOSE(error, 0, 1e-8);
    }

    test_data.close();
}

TEST(TestElasticTensorKernels2DDisplacement) {
    test_elastic_kernel<2>("disp");
}

TEST(TestElasticTensorKernels2DTraction) {
    test_elastic_kernel<2>("trac");
}

TEST(TestElasticTensorKernels2DAdjointTraction) {
    test_elastic_kernel<2>("adj_trac");
}

TEST(TestElasticTensorKernels2DHypersingular) {
    test_elastic_kernel<2>("hyp");
}

TEST(TestElasticTensorKernels3DDisplacement) {
    test_elastic_kernel<3>("disp");
}

TEST(TestElasticTensorKernels3DTraction) {
    test_elastic_kernel<3>("trac");
}

TEST(TestElasticTensorKernels3DAdjointTraction) {
    test_elastic_kernel<3>("adj_trac");
}

TEST(TestElasticTensorKernels3DHypersingular) {
    test_elastic_kernel<3>("hyp");
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
