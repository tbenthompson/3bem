#include "catch.hpp"
#include "elastic_kernels.h"
#include "gravity_kernels.h"

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace tbem;


template <size_t dim>
Vec<double,dim> vec_from_indices(std::vector<std::string> row_entries,
                                 const size_t indices[]) {
    Vec<double,dim> out;
    for (size_t d = 0; d < dim; d++) {
        out[d] = std::stod(row_entries[indices[d]]);
    }
    return out;
}

template <size_t dim>
Vec<Vec<double,dim>,dim> call_kernel(const Vec<double,dim>& obs_loc,
                   const Vec<double,dim>& src_loc,
                   const Vec<double,dim>& obs_n,
                   const Vec<double,dim>& src_n,
                   const std::string& name,
                   double shear_modulus,
                   double poisson_ratio,
                   double b0, double b1) {
    if (name == "disp") {
         return ElasticDisplacement<dim>(shear_modulus, poisson_ratio)
            (obs_loc, src_loc, obs_n, src_n);
    } else if (name == "trac") {
         return ElasticTraction<dim>(shear_modulus, poisson_ratio)
            (obs_loc, src_loc, obs_n, src_n);
    } else if (name == "adj_trac") {
         return ElasticAdjointTraction<dim>(shear_modulus, poisson_ratio)
            (obs_loc, src_loc, obs_n, src_n);
    } else if (name == "hyp") {
         return ElasticHypersingular<dim>(shear_modulus, poisson_ratio)
            (obs_loc, src_loc, obs_n, src_n);
    } else if (name == "grav_disp") {
        return GravityDisplacement<dim>(shear_modulus, poisson_ratio, {b0, b1})
            (obs_loc, src_loc, obs_n, src_n);
    } else if (name == "grav_trac") {
        return GravityTraction<dim>(shear_modulus, poisson_ratio, {b0, b1})
            (obs_loc, src_loc, obs_n, src_n);
    }
    throw std::invalid_argument("name must be one of \
                                ['disp', 'trac', 'adj_trac', 'hyp']");
}

template <size_t dim>
void test_elastic_kernel(std::string name) {
    size_t src_loc_indices[3] = {8, 9, 10};
    size_t obs_loc_indices[3] = {11, 12, 13};
    size_t src_n_indices[3] = {14, 15, 16};
    size_t obs_n_indices[3] = {17, 18, 19};

    std::ifstream test_data;
    std::string line;
    std::string filename = "unit_tests/test_data_elastic" + std::to_string(dim);
    test_data.open(filename, std::ios::in);
    REQUIRE(test_data.is_open()); 
    
    while (std::getline(test_data, line)) {
        std::vector<std::string> es;
        std::string entry;
        std::istringstream iss(line);
        while (std::getline(iss, entry,',')) {
            es.push_back(entry);
        }

        size_t k = std::stoi(es[2]);
        size_t j = std::stoi(es[3]);
        if (es[1] != name) {
            continue;
        }
        double shear_mod = std::stod(es[4]);
        double poisson_ratio = std::stod(es[5]);
        double b0 = std::stod(es[6]);
        double b1 = std::stod(es[7]);

        auto src_loc = vec_from_indices<dim>(es, src_loc_indices);
        auto obs_loc = vec_from_indices<dim>(es, obs_loc_indices);
        auto src_n = vec_from_indices<dim>(es, src_n_indices);
        auto obs_n = vec_from_indices<dim>(es, obs_n_indices);

        double exact = std::stod(es[0]);
        double attempt = call_kernel<dim>(obs_loc, src_loc, obs_n, src_n,
                                          name, shear_mod, poisson_ratio, b0, b1)[k][j];
        double error = std::fabs(exact - attempt) /
                       ((std::fabs(exact) + std::fabs(attempt)) / 2);
        if (exact != 0.0) {
            REQUIRE_CLOSE(error, 0, 1e-8);
        } else {
            REQUIRE(attempt == 0.0);
        }
    }

    test_data.close();
}

TEST_CASE("TestElasticTensorKernels2DDisplacement", "[elastic_kernels]") {
    test_elastic_kernel<2>("disp");
}

TEST_CASE("TestElasticTensorKernels2DTraction", "[elastic_kernels]") {
    test_elastic_kernel<2>("trac");
}

TEST_CASE("TestElasticTensorKernels2DAdjointTraction", "[elastic_kernels]") {
    test_elastic_kernel<2>("adj_trac");
}

TEST_CASE("TestElasticTensorKernels2DHypersingular", "[elastic_kernels]") {
    test_elastic_kernel<2>("hyp");
}

TEST_CASE("TestElasticTensorKernels2DGravDisp", "[elastic_kernels]") {
    test_elastic_kernel<2>("grav_disp");
}

TEST_CASE("TestElasticTensorKernels2DGravTrac", "[elastic_kernels]") {
    test_elastic_kernel<2>("grav_trac");
}

TEST_CASE("TestElasticTensorKernels3DDisplacement", "[elastic_kernels]") {
    test_elastic_kernel<3>("disp");
}

TEST_CASE("TestElasticTensorKernels3DTraction", "[elastic_kernels]") {
    test_elastic_kernel<3>("trac");
}

TEST_CASE("TestElasticTensorKernels3DAdjointTraction", "[elastic_kernels]") {
    test_elastic_kernel<3>("adj_trac");
}

TEST_CASE("TestElasticTensorKernels3DHypersingular", "[elastic_kernels]") {
    test_elastic_kernel<3>("hyp");
}
