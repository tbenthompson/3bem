#include "UnitTest++.h"
#include "bem.h"
#include "numerics.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "shared.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"

using namespace tbem;

TEST(IdentityTensor) {
    IdentityTensor<3,3,3> K;
    (void)K;
}

struct IntegrationProb {
    IntegrationProb():
        src_locs({{{0,0,0},{2,0,0},{0,1,0}}}),
        src_vals({1.0, 1.0, 1.0}),
        obs_loc({2.0, 2.0, 2.0}),
        obs_n({1.0, 0.0, 0.0}),
        q(tri_gauss(2)),
        face({src_locs[0], src_locs[1], src_locs[2]})
    { }
    
    template <typename KT>
    void go(const KT& k) {
        auto basis = integrate<Vec3<double>,2>(q, [&] (std::array<double,2> x_hat) {
                return eval_point_influence<3>(x_hat, k, FacetInfo<3>::build(face),
                                          obs_loc, obs_n);
            });
        result = dot_product(basis, src_vals);
    }

    void check() {
        CHECK_CLOSE(result, exact, 1e-5);
    }

    double result;
    double exact;
    std::vector<Vec3<double>> src_locs;
    Vec3<double> src_vals;
    Vec3<double> obs_loc;
    Vec3<double> obs_n;
    QuadRule<2> q;
    Facet<3> face;
};

TEST_FIXTURE(IntegrationProb, IntegralOne) {
    double abc = integrate<double,2>(q, [](std::array<double,2> x_hat){return 1.0;});
    CHECK_CLOSE(abc, 0.5, 1e-12);
    exact = 1.0; go(IdentityScalar<3>()); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceSingle) {
    exact = 0.0269063; go(LaplaceSingle<3>()); check();
}

TEST_FIXTURE(IntegrationProb, IntegralLaplaceDouble) {
    exact = -0.00621003; go(LaplaceDouble<3>()); check();
}

TEST(Richardson) {
    std::vector<double> input = {0.5, 0.3, 0.2, 0.15};
    double result = richardson_limit(input);
    CHECK_CLOSE(result, 0.1, 1e-12);
}

TEST(RichardsonZeros) {
    std::vector<double> input = {0.0, 0.0, 0.0, 0.0};
    double result = richardson_limit(input);
    CHECK_CLOSE(result, 0.0, 1e-12);
}

TEST_FIXTURE(IntegrationProb, RichardsonIntegral) {
    q = tri_gauss(3);
    double offset = 0.5;
    std::vector<double> vals;
    obs_n = {1,0,0};
    for (int i = 0; i < 5; i++) {
        obs_loc = {2.0, 2.0, 2.0 + offset};
        go(LaplaceSingle<3>());
        vals.push_back(result);
        offset /= 2;
    }
    double result = richardson_limit(vals);
    CHECK_CLOSE(result, 0.0269063, 1e-6);
}

struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(sphere_mesh(center, r, refine_level)),
        qs(gauss_order, gauss_order, near_eval, 2.0, 1e-2),
        obs_pt(random_pt()),
        obs_n(random_pt()),
        obs_length_scale(get_len_scale<3>(sphere, 0, gauss_order)),
        src_strength(std::vector<double>(sphere.n_dofs(), 1.0))
    {}

    template <typename KT>
    double go(const KT& k) {
        auto p = make_problem<3>(sphere, sphere, k);
        ObsPt<3> obs{obs_length_scale, obs_pt, obs_n, obs_n};
        auto op = integral_equation_vector(p, qs, obs);
        return bem_mat_mult(op, k, 1, src_strength)[0];
    }

    Mesh<3> sphere;
    QuadStrategy<3> qs;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    double obs_length_scale;
    std::vector<double> src_strength;
}; 

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(5, 3, 2);
    double result = ep.go(IdentityScalar<3>());
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-1);
}

TEST(ConstantLaplace) {
    EvalProb ep(5, 3, 2);
    double result = ep.go(LaplaceDouble<3>());
    CHECK_CLOSE(result, -1.0, 1e-3);
}

TEST(MatrixRowVsEval) {
    EvalProb ep(4, 3, 2);
    double result = ep.go(LaplaceDouble<3>());
    CHECK_CLOSE(result, -1.0, 1e-3);
}

TEST(MassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 0);
    std::vector<double> str(sphere.n_dofs(), 1.0);
    for (std::size_t i = 0; i < sphere.n_facets(); i++) {
        for (int d = 0; d < 3; d++) {
            if (sphere.facets[i].vertices[d][0] > 0.5) {
                str[3 * i + d] = 0.0;
            }
        }
    }
    auto p = make_problem<3>(sphere, sphere, IdentityScalar<3>());
    QuadStrategy<3> qs(2);
    auto res = mass_term(p, qs, str);
    CHECK_EQUAL(res.size(), sphere.n_dofs());
    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f.vertices);
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    CHECK_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}

TEST(FacetInfo2D) {
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{1.0, 0.0}};
    auto face_info = FacetInfo<2>::build(f);
    CHECK_EQUAL(&face_info.face, &f);
    CHECK_EQUAL(face_info.area_scale, 1);
    CHECK_EQUAL(face_info.length_scale, 1);
    CHECK_EQUAL(face_info.jacobian, 0.5);
    CHECK_EQUAL(face_info.normal, (Vec2<double>{0.0, 1.0}));
}

TEST(FacetInfo3D) {
    Facet<3> f{
        Vec3<double>{0.0, 0.0, 0.0},
        Vec3<double>{1.0, 0.0, 0.0},
        Vec3<double>{0.0, 1.0, 0.0}
    };
    auto face_info = FacetInfo<3>::build(f);
    CHECK_EQUAL(&face_info.face, &f);
    CHECK_EQUAL(face_info.area_scale, 0.5);
    CHECK_EQUAL(face_info.length_scale, std::sqrt(1.0 / 2.0));
    CHECK_EQUAL(face_info.jacobian, 1.0);
    CHECK_EQUAL(face_info.normal, (Vec3<double>{0.0, 0.0, 1.0}));
}

TEST(ObsPtFromFace) {
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{1.0, 1.0}};
    auto face_info = FacetInfo<2>::build(f);
    auto obs = ObsPt<2>::from_face({0}, face_info);
    CHECK_EQUAL(obs.len_scale, std::sqrt(2));
    CHECK_EQUAL(obs.loc, (Vec2<double>{0.5, 0.5}));
    CHECK_EQUAL(obs.normal, (Vec2<double>{-1.0 / std::sqrt(2), 1.0 / std::sqrt(2)}));
    CHECK_EQUAL(obs.richardson_dir, obs.normal);
}

TEST(TensorKernel) {
    ElasticDisplacement<2> k(30e9, 0.25);
    auto facet = Facet<2>{{{{-1.0, 0.0}, {1.0, 0.0}}}};
    auto facet_info = FacetInfo<2>::build(facet);
    auto result = eval_point_influence({0.0}, k, facet_info, {0.0, 1.0}, {0.0, 1.0});
    CHECK_CLOSE(result[0][1][1], 8.84194e-13, 1e-17);
    CHECK_CLOSE(result[1][1][1], 8.84194e-13, 1e-17);
    CHECK_EQUAL(result[0][0][0], 0.0); CHECK_EQUAL(result[0][0][1], 0.0);
    CHECK_EQUAL(result[0][1][0], 0.0); CHECK_EQUAL(result[1][0][0], 0.0);
    CHECK_EQUAL(result[1][0][1], 0.0); CHECK_EQUAL(result[1][1][0], 0.0);
}

TEST(TensorMassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 1);
    std::vector<Vec3<double>> str(sphere.n_dofs(), {1.0, 1.0, 1.0});
    auto p = make_problem<3>(sphere, sphere, IdentityTensor<3,3,3>());
    QuadStrategy<3> qs(2);
    auto res = mass_term(p, qs, str);
    CHECK_EQUAL(res[0].size(), 3);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
