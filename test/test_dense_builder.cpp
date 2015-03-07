#include "UnitTest++.h"
#include "dense_builder.h"
#include "numerics.h"
#include "quadrature.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "util.h"
#include "test_shared.h"
#include "laplace_kernels.h"
#include "elastic_kernels.h"
#include "vectorx.h"

using namespace tbem;
struct EvalProb {
    EvalProb(int refine_level, int near_eval, int gauss_order,
             Vec3<double> center = Vec3<double>{0,0,0},
             double r = 3.0):
        sphere(sphere_mesh(center, r, refine_level)),
        qs(gauss_order, gauss_order, near_eval, 2.0, 1e-2),
        obs_pt(random_pt3d()),
        obs_n(random_pt3d()),
        obs_length_scale(std::sqrt(tri_area(sphere.facets[0]))),
        src_strength(std::vector<double>(sphere.n_dofs(), 1.0))
    {}

    double go(const Kernel<3,1,1>& k) {
        auto p = make_boundary_integral<3>(sphere, sphere, k);
        ObsPt<3> obs{obs_length_scale, obs_pt, obs_n, obs_n};
        auto op = mesh_to_points_operator({obs}, p.src_mesh, p.K, qs);
        return op.apply({src_strength})[0][0];
    }

    Mesh<3> sphere;
    QuadStrategy<3> qs;
    Vec3<double> obs_pt;
    Vec3<double> obs_n;
    double obs_length_scale;
    std::vector<double> src_strength;
}; 

TEST(CollectObsPts) {
    auto m = line_mesh({-1,0}, {3,0}).refine_repeatedly(1);    
    auto q = gauss(3);
    auto res = collect_obs_pts(m, q);
    CHECK_EQUAL(res.size(), m.facets.size() * q.size());
    for (size_t i = 0; i < q.size(); i++) {
        CHECK_CLOSE(res[i].loc[0], q[i].x_hat[0], 1e-15);
    }
}

TEST(EvalIntegralEquationSphereSurfaceArea) {
    EvalProb ep(5, 3, 2);
    IdentityScalar<3> identity;
    double result = ep.go(identity);
    double exact_surf_area = 4*M_PI*9;
    CHECK_CLOSE(result, exact_surf_area, 1e-1);
}

TEST(ConstantLaplace) {
    EvalProb ep(5, 3, 2);
    LaplaceDouble<3> double_kernel;
    double result = ep.go(double_kernel);
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
            if (sphere.facets[i][d][0] > 0.5) {
                str[3 * i + d] = 0.0;
            }
        }
    }
    IdentityScalar<3> identity;
    auto p = make_boundary_integral<3>(sphere, sphere, identity);
    QuadStrategy<3> qs(2);
    auto mass_op = mass_operator(p, qs);
    auto res = mass_op.apply({str})[0];
    CHECK_EQUAL(res.size(), sphere.n_dofs());
    double true_area = 0.0;
    for (auto f: sphere.facets) {
        true_area += tri_area(f);
    }
    double mass_area = 0.0;
    for (auto r: res) {
        mass_area += r;
    }   
    CHECK_CLOSE(mass_area, (5.0 / 6.0) * true_area, 1e-12);
}

TEST(FacetInfo2D) {
    Facet<2> f{Vec2<double>{0.0, 0.0}, Vec2<double>{3.0, 0.0}};
    auto face_info = FacetInfo<2>::build(f);
    CHECK_EQUAL(face_info.area_scale, 9);
    CHECK_EQUAL(face_info.length_scale, 3);
    CHECK_EQUAL(face_info.jacobian, 1.5);
    CHECK_EQUAL(face_info.normal, (Vec2<double>{0.0, 1.0}));
}

TEST(FacetInfo3D) {
    Facet<3> f{
        Vec3<double>{0.0, 0.0, 0.0},
        Vec3<double>{2.0, 0.0, 0.0},
        Vec3<double>{0.0, 2.0, 0.0}
    };
    auto face_info = FacetInfo<3>::build(f);
    CHECK_EQUAL(face_info.area_scale, 2.0);
    CHECK_EQUAL(face_info.length_scale, std::sqrt(2.0));
    CHECK_EQUAL(face_info.jacobian, 4.0);
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

TEST(TensorMassTerm) {
    auto sphere = sphere_mesh({0,0,0}, 1.0, 1);
    BlockVectorX str(3, VectorX(sphere.n_dofs(), 1.0));
    IdentityTensor<3,3,3> identity;
    auto p = make_boundary_integral<3>(sphere, sphere, identity);
    QuadStrategy<3> qs(2);
    auto mass_op = mass_operator(p, qs);
    CHECK_EQUAL(mass_op.ops.size(), 9);
    auto res = mass_op.apply({str});
    CHECK_EQUAL(res.size(), 3);
}

TEST(MakeProblem) {
    auto sphere1 = sphere_mesh({0,0,0}, 1.0, 0);
    auto sphere2 = sphere_mesh({0,0,0}, 1.0, 1);
    IdentityTensor<3,3,3> identity;
    auto p = make_boundary_integral<3>(sphere1, sphere2, identity);
    CHECK_EQUAL(&p.obs_mesh, &sphere1);
    CHECK_EQUAL(&p.src_mesh, &sphere2);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
