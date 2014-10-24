
#include "bem.h"
#include "numerics.h"
#include "mesh.h"
#include "util.h"
#include "petsc_interface.h"
#include <iostream>

//TODO: This could be the start of a Vec3 class. 
class Pt {
public:
    std::array<double,3> x;
    friend std::ostream& operator<<(std::ostream& os, const Pt& obj)
    {
        os << "(" << obj.x[0] << ", " << obj.x[1] << ", " << obj.x[2] << ")";
        return os;
    }
};

double harmonic_u(std::array<double,3> x) {
    return 1.0 / hypot(x);
}

double harmonic_dudn(std::array<double,3> x,
                     std::array<double,3> center) {
    double nx = x[0] - center[0];
    double ny = x[1] - center[1];
    double nz = x[2] - center[2];
    double n_mag = hypot(nx, ny, nz);
    nx /= n_mag; ny /= n_mag; nz /= n_mag;
    double r = hypot(x);
    return -(x[0] * nx + x[1] * ny + x[2] * nz) / (r * r * r);
}

//TODO: This should be in util.h i think
double error(std::vector<double> a, 
             std::vector<double> b) {
    double error = 0.0;
    for (unsigned int i = 0; i < a.size(); i++) {
        error += std::fabs((a[i] - b[i]) / b[i]) / b.size();
    }
    return error;
}

int main() {
    //THIS IS HOT!
    const std::array<double,3> center = {5, 0, 0};
    double r = 3.0;
    double obs_radius = 2.9;
    int refine_level = 2;
    int near_field = 3;
    int src_quad_pts = 4;
    int obs_quad_pts = 4;

    auto sphere = sphere_mesh(center, r);
    for (int i = 0; i < refine_level; i++) {
        sphere = refine_mesh(sphere);
    }
    sphere = clean_mesh(sphere);
    auto q_src = tri_gauss(src_quad_pts);
    auto q_obs = tri_double_exp(obs_quad_pts, 0.3);
    KernelFnc K = BEMlaplace_single;
    KernelFnc Kdn = BEMlaplace_double;
    NearEval ne(near_field);

    int n_verts = sphere.vertices.size();
    std::vector<double> u(n_verts);
    for (int i = 0; i < n_verts; i++) {
        u[i] = harmonic_u(sphere.vertices[i]);
    }

    std::vector<double> dudn(n_verts);
    for (int i = 0; i < n_verts; i++) {
        dudn[i] = harmonic_dudn(sphere.vertices[i], center);
    }

    auto rhs = direct_interact(sphere, sphere, q_src, q_obs, 
                               Kdn, u, near_field);
    auto rhs_mass = mass_term(sphere, q_src, u);
    for (unsigned int i = 0; i < rhs.size(); i++){
        rhs[i] = rhs[i] + rhs_mass[i];
    }

    int count = 0;
    auto dudn_solved = solve_system(rhs, 1e-4,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            // TIC
            auto y_temp = direct_interact(sphere, sphere, q_src, q_obs,
                                          K, x, near_field);
            // TOC("Direct interact on " + std::to_string(sphere.faces.size()) + " segments");
            std::copy(y_temp.begin(), y_temp.end(), y.begin());
        });
    std::cout << error(dudn_solved, dudn) << std::endl;
// 
//     double obs_len_scale = get_len_scale(sphere, 0, obs_quad_pts);
//     for(int i = 0; i < 100; i++) {
//         std::array<double,3> obs_pt = random_pt_sphere(center, obs_radius);
// 
//         auto obs_normal = normalize(diff(center, obs_pt));
// 
//         double result = eval_integral_equation(sphere, q_src, Kdn, ne, obs_pt,
//                                                obs_normal, obs_len_scale, u);
//         result += eval_integral_equation(sphere, q_src, K, ne, obs_pt,
//                                          obs_normal, obs_len_scale, dudn);
//         double exact = 1.0 / hypot(obs_pt);
//         double error = std::fabs(exact - result);
//         if (error > 5e-2) {
//             std::cout << "Failed with point: " << Pt{obs_pt} << std::endl;
//         }
//     }
}
