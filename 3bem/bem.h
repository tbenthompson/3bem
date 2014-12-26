#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <cassert>
#include <functional>
#include <array>
#include <vector>
#include "vec.h"
#include "numerics.h"
#include "kernels.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "quadrature.h"


namespace tbem {

template <int dim>
struct Problem {
    const Mesh<dim>& src_mesh;
    const Mesh<dim>& obs_mesh;
    const Kernel<dim>& K;
    const std::vector<double>& src_strength;
};

template <int dim>
class FaceInfo {
public:
    //The responsibility is on the user to maintain the lifetime of the facet.
    const Facet<dim>& face;

    const double area_scale;
    const double length_scale;
    const double jacobian;
    const Vec<double,dim> normal;

    static const double inv_ref_facet_area;
    static FaceInfo<dim> build(const Facet<dim>& facet);
};

template <> const double FaceInfo<3>::inv_ref_facet_area = 2.0;
template <> const double FaceInfo<2>::inv_ref_facet_area = 0.5;

template <>
FaceInfo<3> FaceInfo<3>::build(const Facet<3>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area = tri_area(unscaled_n);
    auto length_scale = std::sqrt(area);
    auto jacobian = area * inv_ref_facet_area;
    auto normal = unscaled_n / jacobian;
    return FaceInfo<3>{facet, area, length_scale, jacobian, normal};
}

template <>
FaceInfo<2> FaceInfo<2>::build(const Facet<2>& facet){
    auto unscaled_n = unscaled_normal(facet.vertices);
    auto area_scale = hypot2(unscaled_n);
    auto length = std::sqrt(area_scale);
    auto jacobian = length * inv_ref_facet_area;
    auto normal = unscaled_n / length;
    return FaceInfo<2>{facet, area_scale, length, jacobian, normal};
}

template <int dim>
struct ObsPt {
    static ObsPt<dim> from_face(const QuadRule<dim-1>& obs_quad,
                                const FaceInfo<dim>& obs_face,
                                int idx) {
        return {
            obs_face.length_scale,
            ref_to_real(obs_quad[idx].x_hat, obs_face.face.vertices),
            obs_face.normal,
            obs_face.normal 
        };
    }

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
    const Vec<double,dim> richardson_dir;
};

template <int dim>
Vec<double,dim> eval_quad_pt(const Vec<double,dim-1>& x_hat,
                          const Kernel<dim>& kernel,
                          const FaceInfo<dim>& face,
                          const Vec<double,dim>& obs_loc,
                          const Vec<double,dim>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.face.vertices);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return (kernel_val * face.jacobian) * linear_basis(x_hat);
}

template <typename T>
T adaptlobstp2(const double a, const double b, 
              const T& fa, const T& fb, const T& is, double outer_x, 
              const Kernel<3>& kernel, const FaceInfo<3>& face, const Vec3<double>& obs_loc,
                const Vec3<double>& obs_n)
{
    // std::cout << a << " " << b << std::endl;
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    double ah = lobatto_alpha * h;
    double bh = lobatto_beta * h;
    double mll = m - ah;
    double ml = m - bh;
    double mr = m + bh;
    double mrr = m + ah;

    T fmll = eval_quad_pt<3>(Vec2<double>{outer_x, mll}, kernel, face, obs_loc, obs_n);
    T fml = eval_quad_pt<3>(Vec2<double>{outer_x, ml}, kernel, face, obs_loc, obs_n);
    T fm = eval_quad_pt<3>(Vec2<double>{outer_x, m}, kernel, face, obs_loc, obs_n);
    T fmr = eval_quad_pt<3>(Vec2<double>{outer_x, mr}, kernel, face, obs_loc, obs_n);
    T fmrr = eval_quad_pt<3>(Vec2<double>{outer_x, mrr}, kernel, face, obs_loc, obs_n);

    T i2 = (h / 6.) * (fa + fb + 5.0 * (fml + fmr));    
    T i1 = (h / 1470.) * (
            77.0 * (fa + fb) + 
            432.0 * (fmll + fmrr) +
            625.0 * (fml + fmr) +
            672.0 * fm);

    if (all(is + (i1 - i2) == is)) {
        return i1;
    } else if (mll <= a or b <= mrr) {
        std::cout << "YIKES FROM ADAPTIVE!" << std::endl;
        return i1;
    } else {
        return adaptlobstp2(a, mll, fa, fmll, is, outer_x, kernel, face, obs_loc, obs_n)
             + adaptlobstp2(mll, ml, fmll, fml, is, outer_x, kernel, face, obs_loc, obs_n)
             + adaptlobstp2(ml, m, fml, fm, is, outer_x, kernel, face, obs_loc, obs_n)
             + adaptlobstp2(m, mr, fm, fmr, is, outer_x, kernel, face, obs_loc, obs_n)
             + adaptlobstp2(mr, mrr, fmr, fmrr, is, outer_x, kernel, face, obs_loc, obs_n)
             + adaptlobstp2(mrr, b, fmrr, fb, is, outer_x, kernel, face, obs_loc, obs_n);
    }
}

template <typename T>
T adaptive_integrate2(double a, double b, 
                      double p_tol, double outer_x, const Kernel<3>& kernel,
                          const FaceInfo<3>& face,
                          const Vec3<double>& obs_loc,
                          const Vec3<double>& obs_n)
{
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    const T y[13] = {
        eval_quad_pt<3>(Vec2<double>{outer_x, a}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m - lobatto_x1 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m - lobatto_alpha * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m - lobatto_x2 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m - lobatto_beta * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m - lobatto_x3 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m + lobatto_x3 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m + lobatto_beta * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m + lobatto_x2 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m + lobatto_alpha * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, m + lobatto_x1 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt<3>(Vec2<double>{outer_x, b}, kernel, face, obs_loc, obs_n)
    };
    
    const T fa = y[0];
    const T fb = y[12];
    
    const T i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));
    const T i1 = (h / 1470.0) * (
            77.0 * (y[0] + y[12]) +
            432.0 * (y[2] + y[10]) +
            625.0 * (y[4] + y[8]) +
            672.0 * y[6]);
    const T is = h * (
        0.0158271919734802 * (y[0] + y[12]) + 
        0.0942738402188500 * (y[1] + y[11]) + 
        0.155071987336585  * (y[2] + y[10]) +
        0.188821573960182  * (y[3] + y[9]) + 
        0.199773405226859  * (y[4] + y[8]) +
        0.224926465333340  * (y[5] + y[7]) + 
        0.242611071901408  * y[6]);    
   
    const T erri1 = fabs(i1 - is);
    const T erri2 = fabs(i2 - is);
    
    const T err_is = get_error_is(p_tol, erri1, erri2, is, a, b);

    return adaptlobstp2(a, b, fa, fb, err_is, outer_x, kernel, face, obs_loc, obs_n);
}


template <typename T>
T richardson_step(const std::vector<T>& values);

std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x);


template <int dim> 
Vec<double,dim> adaptive_nearfield(const Problem<dim>& p,
                                    const QuadStrategy<dim>& qs,
                                    const ObsPt<dim>& obs,
                                    const FaceInfo<dim>& src_face,
                                    const Vec<double,dim>& nf_obs_pt);

template <>
Vec<double,3> adaptive_nearfield<3>(const Problem<3>& p,
                                    const QuadStrategy<3>& qs,
                                    const ObsPt<3>& obs,
                                    const FaceInfo<3>& src_face,
                                    const Vec<double,3>& nf_obs_pt) {
    return adaptive_integrate<Vec<double,3>>(
        [&] (double x_hat) {
            if (x_hat == 1.0) {
                return zeros<Vec<double,3>>();
            }
            return adaptive_integrate2<Vec<double,3>>(
                        0.0, 1 - x_hat, qs.near_tol, x_hat, p.K,
                        src_face, nf_obs_pt, obs.normal);
        }, 0.0, 1.0, qs.near_tol);
}

template <>
Vec<double,2> adaptive_nearfield<2>(const Problem<2>& p,
                                    const QuadStrategy<2>& qs,
                                    const ObsPt<2>& obs,
                                    const FaceInfo<2>& src_face,
                                    const Vec<double,2>& nf_obs_pt) {
    return adaptive_integrate<Vec<double,2>>(
        [&] (double x_hat) {
            return eval_quad_pt<2>(Vec<double,1>{x_hat}, p.K, src_face,
                                   nf_obs_pt, obs.normal);
        }, -1.0, 1.0, qs.near_tol);
}

template <int dim>
Vec<double,dim> near_field(const Problem<dim>& p, const QuadStrategy<dim>& qs,
                        const ObsPt<dim>& obs, const FaceInfo<dim>& src_face,
                        const double dist2) {
    std::vector<Vec<double,dim>> near_steps(qs.n_singular_steps, 
                                            zeros<Vec<double, dim>>());
    const double singular_threshold = 3.0;
    const double safe_dist_ratio = 5.0;
    if (dist2 < singular_threshold * src_face.area_scale) { 
        for (int nf = 0; nf < qs.n_singular_steps; nf++) {
            double nfdn = safe_dist_ratio * obs.len_scale * qs.singular_steps[nf];
            auto nf_obs_pt = obs.loc + nfdn * obs.richardson_dir;
            auto ns = adaptive_nearfield<dim>(p, qs, obs, src_face, nf_obs_pt);
            near_steps[nf] += ns;
        }
        return richardson_step(near_steps);
    } else {
        return adaptive_nearfield<dim>(p, qs, obs, src_face, obs.loc);
    }
}

/* It may be that the exact distance to an element is not required,
 * but a low precision approximate distance is useful.
 * This function simply approximates the distance by the distance from
 * the given point to each vertex of the element.
 */
template <int dim>
double appx_face_dist2(const Vec<double,dim>& pt,
                       const std::array<Vec<double,dim>,dim>& vs) {
    double res = dist2(pt, vs[0]);
    for (int d = 1; d < dim; d++) {
        res = std::min(res, dist2(pt, vs[d])); 
    }
    return res;
}

/* Given an unknown function defined in terms of a polynomial basis:
 * u(x) = \sum_i U_i \phi_i(x)
 * Determine the coefficients for the expansion of an integral term in
 * terms of U_i.
 * In other words, this function calculates the vector Q_i where
 * Q_i = \int_{S_{src}} K(x,y) \phi_i(y) dy
 */
template <int dim>
std::vector<double> integral_equation_vector(const Problem<dim>& p,
                                             const QuadStrategy<dim>& qs,
                                             const ObsPt<dim>& obs) {
    int n_out_dofs = dim * p.src_mesh.facets.size();
    std::vector<double> result(n_out_dofs);
    for (std::size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        auto src_face = FaceInfo<dim>::build(p.src_mesh.facets[i]);
        const double dist2 = appx_face_dist2<dim>(obs.loc, src_face.face.vertices);

        Vec<double,dim> integrals;
        if (dist2 < pow(qs.far_threshold, 2) * src_face.area_scale) {
            integrals = near_field(p, qs, obs, src_face, dist2);
        } else {
            // farfield
            integrals = zeros<Vec<double,dim>>();
            for (std::size_t i = 0; i < qs.src_far_quad.size(); i++) {
                integrals += qs.src_far_quad[i].w *
                    eval_quad_pt<dim>(qs.src_far_quad[i].x_hat, p.K, src_face,
                                 obs.loc, obs.normal);
            }
        }
        for (int b = 0; b < dim; b++) {
            result[dim * i + b] = integrals[b];
        }
    }
    return result;
}

/* Evaluate the integral equation for a specific observation point:
 * \int_{S_{src}} K(x,y) u(y) dy
 * y is given by the ObsPt<dim> obs.
 * The caller provides a quadrature strategy that specifies the order and 
 * tolerance for evaluating the integral equation.
 */
template <int dim>
double eval_integral_equation(const Problem<dim>& p, const QuadStrategy<dim>& qs,
                              const ObsPt<dim>& obs) {
    double result = 0.0;
    auto row = integral_equation_vector(p, qs, obs);
    for (std::size_t i = 0; i < row.size(); i++) {
        result += row[i] * p.src_strength[i];
    }
    return result;
}

/* Given a kernel function and two meshes this function calculates the
 * Galerkin boundary element matrix representing the operator 
 * \int_{S_{obs}} \phi_i(x) \int_{S_{src}} K(x,y) \phi_j(y) dy dx
 * where S_{obs} is the observation mesh, S_{src} is the source mesh,
 * K(x,y) is the kernel function and \phi_i(x) is a basis function.
 */
template <int dim>
std::vector<double> interact_matrix(const Problem<dim>& p,
                                    const QuadStrategy<dim>& qs) {
    std::size_t n_obs_dofs = dim * p.obs_mesh.facets.size();
    std::size_t n_src_dofs = dim * p.src_mesh.facets.size();
    std::vector<double> matrix(n_obs_dofs * n_src_dofs, 0.0);
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FaceInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt<dim>::from_face(qs.obs_quad, obs_face, obs_q);

            const auto row = integral_equation_vector(p, qs, pt);

            const auto basis = linear_basis(qs.obs_quad[obs_q].x_hat);

            for (int v = 0; v < dim; v++) {
                int b = dim * obs_idx + v;
                for (std::size_t i = 0; i < n_src_dofs; i++) {
                    double contribution = basis[v] * row[i] * qs.obs_quad[obs_q].w *
                                      obs_face.jacobian;
                    matrix[b * n_src_dofs + i] += contribution;
                }
            }
        }
    }
    return matrix;
}

template <int dim>
std::vector<double> direct_interact(const Problem<dim>& p,
                                    const QuadStrategy<dim>& qs) {
    auto matrix = interact_matrix(p, qs);
    assert(p.obs_mesh.facets.size() * dim * p.src_strength.size() == matrix.size());
    return bem_mat_mult(matrix, 
                        p.obs_mesh.facets.size() * dim,
                        p.src_strength);
}

/* In many integral equations, one of the functions of interest appears
 * outside of an integration, possibly as the result of integrating against
 * a delta function. 
 * In a galerkin boundary element formulation, these free terms look like:
 * \int_S \phi_i(x) u(x) dx.
 * This function calculates such integrals using gaussian quadrature.
 */
template <int dim>
std::vector<double> mass_term(const Problem<dim>& p,
                              const QuadStrategy<dim>& qs) {
    int n_obs_dofs = dim * p.obs_mesh.facets.size();
    std::vector<double> integrals(n_obs_dofs, 0.0);
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        auto obs_face = FaceInfo<dim>::build(p.obs_mesh.facets[obs_idx]);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];
            int dof = dim * obs_idx;
            Vec<double,dim> face_vals;
            for (int d = 0; d < dim; d++) {
                face_vals[d] = p.src_strength[dof + d];
            }
            double interp_val = linear_interp<dim>(qpt.x_hat, face_vals);

            auto basis = linear_basis(qpt.x_hat);

            for(int v = 0; v < dim; v++) {
                integrals[dof + v] += obs_face.jacobian * basis[v] * 
                                      interp_val * qpt.w;
            }
        }
    }
    return integrals;
}


template <int dim>
double get_len_scale(Mesh<dim>& mesh, int which_face, int q);

template <>
double get_len_scale<3>(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face].vertices)) / q;
}

template <>
double get_len_scale<2>(Mesh<2>& mesh, int which_face, int q) {
    return dist(mesh.facets[which_face].vertices[1],
                mesh.facets[which_face].vertices[0]) / q;
}

} // END NAMESPACE tbem
#endif
