#include <cassert>
#include "bem.h"
#include "mesh.h"
#include "adaptive_quad.h"
#include "quadrature.h"

FaceInfo::FaceInfo(const Facet<3>& facet):
    face(facet),
    unscaled_normal(tri_unscaled_normal(face.vertices)),
    area(tri_area(unscaled_normal)),
    jacobian(area * 2.0),
    normal(unscaled_normal / jacobian)
{}

ObsPt ObsPt::from_face(const QuadRule2d& obs_quad,
                       const FaceInfo& obs_face, int idx) {
    return {
        //TODO: need to divide by q or something like that
        std::sqrt(obs_face.area),
        ref_to_real(obs_quad[idx].x_hat, obs_face.face.vertices),
        obs_face.normal 
    };
}

Vec3<double> eval_quad_pt(const Vec2<double>& x_hat,
                          const Kernel& kernel,
                          const FaceInfo& face,
                          const Vec3<double>& obs_loc,
                          const Vec3<double>& obs_n) {
    const auto src_pt = ref_to_real(x_hat, face.face.vertices);
    const auto d = src_pt - obs_loc;
    const auto r2 = dot(d, d);
    const auto kernel_val = kernel(r2, d, face.normal, obs_n);
    return (kernel_val * face.jacobian) * linear_basis(x_hat);
}

double appx_face_dist2(const Vec3<double>& pt,
                       const std::array<Vec3<double>,3>& vs) {
    double d0 = dist2(pt, vs[0]);
    double d1 = dist2(pt, vs[1]);
    double d2 = dist2(pt, vs[2]);
    return std::min(d0, std::min(d1, d2));
}

/* Perform the richardson extrapolation for the nearfield quadrature. 
 */
template <typename T>
T richardson_step(const std::vector<T>& values) {
    assert(values.size() > 1);
    int n_steps = values.size();
    std::vector<T> last_level = values;
    std::vector<T> this_level;
    for (int m = 1; m < n_steps; m++) {
        this_level.resize(n_steps - m);
        for (int i = 0; i < n_steps - m; i++) {
            const double mult = pow(2, m);
            const double factor = 1.0 / (mult - 1.0);
            T low = last_level[i];
            T high = last_level[i + 1];
            T moreacc = factor * (mult * high - low);
            this_level[i] = moreacc;
        }
        last_level = this_level;
    }
    // std::cout << this_level[0] << std::endl;
    return this_level[0];
}

template double richardson_step(const std::vector<double>&);
template Vec3<double> richardson_step(const std::vector<Vec3<double>>&);

static unsigned long interacts = 0;
static unsigned long far_field_pairs = 0;
static unsigned long near_field_pairs = 0;
static unsigned long adjacent_pairs = 0;  

template <typename T>
T adaptlobstp2(const double a, const double b, 
              const T& fa, const T& fb, const T& is, double outer_x, const Kernel& kernel,
                          const FaceInfo& face,
                          const Vec3<double>& obs_loc,
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

    T fmll = eval_quad_pt(Vec2<double>{outer_x, mll}, kernel, face, obs_loc, obs_n);
    T fml = eval_quad_pt(Vec2<double>{outer_x, ml}, kernel, face, obs_loc, obs_n);
    T fm = eval_quad_pt(Vec2<double>{outer_x, m}, kernel, face, obs_loc, obs_n);
    T fmr = eval_quad_pt(Vec2<double>{outer_x, mr}, kernel, face, obs_loc, obs_n);
    T fmrr = eval_quad_pt(Vec2<double>{outer_x, mrr}, kernel, face, obs_loc, obs_n);
    interacts += 5;

    T i2 = (h / 6.) * (fa + fb + 5.0 * (fml + fmr));    
    T i1 = (h / 1470.) * (
            77.0 * (fa + fb) + 
            432.0 * (fmll + fmrr) +
            625.0 * (fml + fmr) +
            672.0 * fm);

    if (all(is + (i1 - i2) == is)) {
        return i1;
    } else if (mll <= a or b <= mrr) {
        std::cout << "YIKES!" << std::endl;
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

//TODO: Refactor the shit out of this! Super ugly.
template <typename T>
T adaptive_integrate2(double a, double b, 
                      double p_tol, double outer_x, const Kernel& kernel,
                          const FaceInfo& face,
                          const Vec3<double>& obs_loc,
                          const Vec3<double>& obs_n)
{
    double m = (a + b) / 2.; 
    double h = (b - a) / 2.;

    const T y[13] = {
        eval_quad_pt(Vec2<double>{outer_x, a}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m - lobatto_x1 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m - lobatto_alpha * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m - lobatto_x2 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m - lobatto_beta * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m - lobatto_x3 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m + lobatto_x3 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m + lobatto_beta * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m + lobatto_x2 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m + lobatto_alpha * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, m + lobatto_x1 * h}, kernel, face, obs_loc, obs_n),
        eval_quad_pt(Vec2<double>{outer_x, b}, kernel, face, obs_loc, obs_n)
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
    interacts += 13;

    return adaptlobstp2(a, b, fa, fb, err_is, outer_x, kernel, face, obs_loc, obs_n);
}

Vec3<double> near_field(const Problem& p, const QuadStrategy& qs,
                        const ObsPt& obs, const FaceInfo& src_face,
                        const double dist2) {
    std::vector<Vec3<double>> near_steps(qs.n_singular_steps, {0,0,0});
    for (int nf = 0; nf < qs.n_singular_steps; nf++) {
        double nfdn = 5 * obs.len_scale * qs.singular_steps[nf];
        auto nf_obs_pt = obs.loc + nfdn * obs.normal;
         
        if (dist2 > 0.5 * src_face.area) {
            for (std::size_t i = 0; i < qs.src_near_quad.size(); i++) {
                near_steps[nf] += qs.src_near_quad[i].w *
                    eval_quad_pt(qs.src_near_quad[i].x_hat, p.K, src_face,
                                 nf_obs_pt, obs.normal);
            }
            interacts += qs.src_near_quad.size();
            near_field_pairs++;
        } else {
            Vec3<double> ns = adaptive_integrate<Vec3<double>>(
                [&] (double x_hat) {
                    if (x_hat == 1.0) {
                        return Vec3<double>{0.0,0.0,0.0};
                    }
                    return adaptive_integrate2<Vec3<double>>(0.0, 1 - x_hat, qs.singular_tol, x_hat, p.K, src_face, nf_obs_pt, obs.normal);
                }, 0.0, 1.0, qs.singular_tol);
            near_steps[nf] += ns;
            adjacent_pairs++;
        }
    }
    return richardson_step(near_steps);
}

std::vector<double> integral_equation_vector(const Problem& p,
                                             const QuadStrategy& qs,
                                             const ObsPt& obs) {
    std::vector<double> result(3 * p.src_mesh.facets.size(), 0.0);
    for (std::size_t i = 0; i < p.src_mesh.facets.size(); i++) {
        FaceInfo src_face(p.src_mesh.facets[i]);
        const double dist2 = appx_face_dist2(obs.loc, src_face.face.vertices);

        Vec3<double> integrals;
        if (dist2 < pow(qs.far_threshold, 2) * src_face.area) {
            integrals = near_field(p, qs, obs, src_face, dist2);
        } else {
            // farfield
            integrals = {0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < qs.src_far_quad.size(); i++) {
                integrals += qs.src_far_quad[i].w *
                    eval_quad_pt(qs.src_far_quad[i].x_hat, p.K, src_face,
                                 obs.loc, obs.normal);
            }
            interacts += qs.src_far_quad.size();
            far_field_pairs++;
        }
        for (int b = 0; b < 3; b++) {
            result[3 * i + b] += integrals[b];
        }
    }
    return result;
}

/* Evaluate the integral equation for a specific observation point.
 */
double eval_integral_equation(const Problem& p, const QuadStrategy& qs,
                              const ObsPt& obs) {
    double result = 0.0;
    auto row = integral_equation_vector(p, qs, obs);
    for (std::size_t i = 0; i < row.size(); i++) {
        result += row[i] * p.src_strength[i];
    }
    return result;
}

//TODO: Use a sparse matrix storage format here.
std::vector<double> interact_matrix(const Problem& p,
                                    const QuadStrategy& qs) {
    int n_obs_dofs = 3 * p.obs_mesh.facets.size();
    int n_src_dofs = 3 * p.src_mesh.facets.size();
    std::cout << "HI" << std::endl;
    std::cout << n_obs_dofs << " " << n_src_dofs << std::endl;
    std::vector<double> matrix(n_obs_dofs * n_src_dofs, 0.0);
    std::cout << "HI2" << std::endl;
#pragma omp parallel for
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh.facets[obs_idx]);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto pt = ObsPt::from_face(qs.obs_quad, obs_face, obs_q);

            const auto row = integral_equation_vector(p, qs, pt);

            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp(qs.obs_quad[obs_q].x_hat,
                                                      unit<double>(v)); 
                int b = 3 * obs_idx + v;
                for (int i = 0; i < n_src_dofs; i++) {
                    matrix[b * n_obs_dofs + i] += obs_basis_eval * row[i] * 
                                    qs.obs_quad[obs_q].w * obs_face.jacobian;
                }
            }
        }
    }
    std::cout << interacts << " " << far_field_pairs << " " << near_field_pairs << " " << adjacent_pairs << std::endl;
    return matrix;
}

std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x) {
    std::vector<double> res(n_obs_dofs, 0.0);
#pragma omp parallel for
    for (int i = 0; i < n_obs_dofs; i++) {
        for (std::size_t j = 0; j < x.size(); j++) {
            res[i] += A[i * n_obs_dofs + j] * x[j]; 
        }
    }
    return res;
}

std::vector<double> direct_interact(const Problem& p,
                                    const QuadStrategy& qs) {
    return bem_mat_mult(interact_matrix(p, qs), 
                        p.obs_mesh.facets.size() * 3, p.src_strength);
}

std::vector<double> mass_term(const Problem& p,
                              const QuadStrategy& qs) {
    int n_obs_dofs = 3 * p.obs_mesh.facets.size();
    std::vector<double> integrals(n_obs_dofs, 0.0);
    for (std::size_t obs_idx = 0; obs_idx < p.obs_mesh.facets.size(); obs_idx++) {
        FaceInfo obs_face(p.obs_mesh.facets[obs_idx]);
        for (std::size_t obs_q = 0; obs_q < qs.obs_quad.size(); obs_q++) {
            auto qpt = qs.obs_quad[obs_q];
            int dof = 3 * obs_idx;
            Vec3<double> face_vals = {
                p.src_strength[dof], p.src_strength[dof + 1], p.src_strength[dof + 2]
            };
            double interp_val = linear_interp(qpt.x_hat, face_vals);

            for(int v = 0; v < 3; v++) {
                double obs_basis_eval = linear_interp(qpt.x_hat, unit<double>(v)); 
                integrals[dof + v] += obs_face.jacobian * obs_basis_eval * 
                                      interp_val * qpt.w;
            }
        }
    }
    return integrals;
}

double get_len_scale(Mesh<3>& mesh, int which_face, int q) {
    return std::sqrt(tri_area(mesh.facets[which_face].vertices)) / q;
}
