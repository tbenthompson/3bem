#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <functional>
#include <array>
#include <vector>
#include "vec.h"
#include "numerics.h"
#include "kernels.h"

template <int dim>
class QuadStrategy;

template <int dim>
struct QuadPt;

template <int dim>
using QuadRule = std::vector<QuadPt<dim>>;

template <int dim>
class Facet;
template <int dim>
class Mesh;

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
    FaceInfo(const Facet<dim>& facet);
    
    //The responsibility is on the user to maintain the lifetime of the facet.
    const Facet<dim>& face;
    const Vec<double,dim> unscaled_n;
    const double area;
    const double jacobian;
    const Vec<double,dim> normal;
};

template <int dim>
struct ObsPt {
    static ObsPt<dim> from_face(const QuadRule<dim-1>& obs_quad,
                           const FaceInfo<dim>& obs_face, int idx);

    const double len_scale;
    const Vec<double,dim> loc;
    const Vec<double,dim> normal;
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
T richardson_step(const std::vector<T>& values);

std::vector<double> integral_equation_vector(const Problem<3>& p,
                                             const QuadStrategy<3>& qs,
                                             const ObsPt<3>& obs);
double eval_integral_equation(const Problem<3>& p, const QuadStrategy<3>& qs,
                              const ObsPt<3>& obs);

//TODO: Use a sparse matrix storage format here.
std::vector<double> interact_matrix(const Problem<3>& p, const QuadStrategy<3>& qs);
std::vector<double> bem_mat_mult(const std::vector<double>& A, 
                                 int n_obs_dofs,
                                 const std::vector<double>& x);
std::vector<double> direct_interact(const Problem<3>& p, const QuadStrategy<3>& qs);
std::vector<double> mass_term(const Problem<3>& p, const QuadStrategy<3>& qs);


double get_len_scale(Mesh<3>& mesh, int which_face, int q);
#endif
