#ifndef __AAAAAAAAA_BEM_H
#define __AAAAAAAAA_BEM_H

#include <functional>
#include <array>
#include <vector>
#include "vec.h"
#include "numerics.h"
#include "quadrature.h"

class Mesh;

typedef std::function<double(const double&, const Vec3<double>&,
                             const Vec3<double>&, const Vec3<double>&)> Kernel;

struct Problem {
    const Mesh& src_mesh;
    const Mesh& obs_mesh;
    const Kernel& K;
    const std::vector<double>& src_strength;
};

class FaceInfo {
public:
    FaceInfo(const Mesh& mesh, int face_index);
    
    const int face_index;
    const std::array<int,3>& face;
    const std::array<Vec3<double>,3> corners;
    const Vec3<double> unscaled_normal;
    const double area;
    const double jacobian;
    const Vec3<double> normal;
};

struct ObsPt {
    static ObsPt from_face(const QuadRule2d& obs_quad,
                           const FaceInfo& obs_face, int idx);

    const double len_scale;
    const Vec3<double> loc;
    const Vec3<double> normal;
};

struct SrcPointInfo {
    SrcPointInfo(const QuadRule2d& quad_rule,
              const Kernel& kernel,
              const FaceInfo& face,
              const Vec3<double>& obs_loc,
              const Vec3<double>& obs_n,
              int q_index);
    
    //TODO: This seems slightly vacuous...
    Vec3<double> basis;
    double weighted_kernel;
};

template <typename T>
T richardson_step(const std::vector<T>& values);

Vec3<double> basis_integrals(const QuadRule2d& quad_rule,
                        const Kernel& kernel,
                        const FaceInfo& face,
                        const Vec3<double>& obs_loc,
                        const Vec3<double>& obs_n);

double integral(const QuadRule2d& quad_rule,
           const Kernel& kernel,
           const FaceInfo& face,
           const Vec3<double>& src_vals,
           const Vec3<double>& obs_loc,
           const Vec3<double>& obs_n);

std::vector<double> integral_equation_vector(const Problem& p,
                                             const QuadStrategy& qs,
                                             const ObsPt& obs);
double eval_integral_equation(const Problem& p, const QuadStrategy& qs,
                              const ObsPt& obs);

//TODO: Use a sparse matrix storage format here.
typedef std::vector<std::vector<double>> Mat;
Mat interact_matrix(const Problem& p, const QuadStrategy& qs);
std::vector<double> bem_mat_mult(const Mat& A, const std::vector<double>& x);
std::vector<double> direct_interact(const Problem& p, const QuadStrategy& qs);
std::vector<double> mass_term(const Problem& p, const QuadStrategy& qs);


double get_len_scale(Mesh& mesh, int which_face, int q);
#endif
