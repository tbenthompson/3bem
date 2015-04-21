#include "fmm.h"

namespace tbem {
LUDecomposition LU_decompose(const std::vector<double>& matrix) 
{
    int n = std::sqrt(matrix.size()); 
    std::vector<double> A = matrix;
    std::vector<int> pivots(n);
    int info;
    dgetrf_(&n, &n, A.data(), &n, pivots.data(), &info);
    return {A, pivots};
}

std::vector<double> LU_solve(LUDecomposition& lu,
    const std::vector<double>& b) 
{
    std::vector<double> x = b;
    int n = b.size(); 
    int n_rhs = 1;
    int info;
    // type = 'T' means that we solve A^T x = b rather than Ax = b. This is good
    // because blas operates on column major data while this code sticks to 
    // row major data.
    char type = 'T';
    dgetrs_(&type, &n, &n_rhs, lu.LU.data(), &n, lu.pivots.data(), x.data(), &n, &info);
    return x;
}

template <>
TranslationSurface<2> make_surrounding_surface<2>(size_t expansion_order) 
{
    std::vector<Vec<double,2>> pts;
    std::vector<Vec<double,2>> normals;

    auto theta = linspace(0, 2 * M_PI, expansion_order + 1);
    for (size_t i = 0; i < expansion_order; i++) {
        Vec<double,2> pt{std::cos(theta[i]), std::sin(theta[i])}; 
        Vec<double,2> normal = -pt;
        pts.push_back(pt);
    }

    return {pts, pts};
}

TranslationSurface<3> surrounding_surface_sphere(size_t expansion_order)
{
    std::vector<Vec<double,3>> pts;
    double a = 4 * M_PI / expansion_order;
    double d = std::sqrt(a);
    auto M_theta = static_cast<size_t>(std::round(M_PI / d));
    double d_theta = M_PI / M_theta;
    double d_phi = a / d_theta;
    for (size_t m = 0; m < M_theta; m++) {
        double theta = M_PI * (m + 0.5) / M_theta;
        auto M_phi = static_cast<size_t>(
            std::round(2 * M_PI * std::sin(theta) / d_phi)
        );
        for (size_t n = 0; n < M_phi; n++) {
            double phi = 2 * M_PI * n / M_phi;
            double x = std::sin(theta) * std::cos(phi);
            double y = std::sin(theta) * std::sin(phi);
            double z = std::cos(theta);
            pts.push_back({x, y, z});
        }
    }

    return {pts, pts};
}

template <>
TranslationSurface<3> make_surrounding_surface<3>(size_t expansion_order) 
{
    return surrounding_surface_sphere(expansion_order);
}

template <size_t dim>
std::vector<Vec<double,dim>> 
TranslationSurface<dim>::move(const Box<dim>& box, double r_ref) const
{
    auto cell_r = hypot(box.half_width);
    auto out_r = r_ref * cell_r;
    auto out_center = box.center;

    std::vector<Vec<double,dim>> out_pts(pts.size());
    for (size_t i = 0; i < pts.size(); i++) {
        out_pts[i] = out_r * pts[i] + out_center;
    }

    return out_pts;
}

template struct TranslationSurface<2>;
template struct TranslationSurface<3>;

} // END namespace tbem
