#include "sinh_quad.h"
#include "numerics.h"
#include <cmath>

namespace tbem {

QuadRule<2> sinh_sigmoidal_transform(const QuadRule<1>& gauss_theta,
    const QuadRule<1>& gauss_r, double b, bool iterated_sinh) 
{
    auto sinh1d = sinh_transform(gauss_r, -1.0, b, iterated_sinh); 
    QuadRule<2> out;
    for (size_t i = 0; i < gauss_theta.size(); i++) {
        double sigma = (gauss_theta[i].x_hat[0] + 1.0) / 2.0;
        double sig_transform = std::pow(sigma, 2) /
            (std::pow(sigma, 2) + std::pow(1 - sigma, 2));
        double theta = sig_transform * (M_PI / 2.0);
        double theta_jacobian = (M_PI / 2.0) * sigma * (1 - sigma) /
            std::pow(std::pow(sigma, 2) + std::pow(1 - sigma, 2), 2);

        double R_theta = std::sin(M_PI / 4.0) / (std::sin(theta + M_PI / 4.0));

        for (size_t j = 0; j < sinh1d.size(); j++) {
            double s = sinh1d[j].x_hat[0];
            double r = R_theta * (s + 1.0) / 2.0; 
            double r_jacobian = r * (R_theta / 2.0);
            double w = theta_jacobian * r_jacobian * sinh1d[j].w * gauss_theta[i].w;
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);
            out.push_back({{x, y}, w});
        }
    }

    return out;
};

QuadRule<2> transform_to_tri(const QuadRule<2>& q, Vec<Vec<double,3>,3> tri) {
    QuadRule<2> out;
    double jacobian = 2.0 * tri_area(tri);
    for (const auto& unit_facet_pt: q) {
        auto pt = ref_to_real(unit_facet_pt.x_hat, tri); 
        assert(pt[2] == 0.0);
        out.push_back({{pt[0], pt[1]}, unit_facet_pt.w * jacobian});
    }
    return out;
}

QuadRule<2> sinh_sigmoidal_transform(const QuadRule<1>& gauss_theta,
    const QuadRule<1>& gauss_r, double x0,
    double y0, double b, bool iterated_sinh) 
{
    Vec<double,3> pt0{0, 0, 0};
    Vec<double,3> pt1{1, 0, 0};
    Vec<double,3> pt2{0, 1, 0};
    Vec<double,3> singular_pt{x0, y0, 0};

    auto unit_facet_rule =
        sinh_sigmoidal_transform(gauss_theta, gauss_r, b, iterated_sinh);
    QuadRule<2> out;
    
    // upper left tri
    auto ul_pts = transform_to_tri(unit_facet_rule, {singular_pt, pt0, pt2});
    for (const auto& p: ul_pts) {out.push_back(p);}

    // upper right tri
    auto ur_pts = transform_to_tri(unit_facet_rule, {singular_pt, pt2, pt1});
    for (const auto& p: ur_pts) {out.push_back(p);}

    // lower
    auto l_pts = transform_to_tri(unit_facet_rule, {singular_pt, pt0, pt1});
    for (const auto& p: l_pts) {out.push_back(p);}

    return out;
}

} //end namespace tbem
