#include "laplace.h"

using namespace tbem;

double error_inf(const std::vector<double>& a, const std::vector<double>& b) {
    double error = 0.0;
    for (unsigned int i = 0; i < a.size(); i++) {
        error = std::max(std::fabs(a[i] - b[i]), error);
    }
    return error;
}

double random_val() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

Vec3<double> spherify_pt(Vec3<double> pt, Vec3<double> c, double r) {
    double rf = pt[0];
    double tf = pt[1] * 2 * M_PI;
    double pf = pt[2] * M_PI;
    Vec3<double> ret = {
        c[0] + r * rf * std::cos(tf) * std::sin(pf),
        c[1] + r * rf * std::sin(tf) * std::sin(pf),
        c[2] + r * rf * std::cos(pf)
    };
    return ret;
}

Vec3<double> random_pt_sphere(Vec3<double> c, double r) {
    auto pt = random_pt3d();
    return spherify_pt(pt, c, r);
}

Vec2<double> random_pt_sphere(Vec2<double> c, double r) {
    auto sphere_pt = random_pt_sphere({c[0], c[1], 0.0}, r);
    return {sphere_pt[0], sphere_pt[1]};
}

