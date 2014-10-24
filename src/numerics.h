#ifndef __NUMERICS_H
#define __NUMERICS_H

#include <vector>
#include <cmath>
#include <functional>

//TODO: Move some of this out of numerics.h
// Map from [-1, 1] to [0, 1].
inline double from_11_to_01(double x) {
    return 0.5 * x + 0.5;
}

// Map from [0, 1] to [-1, 1].
inline double from_01_to_11(double x) {
    return 2 * x - 1;
}

// Map from [a, b] to [-1, 1]
inline double real_to_ref(double x, double a, double b) {
    return 2.0 * ((x - a) / (b - a)) - 1.0;
}

// Map from [-1, 1] to [a, b]
inline double ref_to_real(double x_hat, double a, double b) {
    return a + (b - a) * ((x_hat + 1.0) / 2.0);
}

inline double linear_interp(double x_hat, double y_hat,
                            const std::array<double,3> corner_vals) {
    return (1 - x_hat - y_hat) * corner_vals[0] + 
           x_hat * corner_vals[1] + 
           y_hat * corner_vals[2];
}

inline double hypot2(double x, double y) {return x * x + y * y;}
inline double hypot2(double x, double y, double z) {return x * x + y * y + z * z;}
inline double hypot2(std::array<double,3> v0) {return hypot2(v0[0], v0[1], v0[2]);}
inline double hypot(double x, double y) {return sqrt(hypot2(x, y));}
inline double hypot(double x, double y, double z) {return sqrt(hypot2(x, y, z));}
inline double hypot(const std::array<double,3> a) {return hypot(a[0], a[1], a[2]);}

inline double dist2(double x0, double y0, double z0,
                     double x1, double y1, double z1) {
    return hypot2(x1 - x0, y1 - y0, z1 - z0);
}
template <int dim>
inline double dist2(const std::array<double,dim> v0, 
                    const std::array<double,dim> v1) {
    double result = 0.0;
    for (int d = 0; d < dim; d++) {
        const double delta = (v0[d] - v1[d]);
        result += delta * delta;
    }
    return result;
}

//TODO: Move some of this out to a 3D vector class. 
inline std::array<double,3> cross(const std::array<double,3> x,
                    const std::array<double,3> y) {
    const std::array<double,3> result = {
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    };
    return result;
}

inline std::array<double,3> diff(const std::array<double,3> x,
                                 const std::array<double,3> y) {
    return {x[0] - y[0], x[1] - y[1], x[2] - y[2]};
}
inline std::array<double,3> normalize(const std::array<double,3> x) {
    double mag = hypot(x);
    return {
        x[0] / mag,
        x[1] / mag,
        x[2] / mag
    };
}

inline std::array<double,3> negate(const std::array<double,3> x) {
    return {-x[0],-x[1],-x[2]};
}

inline std::array<double,3> 
tri_unscaled_normal(const std::array<std::array<double,3>,3> corners) {
    return cross(diff(corners[2], corners[0]), diff(corners[2], corners[1]));
}

inline std::array<double,3> 
tri_normal(const std::array<std::array<double,3>,3> corners) {
    auto unscaled = tri_unscaled_normal(corners);
    return normalize(unscaled);
}

inline double tri_area(const std::array<double,3> unscaled_normal) {
    return 0.5 * hypot(unscaled_normal);
}
inline double tri_area(const std::array<std::array<double,3>,3> corners) {
    return tri_area(tri_unscaled_normal(corners));
}
 

typedef std::vector<std::pair<double, double>> QuadratureRule;

std::vector<int> naturals(int min, int max);
std::vector<int> naturals(int max);

std::vector<double> cheb_polys(double x_hat, int n_max);
double s_n(double x_hat, double y_hat, unsigned int n);
double s_n_fast(double x_hat, double y_hat, unsigned int n);
std::vector<double> cheb_pts_first_kind(unsigned int n_pts);
std::pair<double, double> legendre_and_n_minus_1(unsigned int n, double x);

QuadratureRule double_exp(int n, double h);
QuadratureRule gauss(unsigned int n);
QuadratureRule diligenti_mapping(unsigned int n, double x0, int q);
double integrate(QuadratureRule& qr, std::function<double (double)> fnc);

struct QuadratureRule2D {
    QuadratureRule2D(int n): x_hat(n), y_hat(n), weights(n) {}
    std::vector<double> x_hat;
    std::vector<double> y_hat;
    std::vector<double> weights;
};


QuadratureRule2D tensor_product(QuadratureRule xq, QuadratureRule yq);
QuadratureRule2D tensor_gauss(int n_pts);
QuadratureRule2D tensor_double_exp(int n_pts, double h);
QuadratureRule2D tri_gauss(int n_pts);
QuadratureRule2D tri_double_exp(int n_pts, double h);
QuadratureRule2D square_to_tri(QuadratureRule2D square_quad);

double integrate(QuadratureRule2D& qr, std::function<double (double,double)> fnc);
#endif
