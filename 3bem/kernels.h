#ifndef __YUIQWEOIUQWE_KERNELS_H
#define __YUIQWEOIUQWE_KERNELS_H

#include <array>
#include "vec.h"

namespace tbem { 

template <int dim>
struct One {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;
    double operator()(const double& r2, const Vec<double,dim>& delta,
                 const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs) const {
        return 1.0;
    }
};


/* Laplace/Poisson equation kernels. */
template <int dim>
struct LaplaceSingle;
template <int dim>
struct LaplaceDouble;
template <int dim>
struct LaplaceHypersingular;

template <>
struct LaplaceSingle<3> {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    double operator()(const double& r2, const Vec<double,3>& delta,
                      const Vec<double,3>& nsrc, const Vec<double,3>& nobs) const {
        return 1.0 / (4.0 * M_PI * std::sqrt(r2));
    }
};

template <>
struct LaplaceDouble<3> {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    double operator()(const double& r2, const Vec<double,3>& delta,
                      const Vec<double,3>& nsrc, const Vec<double,3>& nobs) const {
        return dot_product(nsrc, delta) / (4 * M_PI * r2 * std::sqrt(r2));
    }
};

template <>
struct LaplaceSingle<2> {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    double operator()(const double& r2, const Vec<double,2>& delta,
                      const Vec<double,2>& nsrc, const Vec<double,2>& nobs) const {
        return std::log(std::sqrt(r2)) / (2 * M_PI);
    }
};

template <>
struct LaplaceDouble<2> {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    double operator()(const double& r2, const Vec<double,2>& delta,
                      const Vec<double,2>& nsrc, const Vec<double,2>& nobs) const {
        return dot_product(nsrc, delta) / (2 * M_PI * r2);
    }
};

template <>
struct LaplaceHypersingular<2> {
    typedef double OutType;
    typedef double InType;
    typedef double OperatorType;

    double operator()(const double& r2, const Vec<double,2>& delta,
                      const Vec<double,2>& nsrc, const Vec<double,2>& nobs) const {
        return ((-dot_product(nobs, nsrc) / r2) + 
                ((2 * dot_product(nsrc, delta) * dot_product(nobs, delta)) / (r2 * r2)))
            / (2 * M_PI);
    }
};

const double kronecker[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

template <int dim>
struct ElasticDisplacement;
template <int dim>
struct ElasticTraction;
template <int dim>
struct ElasticAdjointTraction;
template <int dim>
struct ElasticHypersingular;

template <>
struct ElasticDisplacement<2> {
    const double disp_C1;
    const double disp_C2;
    typedef Vec2<double> OutType;
    typedef Vec2<double> InType;
    typedef Vec2<Vec2<double>> OperatorType;

    ElasticDisplacement(double shear_modulus, double poisson_ratio):
        disp_C1(1.0 / (8 * M_PI * shear_modulus * (1 - poisson_ratio))),
        disp_C2(3 - 4 * poisson_ratio)
    {}

    OperatorType operator()(double r2, const Vec2<double>& delta, 
                       const Vec2<double>& nsrc, const Vec2<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                out[k][j] = disp_C1 * (-disp_C2 * kronecker[k][j] * std::log(r)  + 
                                       delta[k] * delta[j] / r2);
            }
        }
        return out;
    }
};

template <>
struct ElasticTraction<2> {
    const double trac_C1;
    const double trac_C2;
    typedef Vec2<double> OutType;
    typedef Vec2<double> InType;
    typedef Vec2<Vec2<double>> OperatorType;

    ElasticTraction(double shear_modulus, double poisson_ratio):
        trac_C1(1.0 / (4 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
    {}
    
    OperatorType operator()(double r2, const Vec2<double>& delta, 
                       const Vec2<double>& nsrc, const Vec2<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        const auto drdn = dot_product(delta, nsrc) / r;
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                const auto term1 = (trac_C2 * kronecker[k][j] +
                                    2 * delta[k] * delta[j] / r2);
                const auto term2 = trac_C2 * (nsrc[j] * delta[k] - 
                                              nsrc[k] * delta[j]) / r;
                out[k][j] = -(trac_C1 / r) * (term1 * drdn - term2);
            }
        }
        return out;
    }
};

template <>
struct ElasticAdjointTraction<2> {
    const double trac_C1;
    const double trac_C2;
    typedef Vec2<double> OutType;
    typedef Vec2<double> InType;
    typedef Vec2<Vec2<double>> OperatorType;

    ElasticAdjointTraction(double shear_modulus, double poisson_ratio):
        trac_C1(1.0 / (4 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
    {}
    
    OperatorType operator()(double r2, const Vec2<double>& delta, 
                       const Vec2<double>& nsrc, const Vec2<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        const auto drdm = dot_product(delta, nobs) / r;
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                const auto term1 = (trac_C2 * kronecker[k][j] +
                                    2 * delta[k] * delta[j] / r2);
                const auto term2 = trac_C2 * (nobs[j] * delta[k] -
                                              nobs[k] * delta[j]) / r;
                out[k][j] = (trac_C1 / r) * (term1 * drdm + term2);
            }
        }
        return out;
    }
};

template <>
struct ElasticHypersingular<2> {
    const double shear_modulus;
    const double poisson_ratio;
    const double trac_C2;
    typedef Vec2<double> OutType;
    typedef Vec2<double> InType;
    typedef Vec2<Vec2<double>> OperatorType;

    ElasticHypersingular(double shear_modulus, double poisson_ratio):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        trac_C2(1 - 2 * poisson_ratio)
    {}

    OperatorType operator()(double r2, const Vec2<double>& delta, 
                       const Vec2<double>& nsrc, const Vec2<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        const auto dr = delta / r;
        const auto drdn = dot_product(dr, nsrc);
        const auto drdm = dot_product(dr, nobs);
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                const auto line1 = 2 * drdn * (
                    trac_C2 * nobs[k] * dr[j] +
                    poisson_ratio * (nobs[j] * dr[k] + kronecker[k][j] * drdm) - 
                    4 * dr[k] * dr[j] * drdm);
                const auto line2 = trac_C2 * (
                    2 * nsrc[j] * dr[k] * drdm + kronecker[k][j] * dot_product(nsrc, nobs) 
                    + nsrc[k] * nobs[j]);
                const auto line3 = 2 * poisson_ratio * (
                        nsrc[k] * dr[j] * drdm + dot_product(nsrc, nobs) * dr[k] * dr[j]);
                const auto line4 = -(1 - 4 * poisson_ratio) * nsrc[j] * nobs[k];
                const auto C = (shear_modulus / (2 * M_PI * (1 - poisson_ratio) * r2));
                out[k][j] = C * (line1 + line2 + line3 + line4);
            }
        }
        return out;
    }
};

template <>
struct ElasticDisplacement<3> {
    const double disp_C1;
    const double disp_C2;
    typedef Vec3<double> OutType;
    typedef Vec3<double> InType;
    typedef Vec3<Vec3<double>> OperatorType;

    ElasticDisplacement(double shear_modulus, double poisson_ratio):
        disp_C1(1.0 / (16 * M_PI * shear_modulus * (1 - poisson_ratio))),
        disp_C2(3 - 4 * poisson_ratio)
    {}

    OperatorType operator()(double r2, const Vec3<double>& delta, 
                       const Vec3<double>& nsrc, const Vec3<double>& nobs) const {
        OperatorType out;
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                double r = std::sqrt(r2);
                out[k][j] = (disp_C1 / r) * (disp_C2 * kronecker[k][j] +
                                             delta[k] * delta[j] / r2);
            }
        }
        return out;
    }
};

template <>
struct ElasticTraction<3> {
    const double trac_C1;
    const double trac_C2;
    typedef Vec3<double> OutType;
    typedef Vec3<double> InType;
    typedef Vec3<Vec3<double>> OperatorType;

    ElasticTraction(double shear_modulus, double poisson_ratio):
        trac_C1(1.0 / (8 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
    {}
    
    OperatorType operator()(double r2, const Vec3<double>& delta, 
                       const Vec3<double>& nsrc, const Vec3<double>& nobs) const {
        OperatorType out;
        const double r = std::sqrt(r2);
        const auto drdn = dot_product(delta, nsrc) / r;
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                const auto term1 = (trac_C2 * kronecker[k][j] +
                                    3 * delta[k] * delta[j] / r2);
                const auto term2 = trac_C2 * (nsrc[j] * delta[k] -
                                              nsrc[k] * delta[j]) / r;
                out[k][j] = -(trac_C1 / r2) * (term1 * drdn - term2);
            }
        }
        return out;
    }
};

template <>
struct ElasticAdjointTraction<3> {
    const double trac_C1;
    const double trac_C2;
    typedef Vec3<Vec3<double>> OperatorType;

    ElasticAdjointTraction(double shear_modulus, double poisson_ratio):
        trac_C1(1.0 / (8 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
    {}
    
    OperatorType operator()(double r2, const Vec3<double>& delta, 
                       const Vec3<double>& nsrc, const Vec3<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        const auto drdm = dot_product(delta, nobs) / r;
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                const auto term1 = (trac_C2 * kronecker[k][j] +
                                    3 * delta[k] * delta[j] / r2);
                const auto term2 = trac_C2 * (nobs[j] * delta[k] -
                                              nobs[k] * delta[j]) / r;
                out[k][j] = (trac_C1 / r2) * (term1 * drdm + term2);
            }
        }
        return out;
    }
};

template <>
struct ElasticHypersingular<3> {
    const double poisson_ratio;
    const double trac_C2;
    const double hyp_C1;
    const double hyp_C2;
    const double hyp_C3;
    typedef Vec3<double> OutType;
    typedef Vec3<double> InType;
    typedef Vec3<Vec3<double>> OperatorType;

    ElasticHypersingular(double shear_modulus, double poisson_ratio):
        poisson_ratio(poisson_ratio),
        trac_C2(1 - 2 * poisson_ratio),
        hyp_C1(shear_modulus / (4 * M_PI * (1 - poisson_ratio))),
        hyp_C2(-1 + 4 * poisson_ratio),
        hyp_C3(3 * poisson_ratio)
    {}
    
    OperatorType operator()(double r2, const Vec3<double>& delta, 
                       const Vec3<double>& nsrc, const Vec3<double>& nobs) const {
        OperatorType out;
        double r = std::sqrt(r2);
        const Vec3<double> dr = delta / r;
        const auto drdm = dot_product(delta, nobs) / r;
        const auto drdn = dot_product(dr, nsrc);
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                const auto line1 = 3 * drdn * (
                    trac_C2 * nobs[k] * dr[j] +
                    poisson_ratio * (nobs[j] * dr[k] + kronecker[k][j] * drdm) -
                    5 * dr[k] * dr[j] * drdm);
                const auto line2 = trac_C2 * (
                    3 * nsrc[j] * dr[k] * drdm + kronecker[k][j] * dot_product(nsrc, nobs) 
                    + nsrc[k] * nobs[j]);
                const auto line3 = hyp_C3 * (
                        nsrc[k] * dr[j] * drdm + dot_product(nsrc, nobs) * dr[k] * dr[j]);
                const auto line4 = hyp_C2 * nsrc[j] * nobs[k];
                const auto C = hyp_C1 / (r2 * r);
                out[k][j] = C * (line1 + line2 + line3 + line4);
            }
        }
        return out;
    }
};

} // END namespace tbem
#endif
