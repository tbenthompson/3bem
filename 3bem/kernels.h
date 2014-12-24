#ifndef __YUIQWEOIUQWE_KERNELS_H
#define __YUIQWEOIUQWE_KERNELS_H

#include <array>
#include <functional>

#include "vec.h"

namespace tbem { 

//TODO: Think about separating each kernel into its own class and removing
//some of this fluff.
#define wrap_K2(KNAME, K, J) \
    [&] (const double& r2, const Vec2<double>& delta,\
         const Vec2<double> nsrc, const Vec2<double>& nobs) {\
        return this->KNAME<K,J>(r2, delta, nsrc, nobs);\
    }

#define wrap_K3(KNAME, K, J) \
    [&] (const double& r2, const Vec3<double>& delta,\
         const Vec3<double> nsrc, const Vec3<double>& nobs) {\
        return this->KNAME<K,J>(r2, delta, nsrc, nobs);\
    }

#define KERNELMAT2(KNAME, ARRNAME) \
    const std::array<std::array<Kernel<2>,2>,2> ARRNAME = {{ {{\
            wrap_K2(KNAME, 0, 0), wrap_K2(KNAME, 0, 1)\
        }}, {{\
            wrap_K2(KNAME, 1, 0), wrap_K2(KNAME, 1, 1)\
        }}\
    }}

#define KERNELMAT3(KNAME, ARRNAME) \
    const std::array<std::array<Kernel<3>,3>,3> ARRNAME = {{ {{\
            wrap_K3(KNAME, 0, 0), wrap_K3(KNAME, 0, 1), wrap_K3(KNAME, 0, 2) \
        }}, {{\
            wrap_K3(KNAME, 1, 0), wrap_K3(KNAME, 1, 1), wrap_K3(KNAME, 1, 2)\
        }}, {{\
            wrap_K3(KNAME, 2, 0), wrap_K3(KNAME, 2, 1), wrap_K3(KNAME, 2, 2)\
        }}\
    }}

template <int dim>
using Kernel = std::function<double(const double&, const Vec<double,dim>&,
                             const Vec<double,dim>&, const Vec<double,dim>&)>;

template <int dim>
inline double one(const double& r2, 
             const Vec<double,dim>& delta,
             const Vec<double,dim>& nsrc,
             const Vec<double,dim>& nobs) {
    return 1.0;
}


/* Laplace/Poisson equation kernels. */
template <int dim>
inline double laplace_single(const double& r2, const Vec<double,dim>& delta,
                        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs);
template <int dim>
inline double laplace_double(const double& r2, const Vec<double,dim>& delta,
                        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs);
template <int dim>
inline double laplace_hypersingular(const double& r2, const Vec<double,dim>& delta,
                        const Vec<double,dim>& nsrc, const Vec<double,dim>& nobs);

template <>
inline double laplace_single<3>(const double& r2, const Vec3<double>& delta,
                        const Vec3<double>& nsrc, const Vec3<double>& nobs) {
    return 1.0 / (4.0 * M_PI * std::sqrt(r2));
}

template <>
inline double laplace_double<3>(const double& r2, const Vec3<double>& delta,
                        const Vec3<double>& nsrc, const Vec3<double>& nobs) {
    return dot(nsrc, delta) / (4 * M_PI * r2 * std::sqrt(r2));
}

template <>
inline double laplace_single<2>(const double& r2, const Vec2<double>& delta,
                        const Vec2<double>& nsrc, const Vec2<double>& nobs) {
    return std::log(std::sqrt(r2)) / (2 * M_PI);
}

template <>
inline double laplace_double<2>(const double& r2, const Vec2<double>& delta,
                        const Vec2<double>& nsrc, const Vec2<double>& nobs) {
    return dot(nsrc, delta) / (2 * M_PI * r2);
}

template <>
inline double laplace_hypersingular<2>(const double& r2, const Vec2<double>& delta,
                        const Vec2<double>& nsrc, const Vec2<double>& nobs) {
    return ((-dot(nobs, nsrc) / r2) + 
            ((2 * dot(nsrc, delta) * dot(nobs, delta)) / (r2 * r2)))
        / (2 * M_PI);
}

/* 3D linear isotropic elastic kernels. */
const double kronecker[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

template <int dim>
struct ElasticKernels;

template <>
struct ElasticKernels<3> {
    ElasticKernels(double shear_modulus, double poisson_ratio):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        disp_C1(1.0 / (16 * M_PI * shear_modulus * (1 - poisson_ratio))),
        disp_C2(3 - 4 * poisson_ratio),
        trac_C1(1.0 / (8 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio),
        hyp_C1(shear_modulus / (4 * M_PI * (1 - poisson_ratio))),
        hyp_C2(-1 + 4 * poisson_ratio),
        hyp_C3(3 * poisson_ratio)
    {}

    /* The displacement kernel U* taken from the SGBEM book by 
     * Sutradhar, Paulino, Gray
     */
    template <int k, int j>
    double displacement(double r2, 
                        const Vec3<double>& delta, 
                        const Vec3<double>& nsrc,
                        const Vec3<double>& nobs) const {
        double r = std::sqrt(r2);
        return (disp_C1 / r) *
               (disp_C2 * kronecker[k][j] + delta[k] * delta[j] / r2);
    }

    /* Traction kernel derived by applying the traction operator
     * to the displacement kernel w.r.t the source coords. The form here
     * is from the SGBEM book.
     */
    template <int k, int j>
    double traction(double r2, 
                    const Vec3<double>& delta, 
                    const Vec3<double>& nsrc,
                    const Vec3<double>& nobs) const {
        double r = std::sqrt(r2);
        const auto term1 = (trac_C2 * kronecker[k][j] + 3 * delta[k] * delta[j] / r2);
        const auto drdn = dot(delta, nsrc) / r;
        const auto term2 = trac_C2 * (nsrc[j] * delta[k] - nsrc[k] * delta[j]) / r;
        return -(trac_C1 / r2) * (term1 * drdn - term2);
    }

    /* Adjoint traction kernel derived by applying the traction operator
     * to the displacement kernel w.r.t. the observation coords. From the
     * SGBEM book multiplied by the observation normal vector.
     * Take the expression from the SGBEM book and multiply by the obs normal:
     * D_{lkm} * n_{obs,l} 
     * Then replace m with j to get the formula below.
     */
    template <int k, int j>
    double adjoint_traction(double r2, 
                    const Vec3<double>& delta, 
                    const Vec3<double>& nsrc,
                    const Vec3<double>& nobs) const {
        double r = std::sqrt(r2);
        const auto term1 = (trac_C2 * kronecker[k][j] + 3 * delta[k] * delta[j] / r2);
        const auto drdm = dot(delta, nobs) / r;
        const auto term2 = trac_C2 * (nobs[j] * delta[k] - nobs[k] * delta[j]) / r;
        return (trac_C1 / r2) * (term1 * drdm + term2);
    }

    /* Hypersingular kernel derived by applying the traction operator twice
     * to the displacement kernel w.r.t. both the observation coords and
     * the source coords. From the SGBEM book multiplied by the observation
     * normal vector. 
     * Take the expression from the SGBEM book and multiply by the obs normal:
     * S_{lkm} * n_{obs,l} 
     * Then replace m with j to get the formula below.
     */
    template <int k, int j>
    double hypersingular(const double& r2, 
                    const Vec3<double>& delta, 
                    const Vec3<double>& nsrc,
                    const Vec3<double>& nobs) {
        double r = std::sqrt(r2);
        const Vec3<double> dr = delta / r;
        const auto drdn = dot(dr, nsrc);
        const auto drdm = dot(dr, nobs);
        const auto line1 = 3 * drdn * (
            trac_C2 * nobs[k] * dr[j] +
            poisson_ratio * (nobs[j] * dr[k] + kronecker[k][j] * drdm) -
            5 * dr[k] * dr[j] * drdm);
        const auto line2 = trac_C2 * (
            3 * nsrc[j] * dr[k] * drdm + kronecker[k][j] * dot(nsrc, nobs) 
            + nsrc[k] * nobs[j]);
        const auto line3 = hyp_C3 * (
                nsrc[k] * dr[j] * drdm + dot(nsrc, nobs) * dr[k] * dr[j]);
        const auto line4 = hyp_C2 * nsrc[j] * nobs[k];
        const auto C = hyp_C1 / (r2 * r);
        return C * (line1 + line2 + line3 + line4);
    }
    
    KERNELMAT3(hypersingular, hypersingular_mat);
    KERNELMAT3(displacement, displacement_mat);
    KERNELMAT3(traction, traction_mat);
    KERNELMAT3(adjoint_traction, adjoint_traction_mat);

    const double shear_modulus;
    const double poisson_ratio;
    const double disp_C1;
    const double disp_C2;
    const double trac_C1;
    const double trac_C2;
    const double hyp_C1;
    const double hyp_C2;
    const double hyp_C3;
};

/* See docs for ElasticKernels<3>
 */
template <>
struct ElasticKernels<2> {
    ElasticKernels(double shear_modulus, double poisson_ratio):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        disp_C1(1.0 / (8 * M_PI * shear_modulus * (1 - poisson_ratio))),
        disp_C2(3 - 4 * poisson_ratio),
        trac_C1(1.0 / (4 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
    {}

    template <int k, int j>
    double displacement(double r2, 
                        const Vec2<double>& delta, 
                        const Vec2<double>& nsrc,
                        const Vec2<double>& nobs) const {
        double r = std::sqrt(r2);
        return disp_C1 *
               (-disp_C2 * kronecker[k][j] * std::log(r)  + 
                delta[k] * delta[j] / r2);
    }

    template <int k, int j>
    double traction(double r2, 
                    const Vec2<double>& delta, 
                    const Vec2<double>& nsrc,
                    const Vec2<double>& nobs) const {
        double r = std::sqrt(r2);
        const auto drdn = dot(delta, nsrc) / r;
        const auto term1 = (trac_C2 * kronecker[k][j] + 2 * delta[k] * delta[j] / r2);
        const auto term2 = trac_C2 * (nsrc[j] * delta[k] - nsrc[k] * delta[j]) / r;
        return -(trac_C1 / r2) * (term1 * drdn - term2);
    }

    template <int k, int j>
    double adjoint_traction(double r2, 
                    const Vec2<double>& delta, 
                    const Vec2<double>& nsrc,
                    const Vec2<double>& nobs) const {
        double r = std::sqrt(r2);
        const auto drdm = dot(delta, nobs) / r;
        const auto term1 = (trac_C2 * kronecker[k][j] + 2 * delta[k] * delta[j] / r2);
        const auto term2 = trac_C2 * (nobs[j] * delta[k] - nobs[k] * delta[j]) / r;
        return (trac_C1 / r2) * (term1 * drdm + term2);
    }

    template <int k, int j>
    double hypersingular(const double& r2, 
                    const Vec2<double>& delta, 
                    const Vec2<double>& nsrc,
                    const Vec2<double>& nobs) {
        double r = std::sqrt(r2);
        const auto dr = delta / r;
        const auto drdn = dot(dr, nsrc);
        const auto drdm = dot(dr, nobs);
        const auto line1 = 2 * drdn * (
            trac_C2 * nobs[k] * dr[j] +
            poisson_ratio * (nobs[j] * dr[k] + kronecker[k][j] * drdm) - 
            4 * dr[k] * dr[j] * drdm);
        const auto line2 = trac_C2 * (
            2 * nsrc[j] * dr[k] * drdm + kronecker[k][j] * dot(nsrc, nobs) 
            + nsrc[k] * nobs[j]);
        const auto line3 = 2 * poisson_ratio * (
                nsrc[k] * dr[j] * drdm + dot(nsrc, nobs) * dr[k] * dr[j]);
        const auto line4 = -(1 - 4 * poisson_ratio) * nsrc[j] * nobs[k];
        const auto C = (shear_modulus / (2 * M_PI * (1 - poisson_ratio) * r2));
        return C * (line1 + line2 + line3 + line4);
    }

    KERNELMAT2(displacement, displacement_mat);
    KERNELMAT2(traction, traction_mat);
    KERNELMAT2(adjoint_traction, adjoint_traction_mat);
    KERNELMAT2(hypersingular, hypersingular_mat);

    const double shear_modulus;
    const double poisson_ratio;
    const double disp_C1;
    const double disp_C2;
    const double trac_C1;
    const double trac_C2;
};

} // END namespace tbem
#endif
