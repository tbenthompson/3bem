#ifndef __YUIQWEOIUQWE_KERNELS_H
#define __YUIQWEOIUQWE_KERNELS_H

#include <array>
#include <functional>

#include "vec.h"

template <typename T, typename F = double>
inline T one(const T& r2, 
             const Vec3<T>& delta,
             const Vec3<F>& nsrc,
             const Vec3<F>& nobs) {
    return 1.0;
}

template <typename T, typename F = double>
inline T laplace_single(const T& r2,
                        const Vec3<T>& delta,
                        const Vec3<F>& nsrc,
                        const Vec3<F>& nobs) {
    return 1.0 / (4.0 * M_PI * sqrt(r2));
}

template <typename T, typename F = double>
inline T laplace_double(const T& r2,
                        const Vec3<T>& delta,
                        const Vec3<F>& nsrc,
                        const Vec3<F>& nobs) {
    return -(nsrc[0] * delta[0] + nsrc[1] * delta[1] + nsrc[2] * delta[2]) / 
           (4 * M_PI * pow(r2, 1.5));
}

const double kronecker[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

class ElasticKernels {
public:
    ElasticKernels(double shear_modulus, double poisson_ratio):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        disp_C1(1.0 / (16 * M_PI * shear_modulus * (1 - poisson_ratio))),
        disp_C2(3 - 4 * poisson_ratio),
        trac_C1(1.0 / (8 * M_PI * (1 - poisson_ratio))),
        trac_C2(1 - 2 * poisson_ratio)
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
     */
    //TODO: precompute some of the constants in here.
    template <int k, int j>
    double hypersingular(double r2, 
                    const Vec3<double>& delta, 
                    const Vec3<double>& nsrc,
                    const Vec3<double>& nobs) const {
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
        const auto line3 = 3 * poisson_ratio * (
                nsrc[k] * dr[j] * drdm + dot(nsrc, nobs) * dr[k] * dr[j]);
        const auto line4 = -(1 - 4 * poisson_ratio) * nsrc[j] * nobs[k];
        const auto C = (shear_modulus / (4 * M_PI * (1 - poisson_ratio) * r2 * r));
        return C * (line1 + line2 + line3 + line4);
    }

    const double shear_modulus;
    const double poisson_ratio;
    const double disp_C1;
    const double disp_C2;
    const double trac_C1;
    const double trac_C2;
};

#endif
