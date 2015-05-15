#ifndef __GGGGGGGGJJJJKKKZZZZZZ_GRAVITY_KERNELS_H
#define __GGGGGGGGJJJJKKKZZZZZZ_GRAVITY_KERNELS_H

namespace tbem {

/* Gravity is implemented for a constant gravitational field via the volume
 * to surface integral transformation demonstrated in 
 * "A boundary element formulation of problems in linear isotropic elasticity
 * with body forces" -- D. J. Danson 1981
 */
template <size_t dim>
struct GravityDisplacement;
template <size_t dim>
struct GravityTraction;

template <>
struct GravityDisplacement<2>: public Kernel<2,2,2>
{
    const double grav_C1;
    const double grav_C2;
    // gravitational field multiplied by density of material
    const Vec<double,2> intensity;

    GravityDisplacement(double shear_modulus, double poisson_ratio,
            Vec<double,2> intensity):
        grav_C1(1.0 / (8 * M_PI * shear_modulus)),
        grav_C2(1.0 / (2 * (1 - poisson_ratio))),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec2<double>& delta, 
        const Vec2<double>& nobs, const Vec2<double>& nsrc) const 
    {
        typename Kernel::OperatorType out;
        double r = std::sqrt(r2);
        const auto dr = delta / r;
        auto drdn = dot_product(dr, nsrc);
        auto drdg = dot_product(dr, intensity);
        for (int k = 0; k < 2; k++) {
            out[k][0] = grav_C1 * r * (-2 * std::log(r) - 1) * 
                (intensity[k] * drdn - (nsrc[k] * drdg) * grav_C2);
            out[k][1] = 0.0;
        }
        return out;
    }
};
template <>
struct GravityDisplacement<3>: public Kernel<3,3,3>
{
    // gravitational field multiplied by density of material
    const Vec<double,3> intensity;

    GravityDisplacement(double shear_modulus, double poisson_ratio,
            Vec<double,3> intensity):
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec3<double>& delta, 
        const Vec3<double>& nobs, const Vec3<double>& nsrc) const 
    {
    }
};

template <>
struct GravityTraction<2>: public Kernel<2,2,2>
{
    const double grav_C1;
    const double grav_C2;
    const double grav_C3;
    const double grav_C4;
    const double poisson_ratio;
    const Vec<double,2> intensity;

    GravityTraction(double shear_modulus, double poisson_ratio,
            Vec<double,2> intensity):
        grav_C1(1.0 / (2 * M_PI)),
        grav_C2(poisson_ratio / (1 - poisson_ratio)),
        grav_C3(1.0 / (2 * (1 - poisson_ratio))),
        grav_C4((1.0 - 2 * poisson_ratio) / 2),
        poisson_ratio(poisson_ratio),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec2<double>& delta, 
        const Vec2<double>& nobs, const Vec2<double>& nsrc) const 
    {
        typename Kernel::OperatorType out;
        double r = std::sqrt(r2);
        const auto dr = delta / r;
        auto drdn = dot_product(dr, nsrc);
        auto drdg = dot_product(dr, intensity);
        auto drdm = dot_product(dr, nobs);
        auto ndotm = dot_product(nobs, nsrc);
        auto bdotm = dot_product(intensity, nobs);
        auto bdotn = dot_product(intensity, nsrc);
        for (int k = 0; k < 2; k++) {
            out[k][0] = grav_C1 * (
                (1 + std::log(r)) * (
                    drdn * (intensity[k] * drdm + bdotm * dr[k] +
                            grav_C2 * drdg * nobs[k]) -
                    drdg * grav_C3 * (nsrc[k] * drdm + ndotm * dr[k])) -
                grav_C3 * (0.5 + std::log(r)) * (
                    grav_C4 * (intensity[k] * ndotm + bdotm * nsrc[k])
                        + poisson_ratio * bdotn * nobs[k]));
            out[k][1] = 0.0;
        }
        return out;
    }
};

template <>
struct GravityTraction<3>: public Kernel<3,3,3>
{
    const Vec<double,3> intensity;

    GravityTraction(double shear_modulus, double poisson_ratio,
            Vec<double,3> intensity):
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec3<double>& delta, 
        const Vec3<double>& nobs, const Vec3<double>& nsrc) const 

    {
    }
};

}

#endif
