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
    const double poisson_ratio;
    const double shear_modulus;
    // gravitational field multiplied by density of material
    const Vec<double,2> intensity;

    GravityDisplacement(double shear_modulus, double poisson_ratio,
            Vec<double,2> intensity):
        grav_C1(1.0 / (8 * M_PI * shear_modulus)),
        grav_C2(1.0 / (2 * (1 - poisson_ratio))),
        poisson_ratio(poisson_ratio),
        shear_modulus(shear_modulus),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec2<double>& delta, 
        const Vec2<double>& nobs, const Vec2<double>& nsrc) const 
    {
        (void)nobs;
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

    virtual std::unique_ptr<Kernel<2,2,2>> clone() const
    {
        return std::unique_ptr<Kernel<2,2,2>>(
            new GravityDisplacement<2>(shear_modulus, poisson_ratio, intensity)
        );
    }
};

template <>
struct GravityTraction<2>: public Kernel<2,2,2>
{
    const double shear_modulus;
    const double poisson_ratio;
    const Vec<double,2> intensity;

    GravityTraction(double shear_modulus, double poisson_ratio,
            Vec<double,2> intensity):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec2<double>& delta, 
        const Vec2<double>& nobs, const Vec2<double>& nsrc) const 
    {
        typename Kernel::OperatorType out;
        double r = std::sqrt(r2);
        double logr = std::log(r);
        auto deltaxy = delta[0] * delta[1];
        auto deltax2 = delta[0] * delta[0];
        auto deltay2 = delta[1] * delta[1];
        auto denom = (16*M_PI*(-poisson_ratio + 1)*(2*poisson_ratio - 1)*r2);
        double C1 = r2 * (0.25 * logr + 0.125);
        out[0][0] = (nobs[0]*(16.0*poisson_ratio*(-intensity[0]*(2*nsrc[0]*(-poisson_ratio + 1)*(0.25*deltax2 + C1) - 1.0*nsrc[0]*(0.25*deltax2 + C1) + 0.5*nsrc[1]*(-poisson_ratio + 1)*deltaxy) - intensity[1]*(0.5*nsrc[0]*(-poisson_ratio + 1)*deltaxy + 2*nsrc[1]*(-poisson_ratio + 1)*(0.25*deltay2 + C1) - 1.0*nsrc[1]*(0.25*deltay2 + C1)) + 0.25*deltaxy*(intensity[0]*nsrc[1] + intensity[1]*nsrc[0])) + 4*(2*poisson_ratio - 1)*(4*intensity[0]*(2*nsrc[0]*(-poisson_ratio + 1)*(0.25*deltax2 + C1) - 1.0*nsrc[0]*(0.25*deltax2 + C1) + 0.5*nsrc[1]*(-poisson_ratio + 1)*deltaxy) - 1.0*intensity[1]*nsrc[0]*deltaxy)) - 8*nobs[1]*(2*poisson_ratio - 1)*(intensity[0]*nsrc[1]*(0.25*deltax2 + C1) - intensity[0]*(0.5*nsrc[0]*(-poisson_ratio + 1)*deltaxy - 0.25*nsrc[0]*deltaxy + 2*nsrc[1]*(-poisson_ratio + 1)*(0.25*deltay2 + C1)) + intensity[1]*nsrc[0]*(0.25*deltay2 + C1) - intensity[1]*(2*nsrc[0]*(-poisson_ratio + 1)*(0.25*deltax2 + C1) + 0.5*nsrc[1]*(-poisson_ratio + 1)*deltaxy - 0.25*nsrc[1]*deltaxy)))/denom;
        out[1][0] = (-8*nobs[0]*(2*poisson_ratio - 1)*(intensity[0]*nsrc[1]*(0.25*deltax2 + C1) - intensity[0]*(0.5*nsrc[0]*(-poisson_ratio + 1)*deltaxy - 0.25*nsrc[0]*deltaxy + 2*nsrc[1]*(-poisson_ratio + 1)*(0.25*deltay2 + C1)) + intensity[1]*nsrc[0]*(0.25*deltay2 + C1) - intensity[1]*(2*nsrc[0]*(-poisson_ratio + 1)*(0.25*deltax2 + C1) + 0.5*nsrc[1]*(-poisson_ratio + 1)*deltaxy - 0.25*nsrc[1]*deltaxy)) - nobs[1]*(-16.0*poisson_ratio*(-intensity[0]*(2*nsrc[0]*(-poisson_ratio + 1)*(0.25*deltax2 + C1) - 1.0*nsrc[0]*(0.25*deltax2 + C1) + 0.5*nsrc[1]*(-poisson_ratio + 1)*deltaxy) - intensity[1]*(0.5*nsrc[0]*(-poisson_ratio + 1)*deltaxy + 2*nsrc[1]*(-poisson_ratio + 1)*(0.25*deltay2 + C1) - 1.0*nsrc[1]*(0.25*deltay2 + C1)) + 0.25*deltaxy*(intensity[0]*nsrc[1] + intensity[1]*nsrc[0])) + 4*(2*poisson_ratio - 1)*(intensity[0]*nsrc[1]*deltaxy - 4*intensity[1]*(0.5*nsrc[0]*(-poisson_ratio + 1)*deltaxy + 2*nsrc[1]*(-poisson_ratio + 1)*(0.25*deltay2 + C1) - 1.0*nsrc[1]*(0.25*deltay2 + C1)))))/denom;
        out[0][1] = 0.0;
        out[1][1] = 0.0;
        return out;
    }

    virtual std::unique_ptr<Kernel<2,2,2>> clone() const
    {
        return std::unique_ptr<Kernel<2,2,2>>(
            new GravityTraction<2>(shear_modulus, poisson_ratio, intensity)
        );
    }
};

template <>
struct GravityDisplacement<3>: public Kernel<3,3,3>
{
    const double poisson_ratio;
    const double shear_modulus;
    // gravitational field multiplied by density of material
    const Vec<double,3> intensity;

    GravityDisplacement(double shear_modulus, double poisson_ratio,
            Vec<double,3> intensity):
        poisson_ratio(poisson_ratio),
        shear_modulus(shear_modulus),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec3<double>& delta, 
        const Vec3<double>& nobs, const Vec3<double>& nsrc) const 
    {
        (void)r2; (void)delta; (void)nobs; (void)nsrc;
        return zeros<Kernel::OperatorType>::make();
    }

    virtual std::unique_ptr<Kernel<3,3,3>> clone() const
    {
        return std::unique_ptr<Kernel<3,3,3>>(
            new GravityDisplacement<3>(shear_modulus, poisson_ratio, intensity)
        );
    }
};

template <>
struct GravityTraction<3>: public Kernel<3,3,3>
{
    const double shear_modulus;
    const double poisson_ratio;
    const Vec<double,3> intensity;

    GravityTraction(double shear_modulus, double poisson_ratio,
            Vec<double,3> intensity):
        shear_modulus(shear_modulus),
        poisson_ratio(poisson_ratio),
        intensity(intensity)
    {}

    typename Kernel::OperatorType call(double r2, const Vec3<double>& delta, 
        const Vec3<double>& nobs, const Vec3<double>& nsrc) const 

    {
        (void)r2; (void)delta; (void)nobs; (void)nsrc;
        return zeros<Kernel::OperatorType>::make();
    }

    virtual std::unique_ptr<Kernel<3,3,3>> clone() const
    {
        return std::unique_ptr<Kernel<3,3,3>>(
            new GravityTraction<3>(shear_modulus, poisson_ratio, intensity)
        );
    }
};

}

#endif
