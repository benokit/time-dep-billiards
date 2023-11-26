#ifndef __DOMAIN_H
#define __DOMAIN_H

#include "billiard.h"

struct Derivatives {
    double f;
    double dfdx;
    double dfdy;
    double dfdt;
};

template <typename C>
struct Domain {
    inline void fdf (const Particle&, double&, double&) const;
    inline void reflection (Particle&) const;
    inline double tangent_velocity (const Particle&) const;
};

template <typename C>
inline void Domain<C>::fdf (const Particle& p, double& f, double& df) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    f  = d.f;
    df = d.dfdx * p.vx + d.dfdy * p.vy + d.dfdt;
}

template <typename C>
inline void Domain<C>::reflection (Particle& p) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    double kick = 2.0 * (d.dfdx * p.vx + d.dfdy * p.vy + d.dfdt) 
                      / (d.dfdx * d.dfdx + d.dfdy * d.dfdy);
    p.vx -= kick * d.dfdx;
    p.vy -= kick * d.dfdy;
}

template <typename C>
inline double Domain<C>::tangent_velocity (const Particle& p) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    return (d.dfdx * p.vy - d.dfdy * p.vx) 
           / sqrt (d.dfdx * d.dfdx + d.dfdy * d.dfdy);
}

#endif