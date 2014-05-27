#ifndef __SINAI2_H
#define __SINAI2_H

#include <cmath>
#include <random>
#include <utility>

namespace Sinai2
{
    inline Derivatives circle_derivatives (const double a, const Particle& p)
    {
        double x = p.x;
        double y = p.y - 2.0 - a;
        Derivatives d;
        d.f = x * x + y * y - 4.0;
        d.dfdx = 2.0 * x;
        d.dfdy = 2.0 * y;
        d.dfdt = 0.0; 
        return d;
    }

    struct Circle {
        Circle (double a_) : a(a_) {};
        inline Derivatives derivatives (const Particle& p) const
                {return circle_derivatives (a, p);}
        private :
        double a;    
    };

    inline Derivatives xaxis_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = p.y;
        d.dfdx = 0.0;
        d.dfdy = 1.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Xaxis {
        inline Derivatives derivatives (const Particle& p) const
                {return xaxis_derivatives (p);}
    };

    inline Derivatives vleft_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = p.x + 1;
        d.dfdx = 1.0;
        d.dfdy = 0.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Vleft {
        inline Derivatives derivatives (const Particle& p) const
                {return vleft_derivatives (p);}
    };

    inline Derivatives vright_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = -p.x + 1;
        d.dfdx = -1.0;
        d.dfdy = 0.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Vright {
        inline Derivatives derivatives (const Particle& p) const
                {return vright_derivatives (p);}
    };
}


#endif
