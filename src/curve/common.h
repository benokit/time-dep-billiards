#ifndef __COMMON_H
#define __COMMON_H

#include <cmath>

namespace Common
{
    inline Derivatives xAxisDerivatives (const Particle& p)
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
                {return xAxisDerivatives (p);}
    };
}


#endif

