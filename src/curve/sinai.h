#ifndef __SINAI_H
#define __SINAI_H

#include <cmath>
#include <random>
#include <utility>

/* x0 = sqrt(2+sqrt(3)) */
/* povrsina S = sqrt(2)*x0 - 2*pi/3*/
/* tezisce xt = yt = (sqrt(98+39*sqrt(3))-2*x0*pi) / (3*S) */
/* < r^2 > = (6*(8+7*sqrt(3))-8*(3+sqrt(3))*pi) / (3*S) */

namespace Sinai
{
    struct Circle : public Curve<Circle> {
        inline Derivatives derivatives (const Particle& p) const {
            static double x0 = sqrt (2.0 + sqrt (3.0));
            double x = p.x - x0;
            double y = p.y - x0;
            Derivatives d;
            d.f = x * x + y * y - 4.0;
            d.dfdx = 2.0 * x;
            d.dfdy = 2.0 * y;
            d.dfdt = 0.0; 
            return d;
        }
    };

    struct Xaxis : public Curve<Xaxis> {
        inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            d.f = p.y;
            d.dfdx = 0.0;
            d.dfdy = 1.0;
            d.dfdt = 0.0; 
            return d;
        }
    };

    struct Yaxis : public Curve<Yaxis> {
        inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            d.f = p.x;
            d.dfdx = 1.0;
            d.dfdy = 0.0;
            d.dfdt = 0.0; 
            return d;
        }
    };
}


#endif
