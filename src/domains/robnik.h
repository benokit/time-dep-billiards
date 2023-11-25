#ifndef __ROBNIK_H
#define __ROBNIK_H

#include "../billiard.h"

class Robnik : public Domain<Robnik> {
    public:
        Robnik (double lam_) : lam(lam_) {}
        inline Derivatives derivatives (const Particle& p) const {
            double w, dwdx, dwdy;
            Derivatives d;
            w  = p.x * p.x + p.y * p.y - lam * lam;
            dwdx = 2.0 * p.x;
            dwdy = 2.0 * p.y;
            d.f = - w * w + w + 2.0 * lam * (lam + p.x);
            d.dfdx = -2.0 * w * dwdx + dwdx + 2.0 * lam;
            d.dfdy = -2.0 * w * dwdy + dwdy;
            d.dfdt = 0.0;
            return d;
        }
    private:
        const double lam;
};

#endif
