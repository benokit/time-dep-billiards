/*
  Copyright (C) 2014 Benjamin Batistic

  This file is part of Billiards Numerical Library.

  Billiards Numerical Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Billiards Numerical Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Billiards Numerical Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ROBNIK_H
#define __ROBNIK_H

#include <cmath>
#include <random>

namespace Robnik
{
    inline Derivatives derivatives_ (double lam, const Particle& p)
    {
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

    Particle rand_particle (double lam)
    {
        static std::default_random_engine generator;
        static std::uniform_real_distribution<double> distribution (0, 1);
        Particle p;
        double xmin, xmax, ymax;
        xmax = lam + 1.0;
        xmin = lam > 0.25 ? -(lam + 1.0 / (8.0 * lam)) : lam - 1.0;
        ymax = (1.0 - 2.0 * lam * lam) * (1.0 + 4.0 * lam * lam);
        Derivatives d;
        do {
            p.x = distribution (generator);
            p.x = xmin + (xmax - xmin) * (p.x);
            p.y = distribution (generator);
            p.y = ymax * (p.y);
            d  = derivatives_ (lam, p);
        } while (d.f <= 0.0);
        double phi = 2 * M_PI * distribution (generator);
        p.vx = cos (phi);
        p.vy = sin (phi);
        return p;
    }

    class Curve {
        public:
            Curve (double lam) : lambda(lam) {}
            inline Derivatives derivatives (const Particle& p) const
                {return derivatives_ (lambda, p);}
        private:
            const double lambda;
    };
}

#endif
