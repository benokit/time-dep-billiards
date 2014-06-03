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

class Robnik : public Curve<Robnik> {
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
