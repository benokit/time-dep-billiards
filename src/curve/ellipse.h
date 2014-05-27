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

#ifndef __ELLIPSE_H
#define __ELLIPSE_H

class Ellipse : public Curve<Ellipse> {
    public:
        Ellipse (double bb) : b(bb) {}
        inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            d.f = 1.0 - p.x * p.x - b * p.y * p.y;
            d.dfdx = -2.0 * p.x;
            d.dfdy = -2.0 * b * p.y;
            d.dfdt = 0.0; 
            return d;
        }
    private:
        const double b;
};

#endif
