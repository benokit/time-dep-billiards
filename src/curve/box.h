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

#ifndef __BOX_H
#define __BOX_H

#include <cmath>
#include <random>
#include <utility>
#include "../billiard.h"

namespace Box
{
    inline Derivatives up_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = -p.y + 1.0;
        d.dfdx = 0.0;
        d.dfdy = 1.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Up {
        inline Derivatives derivatives (const Particle& p) const
                {return up_derivatives (p);}
    };

    inline Derivatives down_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = p.y;
        d.dfdx = 0.0;
        d.dfdy = 1.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Down {
        inline Derivatives derivatives (const Particle& p) const
                {return down_derivatives (p);}
    };

    inline Derivatives left_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = p.x + 1;
        d.dfdx = 1.0;
        d.dfdy = 0.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Left {
        inline Derivatives derivatives (const Particle& p) const
                {return left_derivatives (p);}
    };

    inline Derivatives right_derivatives (const Particle& p)
    {
        Derivatives d;
        d.f = -p.x + 1;
        d.dfdx = -1.0;
        d.dfdy = 0.0;
        d.dfdt = 0.0; 
        return d;
    }

    struct Right {
        inline Derivatives derivatives (const Particle& p) const
                {return right_derivatives (p);}
    };

    Particle rand_particle ()
    {
        static std::default_random_engine generator;
        static std::uniform_real_distribution<double> distribution (0, 1);
        Particle p;
        p.x = 2.0 * (distribution (generator) - 0.5);
        p.y = distribution (generator);
        double phi = 2.0 * M_PI * distribution (generator);
        p.vx = cos (phi);
        p.vy = sin (phi);
        return p;
    }
}


#endif

