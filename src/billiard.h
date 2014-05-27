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

#ifndef __BILLIARD_H
#define __BILLIARD_H

#include <cmath>
#include <iomanip>
#include <tuple>
#include "froot.h"

struct Particle {
    double x;
    double y;
    double vx;
    double vy;
    double t;
    inline double velocity () const {return sqrt (vx * vx + vy * vy);}
    inline void set_velocity (double v0) {double f = v0 / velocity(); vx *= f; vy *= f;}
    inline void set_time (double t0) {t = t0;} 
    inline void time_modulo (double period) {t = fmod (t, period);}
    template<typename S> inline void print (S&);
    inline void print () {print (std::cout);}
};

template<typename S>
void Particle::print (S& stream) {
    stream << std::setw(25) << std::setprecision(16) << x; 
    stream << std::setw(25) << std::setprecision(16) << y; 
    stream << std::setw(25) << std::setprecision(16) << vx; 
    stream << std::setw(25) << std::setprecision(16) << vy; 
    stream << std::setw(25) << std::setprecision(16) << t; 
    stream << std::endl;
}

struct FreeFlight {
    inline Particle operator () (const Particle& p, double dt) const {
        return (Particle) 
            {p.x + p.vx * dt, p.y + p.vy * dt, p.vx, p.vy, p.t + dt};
    }
};

struct Periodic {
    inline void operator () (Particle& p) const {p.time_modulo (2 * M_PI);}
};

struct Aperiodic {
    inline void operator () (const Particle& p) const {}
};

struct Static {
    inline void operator () (Particle& p) const {p.set_time (0);}
};

////////////////////////////////////////////////////////////////////////////////

struct Derivatives {
    double f;
    double dfdx;
    double dfdy;
    double dfdt;
};

template <typename C>
struct Curve {
    inline void fdf (const Particle&, double&, double&) const;
    inline void reflection (Particle&) const;
    inline double tangent_velocity (const Particle&) const;
};

template <typename C>
inline void Curve<C>::fdf (const Particle& p, double& f, double& df) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    f  = d.f;
    df = d.dfdx * p.vx + d.dfdy * p.vy + d.dfdt;
}

template <typename C>
inline void Curve<C>::reflection (Particle& p) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    double kick = 2.0 * (d.dfdx * p.vx + d.dfdy * p.vy + d.dfdt) 
                      / (d.dfdx * d.dfdx + d.dfdy * d.dfdy);
    p.vx -= kick * d.dfdx;
    p.vy -= kick * d.dfdy;
}

template <typename C>
inline double Curve<C>::tangent_velocity (const Particle& p) const
{
    Derivatives d = static_cast<const C*>(this) -> derivatives (p);
    return (d.dfdx * p.vy - d.dfdy * p.vx) 
           / sqrt (d.dfdx * d.dfdx + d.dfdy * d.dfdy);
}

////////////////////////////////////////////////////////////////////////////////

struct GenericTimeStep {
    GenericTimeStep () : geometricScale(0.1), timeScale(0.1), vmin(0.01) {}
    GenericTimeStep (double s, double t, double v) : 
        geometricScale(s), timeScale(t), vmin(v) {}
    inline double operator () (const Particle& p) const {
        double v = p.velocity();
        if (v > vmin) return v < (geometricScale / timeScale) ? timeScale : (geometricScale / v);
        else return (geometricScale / v);
    }
    private:
    const double geometricScale, timeScale, vmin;
};

////////////////////////////////////////////////////////////////////////////////

template <typename F, typename Z, typename T, typename ...Cs>
class Billiard {
    public:
        inline void collision (Particle& p) const {
            base_collision (p, typename genseq<sizeof...(Cs)>::type());
        }
        inline bool is_inside (const Particle& p) const {
            return base_is_inside (p, typename genseq<sizeof...(Cs)>::type());
        }
        T time_fold;
        F fly;
    private:
        Z time_step;
        std::tuple<Cs...> curves;

        template<int ...>  struct seq {};
        template<int N, int ...S> struct genseq : genseq<N-1, N-1, S...> {};
        template<int ...S> struct genseq<0, S...>{ typedef seq<S...> type; };
        
        template<int ...S>
        inline void base_collision (Particle&, seq<S...>) const;

        template<int ...S>
        inline bool base_is_inside (const Particle&, seq<S...>) const;

};

template <typename F, typename Z, typename T, typename ...Cs> 
template <int ...S>
inline bool Billiard<F,Z,T,Cs...>::base_is_inside (const Particle& p, seq<S...>) const
{
    bool isInside = true;
    is_inside_aux (p, isInside, std::get<S>(curves) ...);
    return isInside;
}


template <typename F, typename Z, typename T, typename ...Cs>
template <int ...S>
inline void Billiard<F,Z,T,Cs...>::base_collision (Particle& p, seq<S...>) const
{
    double step = time_step (p);
    double ta, tb{0.0};

    Particle p0 = p;

    bool isCollision = false;

    while (!isCollision) {
        ta = tb;
        tb = tb + step;
        is_collision_aux (p0, ta, tb, p, isCollision, fly, std::get<S>(curves) ...);
    }
}

template<typename F>
static inline void is_collision_aux (Particle p, double ta, double tb, Particle& p1, 
                              bool& isCollision, const F& fly) {}

template<typename F, typename C, typename... Cs>
static inline void is_collision_aux (Particle p, double ta, double tb, Particle& p1, bool& isCollision, 
                              const F& fly, const C& curve, const Cs&... curves) 
{
    double tm;
    auto f = [&fly, &curve, &p] (double t, double& f, double& df) {curve.fdf (fly (p, t), f, df);};
    if (find_next_root (f, ta, tb, tm)) {
        p1 = fly (p, tm);
        curve.reflection (p1);
        isCollision = true;
        is_collision_aux (p, ta, tm, p1, isCollision, fly, curves...);
    }
    else {
        is_collision_aux (p, ta, tb, p1, isCollision, fly, curves...);
    }
}

static inline void is_inside_aux (const Particle& p, bool& isInside) {}

template<typename C, typename... Cs>
static inline void is_inside_aux (const Particle& p, bool& isInside, const C& curve, const Cs&... curves) 
{
    double f, df;
    curve.fdf (p, f, df);
    isInside = isInside && f > 0;
    is_inside_aux (p, isInside, curves...);
}

#endif
