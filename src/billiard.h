#ifndef __BILLIARD_H
#define __BILLIARD_H

#include <cmath>
#include <iomanip>
#include <tuple>
#include <iostream>
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

////////////////////////////////////////////////////////////////////////////////

struct ConstantTimeScale {
    ConstantTimeScale () : time_scale(1.0) {}
    ConstantTimeScale (double t) : time_scale(t) {}
    inline double operator () (const Particle& p) const {
        return time_scale;
    }
    private:
    const double time_scale;
};

struct AdaptiveTimeScale {
    AdaptiveTimeScale () : geometric_scale(0.1), time_scale(0.1), too_slow_velocity(0.01) {}
    AdaptiveTimeScale (double s, double t) : AdaptiveTimeScale(s, t, 0.0) {}
    AdaptiveTimeScale (double s, double t, double v) : 
        geometric_scale(s), time_scale(t), too_slow_velocity(v) {}
    inline double operator () (const Particle& p) const {
        double v = p.velocity();
        if (v > too_slow_velocity) return v < (geometric_scale / time_scale) ? time_scale : (geometric_scale / v);
        else return (geometric_scale / v);
    }
    private:
    const double geometric_scale, time_scale, too_slow_velocity;
};

////////////////////////////////////////////////////////////////////////////////

template <typename F, typename Z, typename ...Cs>
class Billiard {
    public:
        inline void collision (Particle& p) const {
            base_collision (p, typename genseq<sizeof...(Cs)>::type());
        }
        inline bool is_inside (const Particle& p) const {
            return base_is_inside (p, typename genseq<sizeof...(Cs)>::type());
        }
        F fly;
    private:
        Z time_step;
        std::tuple<Cs...> domains;

        template<int ...>  struct seq {};
        template<int N, int ...S> struct genseq : genseq<N-1, N-1, S...> {};
        template<int ...S> struct genseq<0, S...>{ typedef seq<S...> type; };
        
        template<int ...S>
        inline void base_collision (Particle&, seq<S...>) const;

        template<int ...S>
        inline bool base_is_inside (const Particle&, seq<S...>) const;

};

template <typename F, typename Z, typename ...Cs> 
template <int ...S>
inline bool Billiard<F,Z,Cs...>::base_is_inside (const Particle& p, seq<S...>) const
{
    bool isInside = true;
    is_inside_aux (p, isInside, std::get<S>(domains) ...);
    return isInside;
}


template <typename F, typename Z, typename ...Cs>
template <int ...S>
inline void Billiard<F,Z,Cs...>::base_collision (Particle& p, seq<S...>) const
{
    double step = time_step (p);
    double ta, tb{0.0};

    Particle p0 = p;

    bool isCollision = false;

    while (!isCollision) {
        ta = tb;
        tb = tb + step;
        is_collision_aux (p0, ta, tb, p, isCollision, fly, std::get<S>(domains) ...);
    }
}

template<typename F>
static inline void is_collision_aux (Particle p, double ta, double tb, Particle& p1, 
                              bool& isCollision, const F& fly) {}

template<typename F, typename C, typename... Cs>
static inline void is_collision_aux (Particle p, double ta, double tb, Particle& p1, bool& isCollision, 
                              const F& fly, const C& domain, const Cs&... domains) 
{
    double tm;
    auto f = [&fly, &domain, &p] (double t, double& f, double& df) {domain.fdf (fly (p, t), f, df);};
    if (find_next_root (f, ta, tb, tm)) {
        p1 = fly (p, tm);
        domain.reflection (p1);
        isCollision = true;
        is_collision_aux (p, ta, tm, p1, isCollision, fly, domains...);
    }
    else {
        is_collision_aux (p, ta, tb, p1, isCollision, fly, domains...);
    }
}

static inline void is_inside_aux (const Particle& p, bool& isInside) {}

template<typename C, typename... Cs>
static inline void is_inside_aux (const Particle& p, bool& isInside, const C& domain, const Cs&... domains) 
{
    double f, df;
    domain.fdf (p, f, df);
    isInside = isInside && f > 0;
    is_inside_aux (p, isInside, domains...);
}

#endif
