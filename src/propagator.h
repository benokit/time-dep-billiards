#ifndef __PROPAGATOR_H
#define __PROPAGATOR_H

#include <vector>
#include "billiard.h"

////////////////////////////////////////////////////////////////////////////////

template <typename B, typename F>
class CollisionsPropagator {
    public:
        inline void propagate(Particle&, unsigned) const;
    private:
        F time_fold;
        B billiard;
};

template <typename B, typename F>
void CollisionsPropagator<B,F>::propagate (Particle& particle, unsigned n_collisions) const
{
    while (n_collisions > 0) {
        billiard.collision (particle);
        time_fold (particle);
        n_collisions -= 1;
    }
}

template <typename B, typename F>
class TimePropagator {
    public:
        inline void propagate(Particle&, const double) const;
    private:
        F time_fold;
        B billiard;
};


template <typename B, typename F>
void TimePropagator<B,F>::propagate (Particle& particle, const double t_step) const
{
    if (t_step <= 0.0) return;

    Particle p0;
    double t = 0.0, dt = 0.0;
    while (t < t_step) {
        p0 = particle;
        billiard.collision (particle);
        dt = particle.t - p0.t;    
        t += dt;
        time_fold (particle);
    }
    dt -= t - t_step;
    particle = billiard.fly (p0, dt);
    time_fold (particle);
}

template <typename B, typename F>
class TimeTracePropagator {
    public:
        inline std::vector<Particle> propagate (Particle&, const double) const;
    private:
        F time_fold;
        B billiard;
};

template <typename B, typename F>
std::vector<Particle> TimeTracePropagator<B,F>::propagate (Particle& particle, const double t_step) const
{
    std::vector<Particle> particle_trace;
    if (t_step <= 0.0) return particle_trace;

    Particle p0;
    double t = 0.0;
    double dt = 0.0;
    do {
        particle_trace.push_back (particle);
        p0 = particle;
        billiard.collision (particle);
        dt = particle.t - p0.t;    
        t += dt;
        time_fold (particle);
    } while (t < t_step);
    dt -= t - t_step;
    particle = billiard.fly (p0, dt);
    time_fold (particle);
    
    return particle_trace;
}

template <typename Callable, typename Arg>
struct return_type_of {
    using type = decltype(std::declval<Callable>()(std::declval<Arg>()));
};

template <typename P, typename Q>
class Observer {
    public:
        using T = typename return_type_of<Q, Particle>::type;
        template <typename S>
        inline std::vector<T> sample_observable (Particle& particle, const S& steps) const {
            std::vector<T> observed_values(steps.n_steps);
            for (int i = 0; i < steps.n_steps; ++i) {
                propagator.propagate(particle, steps.step(i));
                observed_values[i] = observe(particle);
            } 
            return observed_values;
        }
    private:
        Q observe;
        P propagator;
};

////////////////////////////////////////////////////////////////////////////////

struct TimeFoldMod2Pi {
    inline void operator () (Particle& p) const {p.time_modulo (2 * M_PI);}
};

struct TimeFoldNone {
    inline void operator () (const Particle& p) const {}
};

struct TimeFoldToZero {
    inline void operator () (Particle& p) const {p.set_time (0);}
};

////////////////////////////////////////////////////////////////////////////////

struct ObserveEnergy {
    inline double operator () (const Particle& p) const { return p.vx * p.vx + p.vy * p.vy; }
};

struct ObserveVelocity {
    inline double operator () (const Particle& p) const { return sqrt (p.vx * p.vx + p.vy * p.vy); }
};

struct ObserveParticle {
    inline Particle operator () (const Particle& p) const { return p; }
};

#endif
