#ifndef __METHODS_H
#define __METHODS_H

#include <vector>
#include "billiard.h"

////////////////////////////////////////////////////////////////////////////////

template <typename F>
class CollisionsPropagator
{
    public:
        template <typename B>
        inline void propagate(const B&, Particle&, unsigned);
    private:
        F time_fold;
};

template <typename F>
template <typename B>
void CollisionsPropagator<F>::propagate (const B& billiard, Particle& particle, unsigned n_collisions)
{
    while (n_collisions > 0) {
        billiard.collision (particle);
        time_fold (particle);
        n_collisions -= 1;
    }
}

template <typename F>
class TimePropagator
{
    public:
        template <typename B>
        inline void propagate(const B&, Particle&, const double);
    private:
        F time_fold;
};

template <typename F>
template <typename B>
void TimePropagator<F>::propagate (const B& billiard, Particle& particle, const double t_step)
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

template <typename B>
std::vector<Particle> trace_propagate_time (const B& billiard, Particle& particle, const double t_step)
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
        billiard.time_fold (particle);
    } while (t < t_step);
    dt -= t - t_step;
    particle = billiard.fly (p0, dt);
    billiard.time_fold (particle);
    
    return particle_trace;
}

template <typename P, typename S, typename F>
auto sample_observable (P& propagate, Particle& particle, S& steps, F& observe) -> std::vector<decltype (observe(particle))>
{
   std::vector<decltype (observe(particle))> observed_values(steps.n_steps);
   for (int i = 0; i < steps.n_steps; ++i) {
        propagate(particle, steps.step(i));
        observed_values[i] = observe(particle);
   } 
   return observed_values;
}

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

#endif
