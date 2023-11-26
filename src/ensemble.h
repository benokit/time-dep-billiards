#ifndef __ENSEMBLE_H
#define __ENSEMBLE_H

#include <vector>
#include <random>
#include "billiard.h"

struct Frame {
    double x_min;
    double y_min;
    double dx;
    double dy;
};

template <typename B>
std::vector<Particle> generate_ensemble
    (const B& billiard, const Frame& frame, const double v0, const double t0, const int n_particles)
{
    std::vector<Particle> ensemble(n_particles);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution (0, 1);
    std::for_each (ensemble.begin(), ensemble.end(), 
        [&distribution, &generator, &billiard, &frame, t0, v0] (Particle& p)
            { double phi, x, y;
              do {
                x = frame.x_min + frame.dx * distribution (generator);
                y = frame.y_min + frame.dy * distribution (generator);
                p = (Particle) {x, y, 1, 0, t0};
              } while (! billiard.is_inside (p));
              phi = 2 * M_PI * distribution (generator);
              p = (Particle) {x, y, v0 * cos (phi), v0 * sin (phi), t0};
            }
        );

    return ensemble;
}

template <typename P, typename E>
void ensemble_propagate_time (const P& propagator, E& ensemble, const double t_step)
{
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < ensemble.size(); i++) {
            propagator.propagate(ensemble[i], t_step);
        } 
    }
}

template <typename O, typename E, typename S>
std::vector<std::vector<double>> ensemble_sample_observable (O& observer, E& ensemble, S& steps)
{
    std::vector<std::vector<double>> data(ensemble.size());
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < ensemble.size(); ++i) {
            data[i] = observer.sample_observable (ensemble[i], steps);
        } 
    }
    return data;
}

#endif