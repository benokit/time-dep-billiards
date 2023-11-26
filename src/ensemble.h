#ifndef __ENSEMBLE_H
#define __ENSEMBLE_H

#include<vector>

// struct RandomFrame {
//     double xmin;
//     double ymin;
//     double dx;
//     double dy;
// };

// template <typename B>
// std::vector<Particle> generate_ensemble
//     (const B& billiard, const RandomFrame& frame, double v0, double t0, int nensemble)
// {
//     std::vector<Particle> ensemble(nensemble);
//     std::default_random_engine generator;
//     std::uniform_real_distribution<double> distribution (0, 1);
//     std::for_each (ensemble.begin(), ensemble.end(), 
//         [&distribution, &generator, &billiard, &frame, t0, v0] (Particle& p)
//             { double phi, x, y;
//               do {
//                 x = frame.xmin + frame.dx * distribution (generator);
//                 y = frame.ymin + frame.dy * distribution (generator);
//                 p = (Particle) {x, y, 1, 0, t0};
//               } while (! billiard.is_inside (p));
//               phi = 2 * M_PI * distribution (generator);
//               p = (Particle) {x, y, v0 * cos (phi), v0 * sin (phi), t0};
//             }
//         );

//     return ensemble;
// }

template <typename B, typename E>
void ensemble_propagate_time (const B& b, E& e, const double t_step)
{
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < e.size(); i++) {
            map (b, e[i], t_step);
        } 
    }
}

template <typename B, typename E, typename S, typename F>
auto ensemble_sample_observable (const B& b, E& e, S& s, F& f) -> std::vector<std::vector<decltype (f(e[0]))>> 
{
    std::vector<std::vector<decltype (f(e[0]))>> data(e.size());
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < e.size(); ++i) {
            data[i] = sample_observable (b, e[i], s, f);
        } 
    }
    return data;
}

#endif