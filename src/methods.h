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

#ifndef __METHODS_H
#define __METHODS_H

#include <cstdio>
#include <vector>
#include <array>
#include <numeric>
#include <functional>
#include "billiard.h"


////////////////////////////////////////////////////////////////////////////////

struct RandomFrame {
    double xmin;
    double ymin;
    double dx;
    double dy;
};

template <typename B>
std::vector<Particle> generate_ensemble
    (const B& billiard, const RandomFrame& frame, double v0, double t0, int nensemble)
{
    std::vector<Particle> ensemble(nensemble);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution (0, 1);
    std::for_each (ensemble.begin(), ensemble.end(), 
        [&distribution, &generator, &billiard, &frame, t0, v0] (Particle& p)
            { double phi, x, y;
              do {
                x = frame.xmin + frame.dx * distribution (generator);
                y = frame.ymin + frame.dy * distribution (generator);
                p = (Particle) {x, y, 1, 0, t0};
              } while (! billiard.is_inside (p));
              phi = 2 * M_PI * distribution (generator);
              p = (Particle) {x, y, v0 * cos (phi), v0 * sin (phi), t0};
            }
        );

    return ensemble;
}

////////////////////////////////////////////////////////////////////////////////


class Statistics
{
    public:

    Statistics() {};

    template<typename M>
    Statistics (const M& data) {compute (data);};

    template<typename M, typename F>
    Statistics (const M& data, const F& f) {compute (data, f);};
        
    struct Info {
        double mean;
        double var;
    };

    std::vector<Info> statvec;

    template<typename M>
    void compute (const M&);

    template<typename M, typename F>
    void compute (const M&, const F&);

    template<typename S, typename T>
    void print (const S&, T&);

    template<typename S>
    void print (const S& s) {print (s, std::cout);};
};

template <typename M>
void Statistics::compute (const M& m)
{
    int nensemble = m.size();
    int nsteps = m[0].size();
    statvec.clear();
    statvec.assign(nsteps, (Info) {0,0});
    int i, j;
    for (i = 0; i < nensemble; i++)
        for (j = 0; j < nsteps; j++) {
            statvec[j].mean += m[i][j];
            statvec[j].var += m[i][j] * m[i][j];
        }
    for (j = 0; j < nsteps; j++) {
        statvec[j].mean /= (double) nensemble;
        statvec[j].var -= nensemble * statvec[j].mean * statvec[j].mean;
        statvec[j].var /= (double) nensemble - 1;
    }
}

template <typename M, typename F>
void Statistics::compute (const M& m, const F& f)
{
    int nensemble = m.size();
    int nsteps = m[0].size();
    statvec.clear();
    statvec.assign(nsteps, (Info) {0,0});
    int i, j;
    for (i = 0; i < nensemble; i++)
        for (j = 0; j < nsteps; j++) {
            double x = f (m[i][j]);
            statvec[j].mean += x;
            statvec[j].var += x * x;
        }
    for (j = 0; j < nsteps; j++) {
        statvec[j].mean /= (double) nensemble;
        statvec[j].var -= nensemble * statvec[j].mean * statvec[j].mean;
        statvec[j].var /= (double) nensemble - 1;
    }
}

template <typename S, typename T>
void Statistics::print (const S& steps, T& file)
{
    file << std::setw(25) << "step";
    file << std::setw(25) << "mean";
    file << std::setw(25) << "var";
    file << std::endl;
    for (unsigned i = 0; i < steps.nsteps; i++)
    {
        file << std::setw(25) << std::setprecision(8) << steps.cum_step(i);
        file << std::setw(25) << std::setprecision(8) << statvec[i].mean;
        file << std::setw(25) << std::setprecision(8) << statvec[i].var;
        file << std::endl;
    } 
}

////////////////////////////////////////////////////////////////////////////////

template <typename T, unsigned N> 
class ConstSteps {
  public: 
    ConstSteps (T st) : stepconst(st) {}
    T step (int i) const {return stepconst;}
    T cum_step (int i) const {return (i + 1) * stepconst;}
    unsigned nsteps = N;
  private:
    T stepconst;
};

template <unsigned N>
class LogSteps {
  public:
    LogSteps (double b) 
    {
        const double l10 = log (10.0);
        double a, d, i = 0, f = 0;
        while (f < 1) {
          i += 1;
          a = log10 (i);
          d = (b - a) / double (N - 1);
          f = exp (l10 * (a + d)) - exp (l10 * a);
        }
        for (int j = 0; j < N; j++)
            csteps[j] = (unsigned) round (exp (l10 * (a + j * d)));
        std::adjacent_difference (csteps.begin(), csteps.end(), steps.begin());
    }
    unsigned step (int i) const {return steps[i];}
    unsigned cum_step (int i) const {return csteps[i];}
    unsigned nsteps = N;
  private:
    std::array<unsigned,N> steps;
    std::array<unsigned,N> csteps;
};

////////////////////////////////////////////////////////////////////////////////

template <typename B>
void map (const B& b, Particle& p, unsigned n)
{
    while (n > 0) {
        b.collision (p);
        b.time_fold (p);
        n -= 1;
    }
}

template <typename B>
void map (const B& b, Particle& p, const double t_step)
{
    if (t_step <= 0.0) return;

    Particle p0;
    double t = 0.0, dt = 0.0;
    while (t < t_step) {
        p0 = p;
        b.collision (p);
        dt = p.t - p0.t;    
        t += dt;
        b.time_fold (p);
    }
    dt -= t - t_step;
    p = b.fly (p0, dt);
    b.time_fold (p);
}

template <typename B, typename E>
void ensemble_map (const B& b, E& e, const double t_step)
{
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < e.size(); i++) {
            map (b, e[i], t_step);
        } 
    }
}

template <typename B>
std::vector<Particle> map_trace (const B& b, Particle& p, const double t_step)
{
    std::vector<Particle> vec;
    if (t_step <= 0.0) return vec;

    Particle p0;
    double t = 0.0;
    double dt = 0.0;
    do {
        vec.push_back (p);
        p0 = p;
        b.collision (p);
        dt = p.t - p0.t;    
        t += dt;
        b.time_fold (p);
    } while (t < t_step);
    dt -= t - t_step;
    p = b.fly (p0, dt);
    b.time_fold (p);
    
    return vec;
}

template <typename B, typename S, typename F>
auto sample_observable (const B& b, Particle& p, S& s, F& f) -> std::vector<decltype (f(p))>
{
   std::vector<decltype (f(p))> data(s.nsteps);
   for (int i = 0; i < s.nsteps; ++i) {
        map (b, p, s.step(i));
        data[i] = f(p);
   } 
   return data;
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

template <typename B, typename F>
double averaged_observable (const B& b, Particle& p, F& f, unsigned n)
{
    double x = 0;
    for (int i = 0; i < n; i++) {
        b.collision (p);
        x += f(p);
    }
    return x / (double) n;
}

template <typename B, typename S, typename F>
std::vector<double> sample_averaged_observable (const B& b, Particle& p, S& s, F& f)
{
   std::vector<double> data(s.nsteps);
   unsigned st, tr;
   for (int i = 0; i < s.nsteps; i++) {
        tr = 1 + floor (0.05 * s.step(i));
        tr = tr > 2000 ? 2000 : tr;
        st = s.step(i) - tr;
        map (b, p, st);
        data[i] = averaged_observable (b, p, f, tr);
   } 
   return data;
}

template <typename B, typename E, typename S, typename F>
std::vector<std::vector<double>> ensemble_sample_averaged_observable (const B& b, E& e, S& s, F& f)
{
    std::vector<std::vector<double>> data(e.size());
    #pragma omp parallel
    { 
        #pragma omp for schedule (runtime)
        for (int i = 0; i < e.size(); i++) {
            data[i] = sample_averaged_observable (b, e[i], s, f);
        } 
    }
    return data;
}

#endif
