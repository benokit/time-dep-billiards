#ifndef __STEPS_H
#define __STEPS_H

#include <array>
#include <numeric>

template <typename T, unsigned N> 
class ConstSteps {
  public: 
    ConstSteps (T st) : st(st) {}
    T step (int i) const {return st;}
    T cum_step (int i) const {return (i + 1) * st;}
    unsigned n_steps = N;
  private:
    T st;
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
    unsigned n_steps = N;
  private:
    std::array<unsigned,N> steps;
    std::array<unsigned,N> csteps;
};

#endif