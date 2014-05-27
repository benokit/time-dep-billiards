#ifndef __PARTICLE_H
#define __PARTICLE_H

#include <cmath>
#include <iomanip>

template<typename T>
struct ParticleT {
    double x;
    double y;
    double vx;
    double vy;
    double t;
    inline void fly (double dt) {T (this, dt);}
    inline double velocity () const {return sqrt (vx * vx + vy * vy);}
    inline void set_velocity (double v0) {double f = v0 / velocity(); vx *= f; vy *= f;}
    inline void set_time (double t0) {t = t0;} 
    inline void time_modulo (double period) {t = fmod (t, period);}
    void print () const {
        std::cout << std::setw(25) << std::setprecision(16) << x; 
        std::cout << std::setw(25) << std::setprecision(16) << y; 
        std::cout << std::setw(25) << std::setprecision(16) << vx; 
        std::cout << std::setw(25) << std::setprecision(16) << vy; 
        std::cout << std::setw(25) << std::setprecision(16) << t; 
        std::cout << std::endl;
    }
};

struct FreeFlight {
    FreeFlight (ParticleT<FreeFlight> * p, double dt) {
        p -> x += p -> vx * dt;
        p -> y += p -> vy * dt;
        p -> t += dt;
    }
};

#endif
