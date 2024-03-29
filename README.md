# Library for numerical simulations of billiard dynamical systems

- Written in C++
- Highly efficient
- Suitable for massive numerical computations
- Works for arbitrary time-dependent billiard domains
- Easy to use

## Minimal example

In this example it is shown how to compute 10 consecutive collisions of a particle in a static elliptical billiard whose domain is an ellipse with `a = 1` and `b = 2`. Between the collisions the particle is moving freely - with constant velocity.

Create a file `test.cpp`. Copy the library `src` folder in the same location as the file. Put the following content in the file:

```c++
#include "billiard.h"
#include "domains/ellipse.h"

struct EllipseDomain : public Ellipse {
        EllipseDomain () : Ellipse (2.0) {}
};

// Billiard is defined with a set of template parameters
Billiard<FreeFlight,ConstantTimeScale,EllipseDomain> billiard;

int main() {
    // state of a particle is represented with a vector (x, y, vx, vy, t)
    Particle particle = (Particle) {0.0, 0.0, 1.0, 1.0, 0};

    for (int i = 0; i < 10; i++) {
        // billiard.collision() moves the particle state immediately after the next collision
        billiard.collision(particle);

        // print the particle state: x y vx vy t
        particle.print();
    }
}
```

compile with `g++ -std=c++20 test.cpp -I ./src -o test` (or other compiler).
Running `./test` prints states of the particle after 10 consecutive collisions:

```text
       0.5773502691896257       0.5773502691896257                     -0.2                     -1.4       0.5773502691896257
       0.4023956421624665      -0.6473321200004893       -0.958368026644463        1.039966694421315        1.452123404325422
      -0.6749534173542974        0.521746051447328        0.555520778346515       -1.300537067839776        2.576273031153948
      -0.1536384246173061      -0.6987113976747168       0.8248085734082959        1.148777966898814        3.514698762570627
       0.7070503770300063       0.5000398805804025      -0.8079675037546521       -1.160684501867961         4.55820004201387
      -0.1291384196177312      -0.7011858771319608      -0.5824124391510305        1.288718646843503        5.593128789472079
      -0.6802579052490763       0.5182900647061155       0.9505091817192469       -1.047154379959043        6.539398904229071
       0.3830800733073194      -0.6531652384484568        0.234454388266455        1.394643732220743        7.658102435985668
       0.5889526751409472       0.5714607363784157      -0.9996808917424155       -1.000319006459983        8.536194773501039
      -0.5651339088957582      -0.5833625223718957       0.1651731605513426        1.404534736855405        9.690649753654922
```

## Define billiard domain

A billiard domain is defined as a domain of a continuous function `f(x, y, t)` where `f > 0`. A boundary of a billiard is a curve `f = 0`.

A boundary is defined as a type derived from a generic `Domain` type and implements derivatives function which provides a value of `f` and its first derivatives at any given point of space and time.

Here is an equivalent of the above example but with a static billiard domain defined with `f = 1 - x^4 - y^2`:

```c++
#include "billiard.h"
#include "domain.h"

struct FlattenedCircle : public Domain<FlattenedCircle> {
    inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            double x2 = p.x * p.x; 
            d.f = 1.0 - x2 * x2 - p.y * p.y; // f = 1 - x^4 - y^3
            d.dfdx = 4.0 * x2 * p.x;
            d.dfdy = 2.0 * p.y;
            d.dfdt = 0.0; 
            return d;
        }
};

// Billiard domain is specified with the third template parameter
Billiard<FreeFlight,ConstantTimeScale,FlattenedCircle> billiard;

int main() {
    Particle particle = (Particle) {0.0, 0.0, 1.0, 1.0, 0};

    for (int i = 0; i < 10; i++) {
        billiard.collision(particle);
        particle.print();
    }
}
```

## Define complex billiard domains with composition

A billiard domain can be defined as a an intersection of positive domains of arbitrarily many functions `f1(x, y, t), f2(x, y, t), ...`.

For example, we can cut the billiard domain in the above example with a half plane:

```c++
#include "billiard.h"
#include "domain.h"

struct FlattenedCircle : public Domain<FlattenedCircle> {
    inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            double x2 = p.x * p.x; 
            d.f = 1.0 - x2 * x2 - p.y * p.y;
            d.dfdx = 4.0 * x2 * p.x;
            d.dfdy = 2.0 * p.y;
            d.dfdt = 0.0; 
            return d;
        }
};

struct MyHalfPlane : public Domain<MyHalfPlane> {
    inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            d.f = p.x + p.y + 0.1;
            d.dfdx = 1.0;
            d.dfdy = 1.0;
            d.dfdt = 0.0; 
            return d;
        }
};

// domains are simply added to billiard's template parameters from position three onward
Billiard<FreeFlight,ConstantTimeScale,FlattenedCircle,MyHalfPlane> billiard;

int main() {
    Particle particle = (Particle) {0.0, 0.0, -1.0, -1.0, 0};

    for (int i = 0; i < 10; i++) {
        billiard.collision(particle);
        particle.print();
    }
}
```

## Transform billiard domains

Billiard domains can be easily arbitrarily transformed. Some basic transformations (rotation, translation, scaling, ...) are already defined.

Here is an example how easy it is to define a rotating billiard:

```c++
#include "billiard.h"
#include "domain.h"
#include "transform.h"

struct FlattenedCircle : public Domain<FlattenedCircle> {
    inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            double x2 = p.x * p.x; 
            d.f = 1.0 - x2 * x2 - p.y * p.y;
            d.dfdx = 4.0 * x2 * p.x;
            d.dfdy = 2.0 * p.y;
            d.dfdt = 0.0; 
            return d;
        }
};

// Driver defines an angle q as a function of time an its time derivative dq 
struct Driver {
    Drive operator() (double t) const {
        Drive d;
        d.q = t;
        d.dq = 1.0;
        return d;
    }
};

// Apply rotation to the domain thus yielding a new domain 
using RotatingFlattenCircle = TransformDomain<Rotation<Driver>,FlattenedCircle>;

// characteristic time scale
struct TimeScale : public AdaptiveTimeScale {
    TimeScale() : AdaptiveTimeScale(1.0, 0.1) {}
};

Billiard<FreeFlight,TimeScale,RotatingFlattenCircle> billiard;

int main() {
    Particle particle = (Particle) {0.0, 0.0, -1.0, -1.0, 0};

    for (int i = 0; i < 10; i++) {
        billiard.collision(particle);
        particle.print();
    }
}
```
