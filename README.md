# Library for numerical simulations of billiard systems

- Written in C++
- Highly efficient
- Generic
- Easy to use

## Minimal example

Create a file `test.cpp`. Copy the library `src` folder in the same location as the file. Put the following content in the file:

```C++
#include "billiard.h"
#include "curve/ellipse.h"

struct EllipseCurve : public Ellipse {
        EllipseCurve () : Ellipse (2.0) {}
};

int main() {
    // define a billiard dynamical system with a static domain
    // bounded with an ellipse inside which a particle moves freely
    Billiard<FreeFlight,GenericTimeStep,Static,EllipseCurve> billiard;

    Particle particle = (Particle) {0.0, 0.0, 1.0, 1.0, 0};

    for (int i = 0; i < 10; i++) {
        billiard.collision(particle);
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
