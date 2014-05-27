#ifndef __ELLIPSE_H
#define __ELLIPSE_H

class Ellipse : public Curve<Ellipse> {
    public:
        Ellipse (double bb) : b(bb) {}
        inline Derivatives derivatives (const Particle& p) const {
            Derivatives d;
            d.f = 1.0 - p.x * p.x - b * p.y * p.y;
            d.dfdx = -2.0 * p.x;
            d.dfdy = -2.0 * b * p.y;
            d.dfdt = 0.0; 
            return d;
        }
    private:
        const double b;
};

#endif
