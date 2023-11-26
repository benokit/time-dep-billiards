#ifndef __TRANSFORM_H
#define __TRANSFORM_H

#include "billiard.h"

struct Jacobian {
    double xp;
    double yp;
    double dxpdx;
    double dxpdy;
    double dxpdt;
    double dypdx;
    double dypdy;
    double dypdt;
};

////////////////////////////////////////////////////////////////////////////////

template <typename T>
class Transform {
    public:
        inline Particle operator () (const Particle& p) const;
        inline Particle inverse (const Particle& p) const;
    private:
        inline Particle transform (const Jacobian&, const Particle&) const;
};

template <typename T>
inline Particle Transform<T>::transform (const Jacobian& j, const Particle& p) const
{
    Particle pt;
    pt.x = j.xp;
    pt.y = j.yp;
    pt.vx = j.dxpdx * p.vx + j.dxpdy * p.vy + j.dxpdt;
    pt.vy = j.dypdx * p.vx + j.dypdy * p.vy + j.dypdt;
    pt.t = p.t;
    return pt;
}

template <typename T>
inline Particle Transform<T>::operator () (const Particle& p) const
{
    return transform(static_cast<const T*>(this) -> jacobian (p), p);
}

template <typename T>
inline Particle Transform<T>::inverse (const Particle& p) const
{
    return transform(static_cast<const T*>(this) -> inverse_jacobian (p), p);
}

////////////////////////////////////////////////////////////////////////////////

template <typename T, typename C>
class TransformDomain : public Domain<TransformDomain<T,C>> {
    public:
        inline Derivatives derivatives (const Particle&) const;
    private:
         C domain;
         T transform;
};

template <typename T, typename C>
inline Derivatives TransformDomain<T,C>::derivatives (const Particle& p) const
{
    Particle pt = transform (p);
    Jacobian j  = transform.jacobian (p);
    Derivatives da = domain.derivatives (pt);
    Derivatives db;
    db.f = da.f;
    db.dfdx = da.dfdx * j.dxpdx + da.dfdy * j.dypdx;
    db.dfdy = da.dfdx * j.dxpdy + da.dfdy * j.dypdy;
    db.dfdt = da.dfdx * j.dxpdt + da.dfdy * j.dypdt + da.dfdt;
    return db;
}

////////////////////////////////////////////////////////////////////////////////

struct Drive {
    double q;
    double dq;
};

struct Drive2 {
    double c;
    double dc;
    double s;
    double ds;
};

////////////////////////////////////////////////////////////////////////////////

template <typename Q>
class Translation : public Transform<Translation<Q>> {
    public:
        inline Jacobian jacobian (const Particle&) const;
        inline Jacobian inverse_jacobian (const Particle&) const;
    private:
        Q driver;
};

template <typename Q>
inline Jacobian Translation<Q>::jacobian (const Particle& p) const
{
    Drive2 d = driver (p.t);
    Jacobian j;
    j.xp = p.x - d.c;
    j.yp = p.y - d.s;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = -d.dc;
    j.dypdx = 0.0;
    j.dypdy = 1.0;
    j.dypdt = -d.ds;
    return j;
}

template <typename Q>
inline Jacobian Translation<Q>::inverse_jacobian (const Particle& p) const
{
    Drive2 d = driver (p.t);
    Jacobian j;
    j.xp = p.x + d.c;
    j.yp = p.y + d.s;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = d.dc;
    j.dypdx = 0.0;
    j.dypdy = 1.0;
    j.dypdt = d.ds;
    return j;
}

////////////////////////////////////////////////////////////////////////////////

template <typename Q>
class Rotation : public Transform<Rotation<Q>> {
    public:
        inline Jacobian jacobian (const Particle&) const;
        inline Jacobian inverse_jacobian (const Particle&) const;
    private:
        Q driver;
};

template <typename Q>
inline Jacobian Rotation<Q>::jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double c = cos (d.q);
    double s = sin (d.q);
    double dc = -s * d.dq;
    double ds = c * d.dq;
    j.xp =  c * p.x + s * p.y;
    j.yp = -s * p.x + c * p.y;
    j.dxpdx = c;
    j.dxpdy = s;
    j.dxpdt = dc * p.x + ds * p.y;
    j.dypdx = -s;
    j.dypdy = c;
    j.dypdt = -ds * p.x + dc * p.y;
    return j;
}

template <typename Q>
inline Jacobian Rotation<Q>::inverse_jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double c = cos (d.q);
    double s = sin (d.q);
    double dc = -s * d.dq;
    double ds = c * d.dq;
    j.xp =  c * p.x - s * p.y;
    j.yp =  s * p.x + c * p.y;
    j.dxpdx = c;
    j.dxpdy = -s;
    j.dxpdt = dc * p.x - ds * p.y;
    j.dypdx = s;
    j.dypdy = c;
    j.dypdt = ds * p.x + dc * p.y;
    return j;
}

////////////////////////////////////////////////////////////////////////////////

template <typename Q>
class Scaling : public Transform<Scaling<Q>> {
    public:
        inline Jacobian jacobian (const Particle&) const;
        inline Jacobian inverse_jacobian (const Particle&) const;
    private:
        Q driver;
};

template <typename Q>
inline Jacobian Scaling<Q>::jacobian (const Particle& p) const
{
    Drive2 d = driver (p.t);
    Jacobian j;
    j.xp = d.c * p.x;
    j.yp = d.s * p.y;
    j.dxpdx = d.c;
    j.dxpdy = 0.0;
    j.dxpdt = d.dc * p.x;
    j.dypdx = 0.0;
    j.dypdy = d.s;
    j.dypdt = d.ds * p.y;
    return j;
}

template <typename Q>
inline Jacobian Scaling<Q>::inverse_jacobian (const Particle& p) const
{
    Drive2 d = driver (p.t);
    Jacobian j;
    j.xp = p.x / d.c;
    j.yp = p.y / d.s;
    j.dxpdx = 1 / d.dc;
    j.dxpdy = 0.0;
    j.dxpdt = - d.dc * p.x / (d.c * d.c);
    j.dypdx = 0.0;
    j.dypdy = 1 / d.ds;
    j.dypdt = - d.ds * p.y / (d.s * d.s);
    return j;
}

////////////////////////////////////////////////////////////////////////////////

template <typename Q>
class Deform : public Transform<Deform<Q>> {
    public:
        inline Jacobian jacobian (const Particle&) const;
        inline Jacobian inverse_jacobian (const Particle&) const;
    private:
        Q driver;
};

template <typename Q>
inline Jacobian Deform<Q>::jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double w = 1.0 + d.q * (p.x * p.x - 1.0);
    double dwdx = 2.0 * d.q * p.x;
    double dwdt = d.dq * (p.x * p.x - 1.0);
    j.xp = p.x;
    j.yp = p.y / w;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = 0.0;
    j.dypdx = - dwdx * (p.y / (w * w));
    j.dypdy = 1.0 / w;
    j.dypdt = - dwdt * (p.y / (w * w));
    return j;
}

template <typename Q>
inline Jacobian Deform<Q>::inverse_jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double w = 1.0 + d.q * (p.x * p.x - 1.0);
    double dwdx = 2.0 * d.q * p.x;
    double dwdt = d.dq * (p.x * p.x - 1.0);
    j.xp = p.x;
    j.yp = p.y * w;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = 0.0;
    j.dypdx = dwdx * p.y;
    j.dypdy = w;
    j.dypdt = dwdt * p.y;
    return j;
}

////////////////////////////////////////////////////////////////////////////////

template <typename Q>
class Swing : public Transform<Swing<Q>> {
    public:
        inline Jacobian jacobian (const Particle&) const;
        inline Jacobian inverse_jacobian (const Particle&) const;
    private:
        Q driver;
};

template <typename Q>
inline Jacobian Swing<Q>::jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double w = 1.0 + d.q * p.x;
    double dwdx = d.q;
    double dwdt = d.dq * p.x;
    j.xp = p.x;
    j.yp = p.y / w;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = 0.0;
    j.dypdx = - dwdx * (p.y / (w * w));
    j.dypdy = 1.0 / w;
    j.dypdt = - dwdt * (p.y / (w * w));
    return j;
}

template <typename Q>
inline Jacobian Swing<Q>::inverse_jacobian (const Particle& p) const
{
    Drive d = driver (p.t);
    Jacobian j;
    double w = 1.0 + d.q * p.x;
    double dwdx = d.q;
    double dwdt = d.dq * p.x;
    j.xp = p.x;
    j.yp = p.y * w;
    j.dxpdx = 1.0;
    j.dxpdy = 0.0;
    j.dxpdt = 0.0;
    j.dypdx = dwdx * p.y;
    j.dypdy = w;
    j.dypdt = dwdt * p.y;
    return j;
}

#endif
