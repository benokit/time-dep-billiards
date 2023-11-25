#ifndef __FROOT_H
#define __FROOT_H

#include <cmath>

// Find next root of a function f(t) on the interval (ta, tb) in which
// the first derivative is negative: df/dt < 0.
// Type F mus support operator () (double t, double& f, double& df)
// where t is the independent variable, f is a value of the function at t
// and df is its derivative at t.
// If find_next_root returns true then f(root) = 0 and df(root) < 0.
template <typename F>
inline bool find_next_root (F fdf, double ta, double tb, double& root)
{   
    double fa, dfa, fb, dfb, fm, dfm;
    double tm, dt0, dt1, dtm;
    bool isRoot = false;

    fdf (ta, fa, dfa);
    fdf (tb, fb, dfb);

    if (fb <= 0.0) {
        isRoot = true;
    }
    else if (dfa < 0.0 && dfb > 0.0) {
        // local concaveness, check for pruning in local minimum
        // df = 0, approximately at t = tm
        tm = (dfb * ta - dfa * tb) / (dfb - dfa);
        fdf (tm, fm, dfm);
        if (fm < 0.0) {
            isRoot = true;
            tb = tm;
        }
    }

    if (isRoot) {
        // run hybrid Newton algorithm
        dt0 = tb - ta;
        tm = 0.5 * (ta + tb);
        for(;;) {
            fdf (tm, fm, dfm);
            if (fm < 0.0)
                tb = tm;
            else
                ta = tm;
            dt1 = tb - ta;
            if (dt1 < dt0 && fm != 0.0) {
                dt0 = dt1;
                // check if we are converging to the desired root
                if (dfm < 0.0) {
                    // Newton
                    dtm = -fm / dfm;
                    // interval must shrink
                    if (fabs (dtm) < dt1) 
                        // Newton 
                        tm += dtm; 
                    else 
                        // Bisection
                        tm = 0.5 * (ta + tb);
                } 
                else {
                    // Bisection
                    tm = 0.5 * (ta + tb);
                }
            }
            else {
                root = tm;
                return true;
            }
        }
    } 
    else { 
        return false;
    }
}

#endif
