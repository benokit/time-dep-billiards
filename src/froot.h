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

#ifndef __FROOT_H
#define __FROOT_H

#include <cmath>

// Find next root of a function f(t) on the interval (ta, tb) in which
// the first drerivative is negative: df/dt < 0.
// Type F mus support operator () (double t, double& f, double& df)
// where t is the independnet variable, f is a value of the function at t
// and df is its derivative at t.
// If find_next_root returns true than f(root) = 0 and df(root) < 0.
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
