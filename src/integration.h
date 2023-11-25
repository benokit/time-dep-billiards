#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include <cmath>

namespace GaussKronrod 
{

        struct Rule {
            double node; double kweight; double gweight;
        };

        const Rule rule[8] = 
           {{0.9914553711208126392070, 0.0229353220105292249637, 0},	
            {0.9491079123427585245260, 0.0630920926299785532910, 0.1294849661688696932710},
            {0.8648644233597690727900, 0.1047900103222501838400, 0},
            {0.7415311855993944398639, 0.1406532597155259187450, 0.2797053914892766679010},
            {0.5860872354676911302941, 0.1690047266392679028266, 0},	
            {0.4058451513773971669066, 0.1903505780647854099133, 0.3818300505051189449500},
            {0.2077849550078984676007, 0.2044329400752988924140, 0},	
            {0.0000000000000000000000, 0.2094821410847278280130, 0.4179591836734693877551}};

        template <typename F> 
        double gauss_kronrod (const F f, double a, double b, double& err)
        {
            double xa, ya, xb, yb, ik = 0.0, ig = 0.0;
            double q1 = 0.5 * (b + a);
            double q2 = 0.5 * (b - a);
            int i;
            for (i = 0; i < 7; i++) {
                xa = q1 + q2 * rule[i].node;
                ya = f (xa);
                xb = q1 - q2 * rule[i].node;
                yb = f (xb);
                ik += rule[i].kweight * (ya + yb);
                ig += rule[i].gweight * (ya + yb);
            }
            xa = q1 + q2 * rule[7].node;
            ya = f (xa);
            ik += rule[7].kweight * ya;
            ig += rule[7].gweight * ya;
            err = exp (1.5 * log (200 * q2 * fabs (ik - ig)));
            return q2 * ik;
        }

        template <typename F>
        double integrate (const F f, double xa, double xb)
        {   
            double err;
            double integ = gauss_kronrod (f, xa, xb, err);
            double xm;
            if (err > 1e-15) {
                xm = 0.5 * (xa + xb);    
                integ = integrate (f, xa, xm) + integrate (f, xm, xb);
            }
            return integ;
        }
}


#endif
