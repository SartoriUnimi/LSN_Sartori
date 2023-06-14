#ifndef __FUNCS__
#define __FUNCS__

#include <cmath>
#include "random.h"
#include <cstdlib>
#include <cfloat>

using namespace std;


double error(double AV, double AV2, int n){
    if(n==0)
        return 0;
    double error = sqrt(abs(AV2-AV*AV)/n);
    return error;
};

int lancio(double L, double d, Random& rnd){
    // l'ago lanciato cade ad una certa y in (0,1)
    double y0 = rnd.Rannyu();
    // e l'altra estremità disterà in verticale L*sin(x) in (-L; L)
    double y1 = y0 + sin(rnd.Rannyu()*DBL_MAX)*L;
    // abbiamo colpito se le due estremità cadono in due strisce differenti
    if(floor(y0/d)!=floor(y1/d)){
        return 1;
    } else {
        return 0;
    }
}


#endif // __FUNCS__