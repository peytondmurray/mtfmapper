/*
Copyright 2011 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/

#include "bessel.h"

#include <cmath>

// Rather simple implementation of the Bessel J function of order 1
// it is accurate within about 1e-6, which is good enough for what
// we use it for here

static double p1r(double x, double k = 0) {
    const double n = 1.0;
    double num = (4*n*n - (4*k-3)*(4*k-3))*(4*n*n - (4*k-1)*(4*k-1));
    double denom = 2*k*(2*k-1)*(8*x)*(8*x);
    
    if (k == 10) {
        return 1 - num/denom;
    }
    
    if (k == 0) {
        return 1 + (3*5)/(2*(8*x*8*x))*p1r(x, k+1);
    }
    
    return 1 - (num/denom)*p1r(x, k+1);
}

static double q1r(double x, double k = 0) {
    const double n = 1.0;
    double num = (4*n*n - (4*k-1)*(4*k-1))*(4*n*n - (4*k+1)*(4*k+1));
    double denom = 2*k*(2*k+1)*(8*x)*(8*x);
    
    if (k == 10) {
        return 1 - num/denom;
    }
    
    if (k == 0) {
        return 3*q1r(x, k+1)/(8*x);
    }
    
    return 1 - num/denom*q1r(x, k+1);
}

static double asymptotic(double x) {
    double sign = 1.0;
    if (x < 0) {
        x = fabs(x);
        sign = -1.0;
    }
    double f =  p1r(x)*cos(x - 3*M_PI/4.0);
    f -= q1r(x)*sin(x - 3*M_PI/4.0);
    return sign * f * sqrt(2.0/(M_PI*x));
}

double besselJ(double x, double k) {
    if (fabs(x) > 24) return asymptotic(x); // when |x| exceeds 24, the error in the asymptotic approximation is lower
    double y = (x*0.5)*(x*0.5);
    double denom = k*(k+1);
    
    if (k == 60) {
        return 1 - y/denom;
    }
    
    if (k == 0) {
        return 0.5*x*besselJ(x, k+1);
    }
    
    return 1 - y/denom*besselJ(x, k+1);
}

