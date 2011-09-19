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
or implied, of Frans van den Bergh.
*/
#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
using std::lower_bound;
using std::upper_bound;

double loess_core(vector<Ordered_point>& ordered, size_t start_idx, size_t end_idx,
    double mid,  Point& sol) {

    double rsq = 0;

    int n = end_idx - start_idx;
    
    if (n < 8) {
        sol.x = 0;
        sol.y = 0;
        return 1e10;
    }
    
    double span = std::max(ordered[end_idx-1].first - mid, mid - ordered[start_idx].first);
    vector<double> sig(n,1.0);
    for (int i=0; i < n; i++) {
        double d = fabs((ordered[i + start_idx].first - mid)/span) * 1.2;
        if (d > 1.0) {
            sig[i] = 2;
        } else {
            sig[i] = 1.0 / ( (1 - d*d)*(1 - d*d)*(1 - d*d) + 1);
        }
    }
    
    double sx  = 0;
    double sy  = 0;
    double ss = n;
    
    ss = 0;
    for (int i=0; i < n; i++) {
        double weight = 1.0/SQR(sig[i]);
        ss += weight;
        sx += ordered[i+start_idx].first * weight;
        sy += ordered[i+start_idx].second * weight;
    }
    double sxoss = sx / ss;
    
    double st2 = 0;
    double b = 0;
    for (int i=0; i < n; i++) {
        double t = (ordered[i+start_idx].first - sxoss) / sig[i];
        st2 += t*t;
        b += t*ordered[i+start_idx].second / sig[i];
    }
    b /= st2;
    double a = (sy - sx*b)/ss;
    sol.x = a;
    sol.y = b;
    for (int i=0; i < n; i++) {
        double r = (ordered[i+start_idx].first*sol.y + sol.x) - ordered[i+start_idx].second;
        rsq += fabs(r); // m-estimate of goodness-of-fit 
    }
    return rsq/double(n);
}

void loess_fit(vector< Ordered_point  >& ordered, double* fft_in_buffer, const int fft_size, double lower, double upper, bool deriv) {
    const int nsteps = fft_size;
    double x_span = upper - lower;
    double step = x_span / double(nsteps);
    
    size_t start_idx = 0;
    size_t end_idx = 0;
    
    int fft_idx = 0;
    for (double step_base = lower; step_base < upper; step_base += step) {
    
        double mid = step_base + 0.5*step;
        
        double min_rsq = 1e50;
        double mid_rsq = 1e50;
        double rsq = 0;
        Point sol;
        Point lsol;
        
        // try symmetric solution
        start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.25) - ordered.begin();
        end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.25) - ordered.begin();
        rsq = loess_core(ordered, start_idx, end_idx, mid, lsol);
        
        if (rsq < min_rsq) {
            min_rsq = rsq;
            mid_rsq = rsq;
            sol = lsol;
        }
        
        double w = 0.54 + 0.46*cos(2*M_PI*(fft_idx - fft_size/2)/double(fft_size-1));  // Hamming window function
        fft_idx++;
        if (deriv) {
            fft_in_buffer[fft_idx-1] = sol.y * w ; // derivative
        } else {
            fft_in_buffer[fft_idx-1] = sol.x + mid*sol.y;
        }      
    }
    if (deriv) {
        for (size_t i=0; i < 3; i++) {
            fft_in_buffer[i] = 0;
            fft_in_buffer[fft_idx-1-i] = 0;
        }
    }
}


