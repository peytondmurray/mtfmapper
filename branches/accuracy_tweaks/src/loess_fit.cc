#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
using std::lower_bound;
using std::upper_bound;

double loess_core(vector<Ordered_point>& ordered, size_t start_idx, size_t end_idx,
    double mid,  Point& sol, int mode) {

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
            //sig[i] = 1.0 / ( (1 - d*d*d)*(1 - d*d*d)*(1 - d*d*d) + 1);
            sig[i] = 1.0 / ( (1 - d*d)*(1 - d*d)*(1 - d*d) + 1);
            //sig[i] = 0.05 / ( (1 - d*d)*(1 - d*d)*(1 - d*d) + 0.05);
        }
    }
    for (int i=0; i < n; i++) {
        if (mode == 0) {
            if (i < n/2) {
                sig[i] *= 1.0;
            } else {
                sig[i] *= 0.05;
            }
        } else if (mode == 2) {
            if (i < n/2) {
                sig[i] *= 0.05;
            } else {
                sig[i] *= 1.0;
            }
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
        rsq = loess_core(ordered, start_idx, end_idx, mid, lsol, 1);
        
        if (rsq < min_rsq) {
            min_rsq = rsq;
            mid_rsq = rsq;
            sol = lsol;
        }
        
        #if 0
        // try right-only solution
        start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.35) - ordered.begin();
        end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.15) - ordered.begin();
        rsq = loess_core(ordered, start_idx, end_idx, mid, lsol, 0);
        
        if (rsq < min_rsq*0.75) {
            min_rsq = rsq;
            sol = lsol;
        }
        
        // try left-only solution
        start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.15) - ordered.begin();
        end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.35) - ordered.begin();
        rsq = loess_core(ordered, start_idx, end_idx, mid, lsol, 2);
            
        if (rsq < mid_rsq*0.75 && rsq < min_rsq) {
          min_rsq = rsq;
          sol = lsol;
        }
        #endif
        
        double w = 0.54 + 0.46*cos(2*M_PI*(fft_idx - fft_size/2)/double(fft_size-1));  // Hamming window function
        fft_idx++;
        if (deriv) {
            fft_in_buffer[fft_idx-1] = sol.y * w ; // derivative
        } else {
            fft_in_buffer[fft_idx-1] = sol.x + mid*sol.y;
        }
        fprintf(stderr, "%lf\n", sol.x + mid*sol.y);      
    }
    if (deriv) {
        for (size_t i=0; i < 3; i++) {
            fft_in_buffer[i] = 0;
            fft_in_buffer[fft_idx-1-i] = 0;
        }
    }
}


