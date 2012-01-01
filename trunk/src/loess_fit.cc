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
#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
using std::lower_bound;
using std::upper_bound;

const int MIN_POINTS_TO_FIT = 8;

double loess_core(vector<Ordered_point>& ordered, size_t start_idx, size_t end_idx,
    double mid,  Point& sol) {

    double rsq = 0;

    int n = end_idx - start_idx;
    
    if (n < MIN_POINTS_TO_FIT) {
        sol.x = 0;
        sol.y = 0;
        return 1e10;
    }
    
    double span = std::max(ordered[end_idx-1].first - mid, mid - ordered[start_idx].first);
    vector<double> sig(n,1.0);
    for (int i=0; i < n; i++) {
        double d = fabs((ordered[i + start_idx].first - mid)/span) / 1.2;
        if (d > 1.0) {
            sig[i] = 20;
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
    fft_in_buffer[0] = 0;
    for (double step_base = lower; step_base < upper; step_base += step) {
    
        double mid = step_base + 0.5*step;
        
        double rsq = 0;
        Point sol;
        
        // try symmetric solution
        start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.35) - ordered.begin();
        end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.35) - ordered.begin();

        bool end_capped = false;
        bool start_capped = false;

        // if we have too few points, expand the range a bit
        while ( (end_idx - start_idx) < MIN_POINTS_TO_FIT && !(end_capped && start_capped)) {
            if (end_idx < ordered.size()-1) {
                end_idx += 1;
            } else {
                end_capped = true;
            }
            if (start_idx > 0) {
                start_idx -= 1;
            } else {
                start_capped = true;
            }
        }
        // it is possible that we still have fewer than MIN_POINTS_TO_FIT, 
        // so this would be a good place to raise a warning
        
        rsq = loess_core(ordered, start_idx, end_idx, mid, sol);
        bool missing = false;
        if (fabs(sol.y) > 65535) {
            sol.y = sol.y < 0 ? -65535 : 65535;
            missing = true;
        }
        if (sol.y*mid + sol.x < 0 || sol.y*mid + sol.x > 65535) {
            sol.x = -sol.y*mid;
            missing = true;
        }
              
        double w = 0.54 + 0.46*cos(2*M_PI*(fft_idx - fft_size/2)/double(fft_size-1));  // Hamming window function
        fft_idx++;
        if (deriv) {
            fft_in_buffer[fft_idx-1] = sol.y * w ; // derivative
            //fft_in_buffer[fft_idx-1] = sol.x + mid*sol.y;
        } else {
            fft_in_buffer[fft_idx-1] = sol.x + mid*sol.y;
        }
        if (missing && fft_idx > 1) {
            fft_in_buffer[fft_idx-1] = fft_in_buffer[fft_idx-2];
        }
        //if (!missing) fprintf(stderr, "%lf\n", sol.x + mid*sol.y);
        //fprintf(stderr, "%lf\n", fft_in_buffer[fft_idx-1]);
    }
    if (deriv) {
        //for (size_t i=fft_idx-1; i > 0; i--) {
        //    fft_in_buffer[i] -= fft_in_buffer[i-1];
        //}
        for (size_t i=0; i < 5; i++) {
            fft_in_buffer[i] = 0;
            fft_in_buffer[fft_idx-1-i] = 0;
        }
        double old_x = fft_in_buffer[84];
        const double alpha = 2.0/(80);
        for (size_t i=83; i > 0; i--) {
            old_x = old_x * (1-alpha) + fft_in_buffer[i] * alpha;
            if (i < 64) {
                fft_in_buffer[i] = old_x;
            }
        }
        for (size_t i=0; i < 64; i++) {
            old_x = old_x * (1-alpha) + fft_in_buffer[fft_idx - 1 - i] * alpha;
            fft_in_buffer[fft_idx - 1 - i] = old_x;
        }
    }
    //for (size_t i=0; i < fft_idx; i++) {
    //    fprintf(stderr, "%lf\n", fft_in_buffer[i]);
    //}
}

void bin_fit(vector< Ordered_point  >& ordered, double* fft_in_buffer, const int fft_size, double lower, double upper, bool deriv) {

    const double missing = -1e7;

    size_t start_idx = 0;
    size_t end_idx = 0;

    for (int fft_idx=0; fft_idx < fft_size; fft_idx++) {
        const double scale = 2;
        double mid = fft_idx*scale*(upper-lower)/double(fft_size) + scale*lower;
        start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.25/2.0) - ordered.begin();
        end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.25/2.0) - ordered.begin();

        if (end_idx - start_idx > 2) {
            vector<double> vals;
            for (size_t k=start_idx; k < end_idx; k++) {
                vals.push_back(ordered[k].second);
            }
            sort(vals.begin(), vals.end());
            double sum = 0;
            int count = 0;
            if (vals.size() > 5) {
                for (int j=int(vals.size()*0.1);  j < int(vals.size()*0.9); j++) {
                    sum += vals[j];
                    count++;
                }
            } else {
                for (int j=0;  j < int(vals.size()); j++) {
                    sum += vals[j];
                    count++;
                }
            }
            fft_in_buffer[fft_idx] = sum / double(count);
        } else {
            fft_in_buffer[fft_idx] = missing;
        }
    }
    if (fft_in_buffer[0] == missing) { // first value is missing
        int j = 0;
        while (j < fft_size && fft_in_buffer[j] == missing) {
            j++;
        }
        fft_in_buffer[0] = fft_in_buffer[j];
    }
    for (int fft_idx=1; fft_idx < fft_size; fft_idx++) {
        if (fft_in_buffer[fft_idx] == missing) {
            fft_in_buffer[fft_idx] = fft_in_buffer[fft_idx - 1];
        }
        //fprintf(stderr, "%lf\n", fft_in_buffer[fft_idx]);
    }
    double old = fft_in_buffer[0];
    for (int fft_idx=1; fft_idx < fft_size-1; fft_idx++) {
        double w = 0.54 + 0.46*cos(2*2*M_PI*(fft_idx - fft_size/2)/double(fft_size-1));  // Hamming window function
        fft_in_buffer[fft_idx] = (fft_in_buffer[fft_idx+1] - fft_in_buffer[fft_idx]) * w * 0.5;
        double temp = fft_in_buffer[fft_idx];
        //fft_in_buffer[fft_idx] = (fft_in_buffer[fft_idx+1] - old) * w;
        old = temp;
    }

    // suppress the endpoints
    for (size_t i=0; i < 5; i++) {
        fft_in_buffer[i] = 0;
        fft_in_buffer[fft_size-1-i] = 0;
    }
    #if 0
    // now perform additional smoothing around ends
    const int sstart = 83;
    double old_x = fft_in_buffer[sstart+1];
    const double alpha = 2.0/(80);
    for (size_t i=sstart; i > 0; i--) {
        old_x = old_x * (1-alpha) + fft_in_buffer[i] * alpha;
        if (i < sstart-20) {
            fft_in_buffer[i] = old_x;
        }
    }
    for (size_t i=0; i < sstart-20; i++) {
        old_x = old_x * (1-alpha) + fft_in_buffer[fft_size - 1 - i] * alpha;
        fft_in_buffer[fft_size - 1 - i] = old_x;
    }
    #endif

    for (int fft_idx=1; fft_idx < fft_size-1; fft_idx++) {
        //fprintf(stderr, "%lf\n", fft_in_buffer[fft_idx]);
    }
}
