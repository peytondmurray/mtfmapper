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


void bin_fit(vector< Ordered_point  >& ordered, double* sampled, const int fft_size, double lower, double upper) {

    const double missing = -1e7;

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }

    const double scale = 2;
    for (int idx=fft_size/4-1; idx <= 3*fft_size/4+1; idx++) {
        
        double mid = idx*scale*(upper-lower)/double(fft_size-1) + scale*lower - 0.5 + 0.125;
        size_t start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.25/2.0) - ordered.begin();
        size_t end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.25/2.0) - ordered.begin();

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
            sampled[idx] = sum / double(count);
        } else {
            sampled[idx] = missing;
        }
    }

    if (sampled[fft_size/4] == missing) { // first value is missing
        int j = fft_size/4;
        while (j < fft_size && sampled[j] == missing) {
            j++;
        }
        sampled[fft_size/4] = sampled[j];
    }
    for (int idx=fft_size/4+1; idx <= 3*fft_size/4; idx++) {
        if (sampled[idx] == missing) {
            sampled[idx] = sampled[idx-1];
        }
    }

    double old = sampled[0];
    for (int idx=fft_size/4; idx <= 3*fft_size/4; idx++) {
        double w = 0.54 + 0.46*cos(scale*2*M_PI*(idx - fft_size/2)/double(fft_size-1));  // Hamming window function
        double temp = sampled[idx];
        sampled[idx] = (sampled[idx+1] - old) * w;
        old = temp;
    }

    // pad surrounding are before fft
    for (int idx=0; idx < fft_size/4+8; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=3*fft_size/4-8; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }
}
