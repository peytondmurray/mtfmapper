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
    double ss  = 0;
    
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


double bin_fit(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, double lower, double upper, vector<double>& esf) {

    const double missing = -1e7;

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }

    const double scale = 2;
    double rightsum = 0;
    int rightcount = 0;
    double leftsum = 0;
    int leftcount = 0;
    int left_non_missing  = 0;
    int right_non_missing = 0;
    for (int idx=fft_size/4-1; idx <= 3*fft_size/4+1; idx++) {
        
        double mid = idx*scale*(upper-lower)/double(fft_size-1) + scale*lower;
        size_t start_idx = lower_bound(ordered.begin(), ordered.end(), mid - 0.5) - ordered.begin();
        size_t end_idx   = lower_bound(ordered.begin(), ordered.end(), mid + 0.5) - ordered.begin();
        const double left  = mid - 0.125/2.0;
        const double right = mid + 0.125/2.0;

        int included = 0;

        if (end_idx - start_idx > 2) {
            const double lpwidth = 0.33333333;
            double weight = 0;
            double sum = 0;
            for (size_t k=start_idx; k < end_idx; k++) {
                double lb = ordered[k].first - lpwidth/2;
                double rb = ordered[k].first + lpwidth/2;
                
                
                
                if (lb < right && rb > left) { // we have non-zero intersection
                    if (lb < left) lb = left;
                    if (rb > right) rb = right;
                    double w = rb - lb;
                    assert(w >= 0);
                    sum += ordered[k].second * w;
                    weight += w;
                    included++;
                }
            }
            if (weight > 0) {
                sampled[idx] = sum / weight;
                if (idx < fft_size/2 - fft_size/8) {
                    leftsum += sampled[idx];
                    leftcount++;
                }
                if (idx > fft_size/2 + fft_size/8) {
                    rightsum += sampled[idx];
                    rightcount++;
                }
                if (!left_non_missing) {
                    left_non_missing = idx; // first non-missing value from left
                }
                right_non_missing = idx; // last non-missing value
            } else {
                sampled[idx] = missing;
            }
        } else {
            sampled[idx] = missing;
        }
        
    }
    
    // now just pad out the ends of the sequences with the last non-missing values
    for (int idx=left_non_missing-1; idx >= fft_size/4-16; idx--) {
        sampled[idx] = sampled[left_non_missing];
    }
    for (int idx=right_non_missing+1; idx < 3*fft_size/4 + 16; idx++) {
        sampled[idx] = sampled[right_non_missing];
    }
    
    // Reverse the ESF if necessary
    leftsum /= double(leftcount);
    rightsum /= double(rightcount);
    if (leftsum < rightsum) {
        for (int idx=0; idx <= fft_size/4 + 15; idx++) {
            double t = sampled[fft_size/2 - idx - 1];
            sampled[fft_size/2 - idx - 1] = sampled[fft_size/2 + idx];
            sampled[fft_size/2 + idx] = t;
        }
    }
    
    int lidx = 0;
    for (int idx=fft_size/4; idx < 3*fft_size/4; idx++) {
        esf[lidx++] = sampled[idx+3];
    }
    
    vector<double> med;
    for (int idx=fft_size/4+1; idx < fft_size/2-32; idx++) {
        med.push_back(fabs(sampled[idx] - sampled[idx-1]));
    }
    for (int idx=fft_size/2+32; idx < 3*fft_size/4; idx++) {
        med.push_back(fabs(sampled[idx] - sampled[idx-1]));
    }
    sort(med.begin(), med.end());
    vector<double> lmed;
    for (int idx=fft_size/4+1; idx < fft_size/4+1+3*8; idx++) {
        lmed.push_back(sampled[idx]);
    }
    sort(lmed.begin(), lmed.end());
    double left_median = lmed[lmed.size()/2];
    vector<double> rmed;
    for (int idx=3*fft_size/4 - (1+3*8); idx < 3*fft_size/4; idx++) {
        rmed.push_back(sampled[idx]);
    }
    sort(rmed.begin(), rmed.end());
    double right_median = rmed[rmed.size()/2];

    double noise_est = std::max(1.0,med[9*med.size()/10]);

    double old = sampled[0];
    const double alpha = 0.6;
    const int tukey_w = fft_size/2; 
    // TODO: the window function really should be statically computed ...
    for (int idx=fft_size/4; idx <= 3*fft_size/4; idx++) {
        double lx = idx - fft_size/4;
        double w = 1.0;
        if (lx < alpha*(tukey_w-1)/2.0) {
            w = 0.5 + 0.5*cos(M_PI*(2*lx/(alpha * (tukey_w-1)) - 1));
        }
        if (lx > (tukey_w - 1)*(1 - alpha/2)) {
            w = 0.5 + 0.5*cos(M_PI*(2*lx/(alpha * (tukey_w-1)) - 2/alpha + 1));
        }
        double temp = sampled[idx];
        sampled[idx] = (sampled[idx+1] - old) * w;
        old = temp;
    }
    
    // pad surrounding area before fft
    for (int idx=0; idx < fft_size/4+1; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=3*fft_size/4-1; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }

    return fabs(right_median - left_median)/noise_est;
}
