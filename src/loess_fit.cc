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
#include "include/logger.h"
#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include "include/mtf_tables.h"
#include "include/fastexp.h"
#include "include/loess_parms.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <array>
#include <thread>
using std::lower_bound;
using std::upper_bound;

const int MIN_POINTS_TO_FIT = 8;

double loess_core(vector<Ordered_point>& ordered, size_t start_idx, size_t end_idx,
    double mid,  Point2d& sol) {

    double rsq = 0;

    int n = end_idx - start_idx;
    
    if (n < MIN_POINTS_TO_FIT) {
        sol.x = 0;
        sol.y = 0;
        return 1e10;
    }
    
    double span = max(ordered[end_idx-1].first - mid, mid - ordered[start_idx].first);
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

inline double sgn(double x) {
    return x < 0 ? -1 : 1;
}

inline double tri(double x) {
    const double alpha = 0.95;
    return alpha*(1 - fabs(x*8))*0.5*(sgn(x+0.125) - sgn(x-0.125)) + (1-alpha)*(1 - fabs(x*1.5))*0.5*(sgn(x+0.6666666) - sgn(x-0.6666666));
}

int bin_fit(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, double lower, double upper, vector<double>& esf, bool allow_peak_shift) {
    
    vector<double> weights(fft_size, 0);
    vector<double> mean(fft_size, 0);

    constexpr double missing = -1e7;

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }
    
    int rval = 0;
    
    const int fft_size2 = fft_size/2;
    double rightsum = 0;
    int rightcount = 0;
    double leftsum = 0;
    int leftcount = 0;
    int left_non_missing  = 0;
    int right_non_missing = 0;
    
    int fft_left = 32;
    int fft_right = fft_size-32;
    
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    auto left_it = ordered.begin();
    auto right_it = ordered.end();
    const double alpha = Loess_parms::get_instance().get_alpha();
    const int target_size = ordered.size() * 0.037; // a LOESS q-fraction of 0.035 was necessary for 26.565 degree edges with noise
    int min_left_bin = ordered.front().first*8 + fft_size/2 + 1;
    int max_right_bin = ordered.back().first*8 + fft_size/2 - 1;
    constexpr double max_width = 5.0;
    for (int b=min_left_bin; b < max_right_bin; b++) {
        weights[b] = 1.0;
        
        double mid = (b - fft_size/2)*0.125;
        
        auto mid_it = lower_bound(ordered.begin(), ordered.end(), mid);
        int left_available = mid_it - ordered.begin();
        int right_available = ordered.end() - mid_it;
        
        if (left_available < target_size/2) {
            left_it = ordered.begin();
            right_it = mid_it + (target_size - left_available);
        } else {
            if (right_available < target_size/2) {
                right_it = ordered.end() - 1;
                left_it = mid_it - (target_size - right_available);
            } else {
                left_it = mid_it - target_size/2;
                right_it = mid_it + target_size/2;
            }
        }
        
        constexpr int order = 4;
        size_t npts = right_it - left_it;
        if (npts < (order+1)) {
            printf("empty interval in bin %d\n", b);
            weights[b] = 0;
        } else {
            Eigen::MatrixXd design(npts, order + 1);
            Eigen::VectorXd v(npts);
            size_t row = 0;
            for (auto it=left_it; it != right_it; it++, row++) {
                double d = fabs(it->first - mid);
                double w = exp(-fabs(d)*alpha);
                
                double x = it->first - left_it->first;
                v[row] = w*it->second;
                design(row, 0) = w*1;
                design(row, 1) = w*x;
                design(row, 2) = w*x*x;
                design(row, 3) = w*x*x*x;
                design(row, 4) = w*x*x*x*x;
            }
            Eigen::VectorXd sol = design.colPivHouseholderQr().solve(v);
            
            double ex = mid - left_it->first;
            double ey = sol[0] + ex*sol[1] + ex*ex*sol[2] + ex*ex*ex*sol[3] + ex*ex*ex*ex*sol[4];
            
            mean[b] = ey;
        }
    }
    
    // some housekeeping to take care of missing values
    for (int i=0; i < fft_size; i++) {
        sampled[i] = 0;
    }
    for (int idx=fft_left-1; idx <= fft_right+1; idx++) {
        if (weights[idx] > 0) {
            sampled[idx] = mean[idx] / weights[idx];
            if (idx < fft_size2 - fft_size/8) {
                leftsum += sampled[idx];
                leftcount++;
            }
            if (idx > fft_size2 + fft_size/8) {
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
    }
    
    #if 1
    // estimate a reasonable value for the non-missing samples
    constexpr int nm_target = 8*3;
    int nm_count = 1;
    double l_nm_mean = sampled[left_non_missing];
    for (int idx=left_non_missing+1; idx < fft_size/2 && nm_count < nm_target; idx++) {
        if (sampled[idx] != missing) {
            l_nm_mean += sampled[idx];
            nm_count++;
        }                                                                                                                                                                   
    }                                                                                                                                                                       
    l_nm_mean /= nm_count;
    
    nm_count = 1;
    double r_nm_mean = sampled[right_non_missing];
    for (int idx=right_non_missing-1; idx > fft_size/2 && nm_count < nm_target; idx--) {
        if (sampled[idx] != missing) {
            r_nm_mean += sampled[idx];
            nm_count++;
        }
    }
    r_nm_mean /= nm_count;
    
    // now just pad out the ends of the sequences with the last non-missing values
    for (int idx=left_non_missing-1; idx >= 0; idx--) {
        sampled[idx] = l_nm_mean;
    }
    for (int idx=right_non_missing+1; idx < fft_size; idx++) {
        sampled[idx] = r_nm_mean;
    }

    #else
    // now just pad out the ends of the sequences with the last non-missing values
    for (int idx=left_non_missing-1; idx >= 0; idx--) {
        sampled[idx] = sampled[left_non_missing];
    }
    for (int idx=right_non_missing+1; idx < fft_size; idx++) {
        sampled[idx] = sampled[right_non_missing];
    }
    #endif    
    
    #if 1
    vector<double> smoothed(fft_size, 0);
    // now find 10% / 90% thresholds
    leftsum /= double(leftcount);
    rightsum /= double(rightcount);
    double bright = max(leftsum, rightsum);
    double dark   = min(leftsum, rightsum);
    int p10idx = fft_left-1;
    int p90idx = fft_left-1;
    double p10err = 1e50;
    double p90err = 1e50;
    for (int idx=fft_left; idx <= fft_right; idx++) {
        double smoothed = (sampled[idx-2] + sampled[idx-1] + sampled[idx] + sampled[idx+1] + sampled[idx+2])/5.0;
        if ( fabs((smoothed - dark) - 0.1*(bright - dark)) <  p10err) {
            p10idx = idx;
            p10err = fabs((smoothed - dark) - 0.1*(bright - dark));
        }
        if ( fabs((smoothed - dark) - 0.9*(bright - dark)) <  p90err) {
            p90idx = idx;
            p90err = fabs((smoothed - dark) - 0.9*(bright - dark));
        }
    }
    // we know that mtf50 ~ 1/(p90idx - p10idx) * (1/samples_per_pixel)
    double rise_dist = max(double(4), fabs(double(p10idx - p90idx))*0.125);
    if (p10idx < p90idx) {
        std::swap(p10idx, p90idx);
    }
    p10idx += 4 + 2*lrint(rise_dist); // advance at least one more full pixel
    p90idx -= 4 + 2*lrint(rise_dist);
    int twidth = max(fabs(double(p10idx - fft_size2)), fabs(double(p90idx - fft_size2)));
    constexpr double bwidth = 1.85;
    int left_trans = std::max(fft_size/2 - bwidth*twidth, fft_left + 2.0);
    int right_trans = std::min(fft_size/2 + bwidth*twidth, fft_right - 3.0);
    
    smoothed[0] = sampled[0];
    for (int idx=1; idx < fft_size; idx++) {
        smoothed[idx] = smoothed[idx-1] + sampled[idx];
    }
    constexpr int tpad = 16;
    constexpr int bhw = 16;
    constexpr int bhw_min = 1;
    for (int idx=std::max(fft_left + bhw, left_trans - tpad); idx < left_trans; idx++) {
        int lbhw = (left_trans - idx)*bhw/tpad + bhw_min;
        sampled[idx] = (smoothed[idx+lbhw] - smoothed[idx-lbhw-1])/double(2*lbhw+1);
    }
    for (int idx=std::min(right_trans + tpad - 1, fft_right - bhw - 1); idx > right_trans; idx--) {
        int lbhw = (idx-right_trans)*bhw/tpad + bhw_min;
        sampled[idx] = (smoothed[idx+lbhw] - smoothed[idx-lbhw-1])/double(2*lbhw+1);
    }
    for (int idx=bhw + 1; idx < left_trans - tpad; idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    for (int idx=std::min(right_trans + tpad, fft_right - bhw - 1); idx < fft_size-bhw-1; idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    #endif
    
    /*
    FILE* fo1 = fopen("raw_esf.txt", "wt");
    for (int i=0; i < fft_size; i++) {
        fprintf(fo1, "%d %lf\n", i, sampled[i]);
    }
    fclose(fo1);
    */
    
    int lidx = 0;
    for (int idx=fft_size/4; idx < 3*fft_size/4; idx++) {
        esf[lidx++] = sampled[idx+3];
    }
    
    double old = sampled[fft_left];
    for (int idx=fft_left; idx <= fft_right; idx++) {
        double temp = sampled[idx];
        sampled[idx] = (sampled[idx+1] - old);
        old = temp;
    }
    for (int idx=0; idx < fft_left; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=fft_right; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }
    
    /*
    FILE* fo2 = fopen("raw_psf.txt", "wt");
    for (int i=0; i < fft_size; i++) {
        fprintf(fo2, "%d %lf\n", i, sampled[i]);
    }
    fclose(fo2);
    */
    
    return rval;
}
