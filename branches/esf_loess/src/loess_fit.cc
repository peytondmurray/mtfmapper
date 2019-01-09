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
constexpr double missing = -1e7;

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

int estimate_esf_clipping(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, bool allow_peak_shift, 
    int effective_maxdot, vector<double>& mean, vector<double>& weights,
    int& fft_left, int& fft_right, int& twidth) {
    
    thread_local vector<double> slopes(fft_size, 0);
    
    constexpr double shift_tolerance = 4;

    int rval = 0;

    const int fft_size2 = fft_size/2;
    const double inv_fft_size_m1 = 1.0/double(fft_size-1);
    const double inv_fft_size = 1.0/double(fft_size);
    double rightsum = 0;
    int rightcount = 0;
    double leftsum = 0;
    int leftcount = 0;
    int left_non_missing  = 0;
    int right_non_missing = 0;
    fill(weights.begin(), weights.end(), 0);
    fill(mean.begin(), mean.end(), 0);
    
    int clipped_count = 0;
    
    fft_left = fft_size2 - 8*effective_maxdot;
    fft_right = fft_size2 + 8*effective_maxdot;
    
    //printf("initial: fl:%d, fr:%d ", fft_left, fft_right);
    
    retry:
    
    for (int i=0; i < int(ordered.size()); i++) {
        int cbin = int(ordered[i].first*8 + fft_size2);
        int left = max(fft_left, cbin-5);
        int right = min(fft_right-1, cbin+5);
        
        for (int b=left; b <= right; b++) {
            double mid = (b - fft_size2)*0.125;
            double w = 1 - abs((ordered[i].first - mid)*1.75) > 0 ? 1 - abs((ordered[i].first - mid)*1.75) : 0;
            mean[b] += ordered[i].second * w;
            weights[b] += w;
        }
    }
    // some housekeeping to take care of missing values
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
    
    // now just pad out the ends of the sequences with the last non-missing values
    for (int idx=left_non_missing-1; idx >= 0; idx--) {
        sampled[idx] = sampled[left_non_missing];
    }
    for (int idx=right_non_missing+1; idx < fft_size; idx++) {
        sampled[idx] = sampled[right_non_missing];
    }
    
    fill(slopes.begin(), slopes.end(), 0);
    int peak_slope_idx = 0;
    constexpr int pw = 16;
    for (int idx=pw+1; idx < fft_size-1-pw; idx++) {
        double sx2 = 0;
        double sxy = 0;
        // sx == 0 because the window is symmetric
        for (int j=-pw; j <= pw; j++) {
            sxy += sampled[idx+j] * j;
            sx2 += j*j;
        }
        slopes[idx] = sxy/(sx2);
        if (fabs(slopes[idx]) > fabs(slopes[peak_slope_idx]) && idx > fft_left+pw && idx < fft_right-pw-1) {
            peak_slope_idx = idx;
        }
    }
    
    if (!allow_peak_shift) {
        if (abs(peak_slope_idx - fft_size / 2) > 2 * 8 &&    // peak is more than 2 pixels from centre
            abs(peak_slope_idx - fft_size / 2) < 12 * 8) { // but not at the edge?
            logger.debug("edge rejected because of shifted peak slope: %lf\n", abs(peak_slope_idx - fft_size / 2) / 8.0);
            return -1;
        }
    }
    
    // compute central peak magnitude and sign
    double central_peak = 0;
    for (int w=-16; w <= 16; w++) {
        if (fabs(slopes[fft_size2+w]) > fabs(central_peak)) {
            central_peak = slopes[fft_size2+w];
        }
    }
    
    double peak_threshold = fabs(central_peak * 0.001); // must be negative by at least a little bit
    // scan for significant slope sign change
    bool clipped = false;
    for (int idx=fft_size2-16; idx >= fft_left+4; idx--) {
        if (slopes[idx]*central_peak < 0 && fabs(slopes[idx]) > peak_threshold) {
            // check if a fair number of the remaining slopes are also negative (relative to peak)
            int below = 0;
            double maxdev = 0;
            int scount = 0;
            for (int j=idx; j >= fft_left; j--) {
                if (slopes[j]*central_peak < 0) {
                    below++;
                    maxdev = max(maxdev, fabs(slopes[j]));
                }
                scount++;
            }
            if ((below > scount*0.4 && maxdev/fabs(central_peak) > 0.25) || (below > 0.9*scount && scount > 16)) {
                fft_left = min(idx, fft_size2 - 2*8);
                clipped = true;
                break;
            } 
        }
    }
    for (int idx=fft_size2+16; idx < fft_right-4; idx++) {
        if (slopes[idx]*central_peak < 0 && fabs(slopes[idx]) > peak_threshold) {
            // check if a fair number of the remaining slopes are also negative (relative to peak)
            int below = 0;
            double maxdev = 0;
            int scount=0;
            for (int j=idx; j < fft_right; j++) {
                if (slopes[j]*central_peak < 0) {
                    below++;
                    maxdev = max(maxdev, fabs(slopes[j]));
                }
                scount++;
            }
            if ((below > scount*0.4 && maxdev/fabs(central_peak) > 0.25) || (below > 0.9*scount && scount > 16)) {
                fft_right = max(idx, fft_size2 + 2*8);
                clipped = true;
                break;
            } 
        }
    }
    
    
    if (clipped) {
        if (fft_size2 - fft_left < shift_tolerance*8 ||
            fft_right  - fft_size2 < shift_tolerance*8) {
            
            logger.debug("probably contamination. tagging edge as dodgy\n");
            rval = 1;
        }
    }
    
    if (clipped && clipped_count < 2) {
        for (size_t idx=0; idx < weights.size(); idx++) {
            sampled[idx] = 0;
            weights[idx] = 0;
        }
        leftsum = 0;
        rightsum = 0;
        rightcount = 0;
        leftcount = 0;
        left_non_missing = 0;
        right_non_missing = 0;
        clipped_count++;
        goto retry;
    }
    
    //printf("final: fl:%d, fr:%d \n", fft_left, fft_right);
    
    leftsum /= double(leftcount);
    rightsum /= double(rightcount);
    
    // now find 10% / 90% thresholds
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
    twidth = max(fabs(double(p10idx - fft_size2)), fabs(double(p90idx - fft_size2)));
    
    return rval;
}

void moving_average_smoother(vector<double>& smoothed, double* sampled, 
    int fft_size, int fft_left, int fft_right, int left_trans, int right_trans) {
    
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
}

inline double kernel(double d, double alpha, double w=0.125) {
    if (fabs(d) < w) return 1.0;
    return fastexp(-fabs(fabs(d)-w)*alpha);
}

int bin_fit(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, double lower, double upper, vector<double>& esf, bool allow_peak_shift) {
    
    thread_local vector<double> weights(fft_size, 0);
    thread_local vector<double> mean(fft_size, 0);

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }
    
    int rval = 0;
    
    int fft_left = 0;
    int fft_right = fft_size-1;
    int twidth = 32;
    
    rval = estimate_esf_clipping(ordered, sampled, fft_size, allow_peak_shift, 
        fabs(upper), mean, weights, fft_left, fft_right, twidth
    );
    
    if (rval < 0) return rval;
    
    constexpr double bwidth = 1.85;
    int left_trans = std::max(fft_size/2 - bwidth*twidth, fft_left + 2.0);
    int right_trans = std::min(fft_size/2 + bwidth*twidth, fft_right - 3.0);
    
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    
    auto left_it = ordered.begin();
    auto right_it = ordered.end();
    const double alpha = Loess_parms::get_instance().get_alpha();
    int min_left_bin = ordered.front().first*8 + fft_size/2 + 1;
    int max_right_bin = ordered.back().first*8 + fft_size/2 - 1;
    
    constexpr int order = 4;
    
    for (int b=min_left_bin; b < max_right_bin; b++) {
        weights[b] = 1.0;
        
        double mid = (b - fft_size/2)*0.125;
        
        int target_size = ordered.size() * 0.037;
        if (fabs(b - fft_size/2) > bwidth*twidth) {
            target_size = ordered.size() * 0.075;
        } else {
            if (fabs(b - fft_size/2) > twidth) {
                target_size = ordered.size() * 0.055;
            }
        }
        
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
        
        size_t npts = right_it - left_it;
        if (npts < (order+1)) {
            printf("empty interval in bin %d\n", b);
            weights[b] = 0;
        } else {
            if (fabs(mid) < 0.125*twidth) {
            
                Eigen::MatrixXd design(npts, order + 1);
                Eigen::VectorXd v(npts);
                
                double fw = 0.125*0.5;
                if (fabs(mid) >= 0.125*0.5*twidth) {
                    fw = 0.125;
                }
                
                size_t row = 0;
                for (auto it=left_it; it != right_it; it++, row++) {
                    double d = fabs(it->first - mid);
                    double w = kernel(d, alpha, fw);
                    
                    double x = it->first - left_it->first;
                    v[row] = w*it->second;
                    design(row, 0) = w*1;
                    design(row, 1) = w*x;
                    design(row, 2) = w*x*x;
                    design(row, 3) = w*x*x*x;
                    design(row, 4) = w*x*x*x*x;
                }
                Eigen::VectorXd sol = (design.transpose() * design).llt().solve(design.transpose() * v);
                
                double ex = mid - left_it->first;
                double ey = sol[0] + ex*sol[1] + ex*ex*sol[2] + ex*ex*ex*sol[3] + ex*ex*ex*ex*sol[4];
                
                mean[b] = ey;
            } else {
                double sum = 0;
                double count = 0;
                for (auto it=left_it; it != right_it; it++) {
                    double d = fabs(it->first - mid);
                    double w = kernel(d, alpha, 0.25);
                    sum += w*it->second;
                    count += w;
                }
                mean[b] = sum / count;
            }
        }
    }
    
    int left_non_missing = 0;
    int right_non_missing = 0;
    
    // some housekeeping to take care of missing values
    for (int i=0; i < fft_size; i++) {
        sampled[i] = 0;
    }
    for (int idx=fft_left-1; idx <= fft_right+1; idx++) {
        if (weights[idx] > 0) {
            sampled[idx] = mean[idx];
            if (!left_non_missing) {
                left_non_missing = idx; // first non-missing value from left
            }
            right_non_missing = idx; // last non-missing value
        } else {
            sampled[idx] = missing;
        }
    }
    
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

    thread_local vector<double> smoothed(fft_size, 0);
    moving_average_smoother(smoothed, sampled, fft_size, fft_left, fft_right, left_trans, right_trans);
    
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
    
    return rval;
}
