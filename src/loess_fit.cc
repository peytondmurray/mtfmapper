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
    
    thread_local vector<double> weights(fft_size);
    thread_local vector<double> mean(fft_size);
    thread_local vector<double> slopes(fft_size);
    thread_local vector<double> smoothed(fft_size);
    
    constexpr double missing = -1e7;
    constexpr double shift_tolerance = 4;

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }
    
    int rval = 0;
    
    constexpr double scale = 2 * 16.0 / 28.0; // maxdot relative to 16
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
    
    int fft_left = fft_size2 + 8*lower;
    int fft_right = fft_size2 + 8*upper;
    
    retry:
    
    for (int i=0; i < int(ordered.size()); i++) {
        int cbin = int(ordered[i].first*8 + fft_size2);
        int left = max(fft_left, cbin-5);
        int right = min(fft_right-1, cbin+5);
        
        for (int b=left; b <= right; b++) {
            double mid = b*scale*(upper-lower)*inv_fft_size_m1 + scale*lower;
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
    int midpoint = fft_size2;
    int left_twidth = fft_size2 - p90idx;
    int right_twidth = p10idx - fft_size2;
    
    // now recompute ESF using box() for tails, and exp() for transition
    rightsum = 0;
    rightcount = 0;
    leftsum = 0;
    leftcount = 0;
    left_non_missing  = 0;
    right_non_missing = 0;
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    int left_trans = fft_size2;
    int right_trans = fft_size2;
    const double nbin_scale = 1.0/std::max(fft_right - fft_size2, fft_size2 - fft_left);
    constexpr double bwidth = 1.85;
    constexpr double lwidth = 0.5;
    left_trans = std::max(midpoint - bwidth*left_twidth, fft_left + 2.0);
    right_trans = std::min(midpoint + bwidth*right_twidth, fft_right - 3.0);
    for (int i=0; i < int(ordered.size()); i++) {
        int cbin = int(ordered[i].first*8 + fft_size2);
        
        int nbins = 5;
        if (cbin > midpoint) {
            if (cbin - midpoint > lwidth*right_twidth) {
                if (cbin - midpoint > 2*lwidth*right_twidth) {
                    nbins = 12;
                } else {
                    nbins = 7;
                }
            }
        } else {
            if (midpoint - cbin > lwidth*left_twidth) {
                if (midpoint - cbin > 2*lwidth*left_twidth) {
                    nbins = 12;
                } else {
                    nbins = 7;
                }
            }
        }
                
        int left = max(fft_left, cbin-nbins);
        int right = min(fft_right-1, cbin+nbins);
        
        if (right < midpoint - bwidth*left_twidth || left > midpoint + bwidth*right_twidth) {
            for (int b=left; b <= right; b++) {
                mean[b] += ordered[i].second;
                weights[b] += 1.0;
            }
        } else {
            for (int b=left; b <= right; b++) {
                double w = 1; // in extreme tails, just plain box filter
                if (b > midpoint) {
                    if (b - midpoint < bwidth*right_twidth) {
                        double mid = b*scale*(upper-lower)*inv_fft_size + scale*lower;
                        if (b - midpoint < lwidth*right_twidth) {
                            // edge transition itself, use preferred low-pass function
                            w = fastexp( -fabs(ordered[i].first - mid)*Mtf_correction::sdev );
                        } else {
                            constexpr double start_factor = 1;
                            constexpr double end_factor =   0.01;
                            double alpha = ((double(b - midpoint))/right_twidth - lwidth)/(bwidth - lwidth);
                            double sfactor = start_factor * (1 - alpha) + end_factor * alpha;
                            // between edge and tail region, use slightly wider low-pass function
                            w = fastexp( -fabs(ordered[i].first - mid)*Mtf_correction::sdev*sfactor );
                        }
                    }
                } else {
                    if (midpoint - b < bwidth*left_twidth) {
                        double mid = b*scale*(upper-lower)*inv_fft_size + scale*lower;
                        if (midpoint - b < lwidth*left_twidth) {
                            // edge transition itself, use preferred low-pass function
                            w = fastexp( -fabs(ordered[i].first - mid)*Mtf_correction::sdev );
                        } else {
                            constexpr double start_factor = 1;
                            constexpr double end_factor =   0.01;
                            double alpha = ((double(midpoint - b))/left_twidth - lwidth)/(bwidth - lwidth);
                            double sfactor = start_factor * (1 - alpha) + end_factor * alpha;
                            // between edge and tail region, use slightly wider low-pass function
                            w = fastexp( -fabs(ordered[i].first - mid)*Mtf_correction::sdev*sfactor );
                        }
                    }
                }
                mean[b] += ordered[i].second * w;
                weights[b] += w;
            }
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
    
    // apply some box smoothing in the tails
    smoothed[fft_left-1] = sampled[fft_left];
    for (int idx=fft_left; idx < fft_right; idx++) {
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
    for (int idx=fft_left + bhw + 1; idx < left_trans - tpad; idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    for (int idx=std::min(right_trans + tpad, fft_right - bhw - 1); idx < fft_right-bhw; idx++) {
        sampled[idx] = (smoothed[idx+bhw] - smoothed[idx-bhw-1])/double(2*bhw+1);
    }
    
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
    
    
    // quickly apply some additional smoothing to the PSF, for good luck
    const int sgh = 2;
    const double sgw[] = {-0.086, 0.343, 0.486, 0.343, -0.086};
    fill(smoothed.begin(), smoothed.end(), 0);
    for (int idx=fft_left+sgh; idx <= fft_right-sgh; idx++) {
        for (int x=-sgh; x <= sgh; x++) {
            smoothed[idx] += sampled[idx+x] * sgw[x+sgh];
        }
    }
    
    for (int idx=0; idx < fft_left; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=fft_left+1; idx < fft_right; idx++) {
        if (abs(idx - fft_size2) <= 8) {
            sampled[idx] = 0.5*(smoothed[idx] + sampled[idx]); // only apply half-strength filtering
        } else {
            sampled[idx] = smoothed[idx];
        }
    }
    for (int idx=fft_right; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }
    for (size_t idx=0; idx < smoothed.size(); idx++) {
        smoothed[idx] = 0;
    }
    
    // apply stronger SG filtering on outer tails
    const int sgh2 = 7;
    const double sgw2[] = {-0.070588, -0.011765, 0.038009, 0.078733, 0.110407, 0.133032, 0.146606, 0.151131, 0.146606, 0.133032, 0.110407, 0.078733, 0.038009, -0.011765, -0.070588}; 
    
    for (int idx=fft_left; idx <= fft_right; idx++) {
        for (int x=-sgh2; x <= sgh2; x++) {
            smoothed[idx] += sampled[idx+x] * sgw2[x+sgh2];
        }
    }
    
    for (int idx=fft_left; idx < left_trans; idx++) {
        sampled[idx] = smoothed[idx];
    }
    for (int idx=left_trans; idx < left_trans+16; idx++) {
        sampled[idx] = 0.5*(smoothed[idx] + sampled[idx]);
    }
    for (int idx=right_trans+1; idx < fft_right; idx++) {
        sampled[idx] = smoothed[idx];
    }
    for (int idx=right_trans+1-16; idx < right_trans; idx++) {
        sampled[idx] = 0.5*(smoothed[idx] + sampled[idx]);
    }
    
    return rval;
}
