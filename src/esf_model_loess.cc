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

#include "include/esf_model_loess.h"
#include "include/fastexp.h"
#include <Eigen/Dense>

// Note: The kernel function is used during ESF construction, but remember that this kernel
// is not used (at all) when calculating MTF corrections
static inline double local_kernel(double d, double alpha, double w=0.125) {
    if (fabs(d) < w) return 1.0;
    return fastexp(-fabs(fabs(d)-w)*alpha);
}

static double interpolate(double cnr, const vector< std::pair<double, double> >& lut) {
    if (cnr < lut.back().first) {
        if (cnr > lut.front().first) {
            size_t i = lut.size() - 1;
            while (i > 0 && cnr < lut[i].first) {
                i--;
            }
            if (i < lut.size() - 1) {
                return lut[i].second + (lut[i+1].second - lut[i].second) * (cnr - lut[i].first) / (lut[i+1].first - lut[i].first);
            } else {
                return lut.back().second;
            }
        } else {
            return lut.front().second;
        }
    } else {
        return lut.back().second;
    }
}

int Esf_model_loess::build_esf(vector< Ordered_point  >& ordered, double* sampled, 
    const int fft_size, double max_distance_from_edge, vector<double>& esf, 
    bool allow_peak_shift) {
        
    thread_local vector<double> weights(fft_size, 0);
    thread_local vector<double> mean(fft_size, 0);

    for (int i=0; i < fft_size; i++) {
        sampled[i] = missing;
    }
    
    int rval = 0;
    
    int fft_left = 0;
    int fft_right = fft_size-1;
    int twidth = 32;
    double cnr = 1;
    double contrast = 1;
    
    rval = estimate_esf_clipping(ordered, sampled, fft_size, allow_peak_shift, 
        max_distance_from_edge, mean, weights, fft_left, fft_right, twidth, cnr, contrast
    );
    
    if (rval < 0) return rval;
    
    constexpr double bwidth = 1.0;
    int left_trans = std::max(fft_size/2 - bwidth*twidth, fft_left + 2.0);
    int right_trans = std::min(fft_size/2 + bwidth*twidth, fft_right - 3.0);
    
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    
    auto left_it = ordered.begin();
    auto right_it = ordered.end();
    double alpha = get_alpha();
    constexpr double model_switch_bias = 3*0.125;
    int min_left_bin = ordered.front().first*8 + fft_size/2 + 1;
    int max_right_bin = ordered.back().first*8 + fft_size/2 - 1;
    
    // TODO: I still see a lot of noise in the transition area (0.75*twidth*0.125 up to 1.0*twidth + model_switch_bias)
    // so some additional smoothing in this area will help
    
    vector< std::pair<double, double> > ridge_lut = {
        {20, 250},{22, 231},{24, 190},{26, 163},{28, 143},{30, 130},{35, 106},{40, 93},{45, 87},
        {50, 78},{60, 67},{70, 60},{80, 55},{90, 51},{99, 47},{124, 42},{149, 37},{174, 34},{199, 32},
        {298, 26},{398, 13},{497, 8},{596, 2},{696, 2},{794, 1},{894, 1},{1000, 5e-8}
    };
    const double ridge_parm = interpolate(cnr, ridge_lut);
    
    constexpr int order = 6;
    for (int b=min_left_bin; b < max_right_bin; b++) {
        double mid = (b - fft_size/2)*0.125;
        
        constexpr double loess_span = 4.5;
        left_it = lower_bound(ordered.begin(), ordered.end(), mid - 0.5*loess_span);
        right_it = lower_bound(left_it, ordered.end(), mid + 0.5*loess_span);
        
        size_t npts = right_it - left_it;
        if (npts < (order+1)) {
            printf("empty interval in bin %d\n", b);
            weights[b] = 0;
        } else {
            if (fabs(mid) < 0.125*twidth + model_switch_bias) {
                Eigen::MatrixXd design(npts, order + 1);
                Eigen::VectorXd v(npts);
                
                double fw = 0.125;
                if (fabs(mid) >= 0.125*0.5*twidth) {
                    fw = 0.250;
                }
                
                size_t row = 0;
                for (auto it=left_it; it != right_it; it++, row++) {
                    double d = fabs(it->first - mid);
                    double w = local_kernel(d, alpha, fw); 
                    
                    double x = (it->first - mid)/(0.5*loess_span);
                    v[row] = w*(it->second - mean[b-1])/contrast;
                    double x2 = x*x;
                    double x3 = x2*x;
                    double x4 = x2*x2;
                    double x5 = x4*x;
                    double x6 = x4*x2;
                    
                    design(row, 0) = w*1;
                    design(row, 1) = w*(x);
                    design(row, 2) = w*(2*x2 - 1);
                    design(row, 3) = w*(4*x3 - 3*x);
                    design(row, 4) = w*(8*x4 - 8*x2 + 1);
                    design(row, 5) = w*(16*x5 - 20*x3 + 5*x);
                    design(row, 6) = w*(32*x6 - 48*x4 + 18*x2 - 1);
                }
                const double phi = (fabs(mid) > 0.75*twidth*0.125) ? ridge_parm : 5e-8;
                Eigen::VectorXd sol = (design.transpose() * design + phi*Eigen::MatrixXd::Identity(order+1, order+1)).llt().solve(design.transpose() * v);
                
                mean[b] = (sol[0] + sol[4])*contrast - (sol[2] + sol[6])*contrast + mean[b-1];
                weights[b] = 1.0;
            } else {
                // blend in the box filter, otherwise f/16 LSFs lose some detail right at the foot
                const double upper_thresh = 2.0*0.125*twidth;
                const double lower_thresh = 0.125*twidth;
                constexpr double min_span = 0.126; // such a narrow interval can produce empty bins, handled below
                constexpr double max_span = 1.0;
                double span = min_span;
                if (fabs(mid) >= upper_thresh) {
                    span = max_span;
                } else {
                    double t = (fabs(mid) - lower_thresh)/(upper_thresh - lower_thresh);
                    span = min_span + t*(max_span - min_span);
                }
                left_it = lower_bound(left_it, right_it, mid - span);
                right_it = lower_bound(left_it, right_it, mid + span);
                
                double sum = 0;
                for (auto it=left_it; it != right_it; it++) {
                    sum += it->second;
                }
                if (right_it - left_it > 0) {
                    mean[b] = sum / double(right_it - left_it);
                    weights[b] = 1.0;
                } else {
                    mean[b] = mean[b-1];
                    weights[b] = weights[b-1];
                }
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
    
    vector< std::pair<double, double> > smoothing_lut = {
        {4, 30},{5, 30},{6, 30},{7, 30},{8, 30},{9, 30},{10, 30},{12, 30},{14, 30},{16, 30},{18, 30},
        {20, 30},{22, 30},{24, 30},{26, 30},{28, 30},{30, 30},{35, 30},{40, 30},{45, 30},{49, 30},
        {59, 30},{69, 30},{79, 30},{89, 21},{99, 21},{124, 19},{148, 14},{173, 14},{198, 14},
        {297, 13},{395, 6},{493, 6},{591, 1},{690, 1},{787, 0},{885, 0}
    };
    
    int sw_width = lrint(interpolate(cnr, smoothing_lut));
    if (sw_width > 0) {
        moving_average_smoother(smoothed, sampled, fft_size, fft_left, fft_right, left_trans, right_trans, sw_width);
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
    for (int idx=0; idx < fft_left+2; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=fft_right-3; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }
    
    return rval;
}

