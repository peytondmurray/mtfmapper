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

#include "include/esf_model_hann.h"
#include "include/afft.h"
#include "include/pava.h"

// Note: we avoid a virtual call by explicitly calling an inline version of the 
// kernel function
static inline double local_kernel(double x) {
    constexpr double c = 0.5/8.0;
    return fabs(x) < c ? 1.0 : 0.0;
}

// Note: this is the generic version of the kernel function that we end up
// using to construct the MTF correction table
double Esf_model_hann::kernel(double x) const {
    return local_kernel(x);
}

int Esf_model_hann::build_esf(vector< Ordered_point  >& ordered, double* sampled, 
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
    
    const int fft_size2 = fft_size/2;
    
    fill(weights.begin(), weights.end(), 0.0);
    fill(mean.begin(), mean.end(), 0.0);
    
    for (int i=0; i < int(ordered.size()); i++) {
        int cbin = int(ordered[i].first*8 + fft_size2);
        
        int nbins = 5;
        int left = max(fft_left, cbin-nbins);
        int right = min(fft_right-1, cbin+nbins);
        
        for (int b=left; b <= right; b++) {
            double mid = (b - fft_size2)*0.125;
            double w = local_kernel(ordered[i].first - mid);
            mean[b] += ordered[i].second * w;
            weights[b] += w;
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
            sampled[idx] = mean[idx] / weights[idx];
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
    
    // apply Masaoka's interpolation method to take care of zero-count bins
    for (int idx=fft_left; idx < fft_right; idx++) {  
        if (fabs(sampled[idx] - missing) < 1e-6) {
            int prev_present = idx - 1;
            int next_present = idx;
            while (fabs(sampled[next_present] - missing) < 1e-6) {
                next_present++;
            }
            int j;
            for (j=prev_present + 1; j < next_present - 1; j++) {
                sampled[j] = sampled[prev_present];
            }
            sampled[j] = 0.5*(sampled[prev_present] +  sampled[next_present]);
        }         
    }
    
    if (apply_monotonic_filter) {
        vector<double> pre_pava(fft_size);
        for (size_t i=0; i < size_t(fft_size); i++) {
            pre_pava[i] = sampled[i];
        }
        
        // first, reverse the vector if necessary
        bool reversed = false;
        if (l_nm_mean > r_nm_mean) {
            reversed = true;
            std::reverse(sampled, sampled + fft_size);
        }
        jbk_pava(sampled, fft_size);
        if (reversed) {
            std::reverse(sampled, sampled + fft_size);
        }
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
    for (int idx=0; idx < fft_left; idx++) {
        sampled[idx] = 0;
    }
    for (int idx=fft_right; idx < fft_size; idx++) {
        sampled[idx] = 0;
    }

    // apply Hann window to LSF
    if (window_function.size() != size_t(fft_size)) {
        fprintf(stderr, "Someone changed the FFT_SIZE, and forgot to tell me!\n");
        exit(1);
    }
    for (int idx=0; idx < fft_size; idx++) {
        sampled[idx] *= window_function[idx];
    }
    
    return rval;        
}

void Esf_model_hann::compute_mtf_corrections(void) {
    // obtain derivative correction from parent
    Esf_model::compute_mtf_corrections();
    
    constexpr int c_fft_size = 512*16;
    vector<double> fft_buffer(c_fft_size);
    AFFT<c_fft_size> afft;
    
    for (int i=0; i < c_fft_size; i++) {
        double x = (i - c_fft_size/2)/(16*8.0);
        fft_buffer[i] = fabs(x) <= 0.625 ? kernel(x) : 0;
    }
    
    afft.realfft(fft_buffer.data());
    
    vector<double> correction(w.size());
    double n0 = fabs(fft_buffer[0]);
    correction[0] = 1.0;
    for (int i=1; i < NYQUIST_FREQ*4; i++) {
        correction[i] = sqrt(SQR(fft_buffer[i]) + SQR(fft_buffer[c_fft_size - i])) / n0;
    }
    
    // combine kernel correction with derivative correction
    for (int i=1; i < NYQUIST_FREQ*4; i++) {
        w[i] *= correction[i];
    }
    
    // just guess that the window function will be 512 samples in length
    constexpr int N = 512;
    window_function.resize(N);
    for (int i=0; i < N; i++) {
        double theta = 2*M_PI*i/double(N);
        window_function[i] = 0.5 * (1.0 - cos(theta));
    }
}