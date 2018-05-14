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

#include "include/imatest_crop.h"

Imatest_cropper::Imatest_cropper(cv::Mat& in_img) : Cropper(in_img) {
    double maxval = -1e20;
    for (int r=2; r < in_img.rows-2; r+=2) {
        for (int c=4; c < in_img.cols-4; c+=2) {
            int_rs[r] += in_img.at<uint16_t>(r, c) / 65536.0;
            maxval = std::max(maxval, int_rs[r]);
            int_rs[r+1] = int_rs[r];
        }
    }
    int_rs[0] = int_rs[1] = int_rs[2];
    int_rs[in_img.rows-1] = int_rs[in_img.rows-2] = int_rs[in_img.rows-3];

    cv::Point_<int> bounds = black_bar_bounds(int_rs);
    rstart = bounds.x;
    height = bounds.y - bounds.x;
    
    cstart = 0;
    width  = in_img.cols;
    max_brightness = lrint(2*maxval/in_img.cols*65535);
}
    
static cv::Point_<size_t> region_grow(const vector<double>& data) {
    class Running {
      public:
        Running(double mean, double min_sdev) : mean(mean), prev_mean(mean), min_sdev(min_sdev) {}
        
        void add(double x) {
            n++;
            mean = prev_mean + (x - prev_mean)/n;
            var = prev_var + (x - prev_mean)*(x - mean);
            
            prev_mean = mean;
            prev_var = var;
        }
        
        double get_mean(void) const { return mean; };
        
        double get_sdev(void) const { 
            return std::max(sqrt(var/(n-1)), min_sdev); 
        };
        
        double mean;
        double prev_mean;
        
        size_t n = 1;
        double prev_var = 0;
        double var = 0;
        double min_sdev;
    };
    
    vector<double> lscopy(data);
    sort(lscopy.begin(), lscopy.end());
    double median = lscopy[lscopy.size()/2];
    
    size_t med_idx = 0;
    for (size_t k=0; k < data.size(); k++) {
        if (data[k] == median) {
            med_idx = k;
        }
    }
    
    Running m(median, 1.0/600.0);
    
    for (size_t k=lscopy.size()/2-2; k <= lscopy.size()/2 + 2; k++) {
        m.add(lscopy[k]);
    }
    
    size_t left_idx = med_idx;
    size_t right_idx = med_idx;
    
    bool done = false;
    while (!done) {
        done = true;
        if (left_idx > 0 && fabs(data[left_idx - 1] - m.get_mean()) < 3*m.get_sdev()) {
            left_idx--;
            m.add(data[left_idx]);
            done = false;
        }
        if (right_idx < data.size()-1 && fabs(data[right_idx + 1] - m.get_mean()) < 3*m.get_sdev()) {
            right_idx++;
            m.add(data[right_idx]);
            done = false;
        }
    }
    
    return cv::Point_<size_t>(left_idx, right_idx);
}
    
cv::Point_<int> Imatest_cropper::black_bar_bounds(const vector<double>& in_data) {

    vector<double> data(in_data);
    double maxv = -1e20;
    for (auto d : in_data) {
        maxv = std::max(maxv, d);
    }
    for (auto& d : data) {
        d /= maxv;
    }

    double otsu = otsu_threshold(data);
    
    vector<double> lsegment;
    size_t idx = 0;
    while (idx < data.size() && data[idx] > otsu) idx++;
    size_t lseg_start = idx;
    while (idx < data.size() && data[idx] < otsu) {
        lsegment.push_back(data[idx]);
        idx += 2;
    }
    size_t lseg_end = idx;
    
    cv::Point_<size_t> lseg = region_grow(lsegment);
    
    vector<double> rsegment;
    idx = data.size()-1;
    while (idx > lseg_end && data[idx] > otsu) idx--;
    size_t rseg_start = idx;
    while (idx > lseg_end && data[idx] < otsu) {
        rsegment.push_back(data[idx]);
        idx -= 2;
    }
    
    cv::Point_<size_t> rseg = region_grow(rsegment);
    
    int even_start = lseg_start + 2*lseg.x;
    even_start += even_start % 2;
    int even_end = rseg_start - 2*rseg.x;
    even_end -= even_end % 2;
    
    return cv::Point_<int>(even_start, even_end);
}

void Imatest_cropper::fill_bars(cv::Mat& X) {
    
    cv::Mat top_block(X, cv::Rect(0, 0, X.cols, rstart));
    top_block = get_max_brightness();
    
    cv::Mat bot_block(X, cv::Rect(0, rstart+height, X.cols, X.rows - (rstart+height)));
    bot_block = get_max_brightness();
}
