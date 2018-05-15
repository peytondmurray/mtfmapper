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
#include <map>
using std::multimap;
using std::pair;
using std::make_pair;

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

static cv::Point2d ransac_line(const vector<double>& data, double il_thresh) {

    double best_score = -1e20;
    cv::Point2d best_sol;
    cv::Point2d sol;
    
    const double s = 1.0/data.size();
    
    multimap<double, cv::Point2d> sols;
    for (size_t i=0; i < data.size()-1; i++) {
        for (size_t j=i+1; j < data.size(); j++) {
        
            sol.x = data.size() * (data[j] - data[i])/(j - i);
            sol.y = data[i];
            
            double sum_x = 0;
            double sum_y = 0;
            double sum_xx = 0;
            double sum_xy = 0;
            double n = 0;
            
            for (size_t k=0; k < data.size(); k++) {
                double pred = k*s*sol.x + sol.y;
                double e = fabs(pred - data[k]);
                if (e < il_thresh) {
                    n++;
                    sum_x += k*s;
                    sum_y += data[k];
                    sum_xy += k*s*data[k];
                    sum_xx += k*s*k*s;
                }
            }
            if (n > 1) {
                double b = (n*sum_xy - sum_x*sum_y)/(n*sum_xx - sum_x*sum_x);
                double a = (sum_y - b*sum_x)/n;
                
                n = 0;
                double rmse = 0;
                for (size_t k=0; k < data.size(); k++) {
                    double pred = k*s*b + a;
                    double e = fabs(pred - data[k]);
                    if (e < il_thresh) {
                        n++;
                        rmse += e*e;
                    }
                }
                rmse = sqrt(rmse/n);
                
                double score = (data.size() - n) + rmse;
                
                if (score > best_score) {
                    best_score = score;
                    best_sol = cv::Point2d(b, a);
                }
                
                sols.insert(make_pair(score, cv::Point2d(b, a)));
                
            } 
        }
    }
    
    return sols.begin()->second;
}
    
static cv::Point_<size_t> region_grow(const vector<double>& data) {
    
    const double il_thresh = 0.005;
    cv::Point2d sol = ransac_line(data, il_thresh);
    
    size_t left_idx = 0;
    size_t right_idx = 0;
    int cstart = -1;
    int cend = -1;
    for (size_t k=0; k < data.size(); k++) {
        double x = k/double(data.size());
        if (fabs(data[k] - (x*sol.x + sol.y)) < 2*il_thresh) {
            if (cstart < 0) {
                cstart = k;
            }
            cend = k;
        } else {
            
            int rlen = cend - cstart;
            if (rlen > int(right_idx - left_idx)) {
                left_idx = cstart;
                right_idx = cend;
                if (rlen > 0.1*data.size()) {
                    break;
                }
            }
            cstart = -1;
            cend = -1;
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
    
    // try using the mean i.s.o Otsu
    double mean = 0;
    for (auto d : data) {
        mean += d;
    }
    double threshold = mean / data.size();
    
    // check the number of transitions
    size_t transitions = 0;
    for (size_t k=0; k < data.size() - 1; k++) {
        if ( (data[k] > threshold && data[k+1] < threshold) ||
             (data[k] < threshold && data[k+1] > threshold) ) {
             
             transitions++;
        }
    }
    
    if (transitions < 5) {
        logger.debug("Imatest_crop: mean thresholding failed. Now trying Otsu\n");
        threshold = otsu_threshold(data);
    }
    
    vector<double> lsegment;
    size_t idx = 0;
    while (idx < data.size() && data[idx] > threshold) idx++;
    size_t lseg_start = idx;
    while (idx < data.size() && data[idx] < threshold) {
        lsegment.push_back(data[idx]);
        idx += 2;
    }
    size_t lseg_end = idx;
    
    if (lsegment.size() < 2) return cv::Point_<int>(0, in_data.size());
    
    cv::Point_<size_t> lseg = region_grow(lsegment);
    logger.debug("Imatest_crop: top segment: [%ld, %ld]\n", lseg_start + 2*lseg.x, lseg_start + 2*lseg.y);
    
    vector<double> rsegment;
    idx = data.size()-1;
    while (idx > lseg_end && data[idx] > threshold) idx--;
    size_t rseg_start = idx;
    while (idx > lseg_end && data[idx] < threshold) {
        rsegment.push_back(data[idx]);
        idx -= 2;
    }
    
    if (rsegment.size() < 2) return cv::Point_<int>(0, in_data.size());
    
    cv::Point_<size_t> rseg = region_grow(rsegment);
    logger.debug("Imatest_crop: bottom segment: [%ld, %ld]\n", rseg_start - 2*rseg.y, rseg_start - 2*rseg.x);
    
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
