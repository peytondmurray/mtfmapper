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

#ifndef UNDISTORT_H
#define UNDISTORT_H

#include "common_types.h"
#include "include/logger.h"
#include <opencv2/imgproc/imgproc.hpp>

#include <map>
using std::map;

class Undistort { 
  public:
    Undistort(const cv::Rect& r) : centre(r.width/2, r.height/2), offset(r.x, r.y) {};
    
    cv::Point2i transform_pixel(int col, int row) {
        cv::Point2d tp = transform_point(double(col), double(row));
        return cv::Point2i(lrint(tp.x), lrint(tp.y));
    }
    
    virtual cv::Point2d transform_point(double px, double py) = 0;
    virtual cv::Point2d inverse_transform_point(double px, double py) = 0;
    
    virtual cv::Mat unmap(const cv::Mat& src, cv::Mat& rawimg) = 0;
    
    bool rectilinear_equivalent(void) const {
        return rectilinear;
    }
    
    void set_rectilinear_equivalent(bool b) {
        rectilinear = b;
    }
    
    
    cv::Point2d centre;
    cv::Point2d offset;
    
    bool rectilinear = false;
};

class Undistort_equiangular : public Undistort {
  public:
    Undistort_equiangular(const cv::Rect& r, double f, double pitch) : Undistort(r), f(f), pitch(pitch), rect_f(f) {}
    
    cv::Point2d transform_point(double col, double row) {
        double px = (col + offset.x - centre.x)*pitch;
        double py = (row + offset.y - centre.y)*pitch;
    
        double rd = sqrt((px)*(px) + (py)*(py)); // radial distance in mm
        double theta = atan(rd/rect_f);
        double ru = theta*f;
        
        px = px*ru/rd / pitch + centre.x - offset.x;
        py = py*ru/rd / pitch + centre.y - offset.y;
        
        return cv::Point2d(px, py);
    }
    
    cv::Point2d inverse_transform_point(double col, double row) {
        double px = (col + offset.x - centre.x)*pitch;
        double py = (row + offset.y - centre.y)*pitch;
    
        double rd = sqrt((px)*(px) + (py)*(py)); // radial distance in mm
        double theta = rd/f;
        double ru = tan(theta)*rect_f;
        
        px = px*ru/rd / pitch + centre.x - offset.x;
        py = py*ru/rd / pitch + centre.y - offset.y;
        
        return cv::Point2d(px, py);
    }
    
    // note: second parameter is the raw Bayer image, which must also be padded out
    cv::Mat unmap(const cv::Mat& in_src, cv::Mat& rawimg) {
        
        Point2d extreme = inverse_transform_point(0, in_src.rows/2);
        int pad_left = extreme.x < 0 ? ceil(-extreme.x) : 0;
        
        extreme = inverse_transform_point(in_src.cols/2, 0);
        int pad_top  = extreme.y < 0 ? ceil(-extreme.y) : 0;
        
        int pad = std::max(pad_left, pad_top);
        pad_left = pad;
        pad_top = pad;
        
        // to preserve Bayer CFA alignment
        pad_left += pad_left % 2;
        pad_top += pad_top % 2;
        
        cv::Mat src;
        copyMakeBorder(in_src, src, pad_top, pad_top, pad_left, pad_left, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        centre.x = src.cols / 2;
        centre.y = src.rows / 2;
        
        logger.debug("Padding with %d left/right pixels, %d top/bottom pixels\n", pad_left, pad_top);
        
        if (pad_left > 0 || pad_top > 0) {
            logger.debug("Padding distorted/Bayer image\n");
            cv::Mat rcopy = rawimg.clone();
            copyMakeBorder(rcopy, rawimg, pad_top, pad_top, pad_left, pad_left, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        }
    
        cv::Mat map_x(src.rows, src.cols, CV_32FC1);
        cv::Mat map_y(src.rows, src.cols, CV_32FC1);
        
        double maxrad = sqrt(src.cols*src.cols/4 + src.rows*src.rows/4);
        vector<double> radmap;
        for (int rad=0; rad <= 2*(maxrad + 2); rad++) { // half-pixel spacing
            double theta = atan(double(rad*0.5)*pitch/rect_f);
            double r_d = theta*f;
            radmap.push_back(r_d/pitch);
        }
        
        for (int r=0; r < src.rows; r++) {
            for (int c=0; c < src.cols; c++) {
                
                double px = (c + offset.x - centre.x);
                double py = (r + offset.y - centre.y);
                
                double rad = 2*sqrt(px*px + py*py); // convert to half-pixel spacing
                int rad_f = (int)rad;
                
                double w = rad - rad_f;
                // quadratic interpolation through half-pixel spaced points
                double fc = radmap[rad_f];
                double fb = -1.5*radmap[rad_f] + 2*radmap[rad_f+1] - 0.5*radmap[rad_f+2];
                double fa =  0.5*radmap[rad_f] - radmap[rad_f+1] + 0.5*radmap[rad_f+2];
                double rad_d = 2*(fa*w*w + fb*w + fc);
                
                map_x.at<float>(r, c) = px/rad*rad_d + centre.x - offset.x;
                map_y.at<float>(r, c) = py/rad*rad_d + centre.y - offset.y;
            }
        }
        cv::Mat timg;
        cv::remap(src, timg, map_x, map_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        
        return timg;
    }
    
    double f;
    double pitch;
    double rect_f;
    
};
    
#endif
    
