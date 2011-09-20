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
#include "include/gradient.h"
#include "include/common_types.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//------------------------------------------------------------------------------
Gradient::Gradient(const cv::Mat& in_img, bool snapshot)
 : _width(in_img.cols), _height(in_img.rows)
{
    
    cv::Mat in_float;
    in_img.convertTo(in_float, CV_32FC1, 1.0/65536.0);
    
    cv::Mat smoothed;
    cv::GaussianBlur(in_float, smoothed, cv::Size(5,5), 1.0, 1.0);
    
    _compute_gradients(smoothed);

    _gradient_m = cv::Mat(in_img.rows, in_img.cols, CV_32FC1);
    for (int y=0; y < in_img.rows; y++) {
        for (int x=0; x < in_img.cols; x++) {
            _gradient_m.at<float>(y,x) = sqrt(SQR(_gradient_x.at<float>(y,x)) + SQR(_gradient_y.at<float>(y,x)));
        }
    }

    if (snapshot) {
        imwrite(std::string("gradient.tiff"), _gradient_m*255.0);
    }
}

//------------------------------------------------------------------------------
Gradient::~Gradient(void) {
}


//------------------------------------------------------------------------------
void Gradient::_compute_gradients(const cv::Mat& smoothed_im) {

    _gradient_x = cv::Mat(smoothed_im.rows, smoothed_im.cols, CV_32FC1);
    _gradient_y = cv::Mat(smoothed_im.rows, smoothed_im.cols, CV_32FC1);

    for(int r=0; r < smoothed_im.rows; r++){
        for(int c=1; c < smoothed_im.cols - 1; c++){
            _gradient_x.at<float>(r,c) = smoothed_im.at<float>(r, c + 1) - smoothed_im.at<float>(r, c - 1);
        }
    }
    for(int r=1; r < smoothed_im.rows - 1; r++){
        for(int c=0; c < smoothed_im.cols; c++){
            _gradient_y.at<float>(r,c) = smoothed_im.at<float>(r + 1, c) - smoothed_im.at<float>(r - 1, c);
        }
    }


}

