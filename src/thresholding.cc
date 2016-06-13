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
#include <math.h>

#include "include/thresholding.h"
#include "include/common_types.h"

#include <stdint.h>

// D. Bradley, G. Roth. ACM Journal of Graphics Tools. 2007. Vol 12, No. 2: 13-21.
//
//------------------------------------------------------------------------------
void bradley_adaptive_threshold(const cv::Mat& cvimg, cv::Mat& out, double threshold, int S) {

    out = cv::Mat(cvimg.rows, cvimg.cols, CV_8UC1);

    uint64_t* integralImg = 0;
    int i, j;
    int64_t sum=0;
    int64_t count=0;
    int index;
    int x1, y1, x2, y2;
    int s2 = S/2;

    // create the integral image
    integralImg = new uint64_t[out.rows*out.cols];


    for (i=0; i < out.cols; i++) {
        // reset this column sum
        sum = 0;

        for (j=0; j < out.rows; j++) {
            index = j*out.cols+i;

            sum += int64_t(cvimg.at<uint16_t>(j, i));
            if (i==0) {
                integralImg[index] = sum;
            } else {
                integralImg[index] = integralImg[index-1] + sum;
            }
        }
    }

    // perform thresholding
    cv::MatConstIterator_<uint16_t> it = cvimg.begin<uint16_t>(); 
    for (j=0; j < out.rows; j++) {
        for (i=0; i < out.cols; i++) {
            index = j*out.cols+i;

            // set the SxS region
            x1=i-s2; x2=i+s2;
            y1=j-s2; y2=j+s2;

            // check the border
            x1 = max(0, x1);
            x2 = min(x2, out.cols - 1);
            y1 = max(0, y1);
            y2 = min(y2, out.rows -1);

            count = (x2-x1)*(y2-y1);

            // I(x,y)=s(x2,y2)-s(x1,y2)-s(x2,y1)+s(x1,x1)
            sum = integralImg[y2*out.cols+x2] -
                  integralImg[y1*out.cols+x2] -
                  integralImg[y2*out.cols+x1] +
                  integralImg[y1*out.cols+x1];

            if (((int64_t)(*it)*count) < (int64_t)(sum*(1.0-threshold))) {
                out.data[index] = 0;
            } else {
                out.data[index] = 255;
            }
            ++it;
        }
    }

    delete [] integralImg;
}

