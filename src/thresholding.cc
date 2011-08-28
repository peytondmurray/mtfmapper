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

        for (j=0; j<out.rows; j++) {
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
    for (i=0; i < out.cols; i++) {
        for (j=0; j < out.rows; j++) {
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

            if (((int64_t)(cvimg.at<uint16_t>(j, i))*count) < (int64_t)(sum*(1.0-threshold))) {
                out.data[index] = 0;
            } else {
                out.data[index] = 255;
            }
        }
    }

    delete [] integralImg;
}

