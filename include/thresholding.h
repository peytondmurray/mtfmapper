#ifndef THRESHOLDING_H
#define THRESHOLDING_H

#include "opencv2/opencv.hpp"

void bradley_adaptive_threshold(const cv::Mat& cvimg, cv::Mat& img, double threshold, int S);

#endif // THRESHOLDING_H


