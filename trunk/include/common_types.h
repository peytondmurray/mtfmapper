#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <vector>
#include <list>
#include <map>

using std::vector;
using std::pair;
using std::list;
using std::map;

#include <iostream>

using std::cout;
using std::endl;

// OpenCV headers
#include <cv.h>
#include <highgui.h>

typedef cv::Point_<double> Point;
typedef vector<Point> Pointlist;
typedef map<int, Pointlist> Boundarylist;

using std::make_pair;

#include "include/scanline.h"

#define SQR(x) ((x)*(x))

#endif //COMMON_TYPES_H


