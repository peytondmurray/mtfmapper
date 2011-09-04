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
#include <string>

using std::cout;
using std::endl;

// OpenCV headers
#include <cv.h>
#include <highgui.h>

typedef cv::Point_<double> Point;
typedef vector<Point> Pointlist;
typedef map<int, Pointlist> Boundarylist;

using std::make_pair;

#ifdef _WIN32
	#define EXE_SUFFIX ".exe"
#else
	#define EXE_SUFFIX ""
#endif

#ifdef _MSC_VER
	typedef unsigned short int uint16_t;

	#define M_PI 3.14159265358979
	#define lrint(x) ( (x < 0) ? int(floor(x-0.5)) : int(floor(x+0.5)) )
#endif

#include "include/scanline.h"

#define SQR(x) ((x)*(x))

#endif //COMMON_TYPES_H


