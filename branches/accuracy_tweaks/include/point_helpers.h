#ifndef POINT_HELPERS_H
#define POINT_HELPERS_H

#include "include/common_types.h"
#include "include/gradient.h"

Point centroid(const Pointlist& p);
Point average_dir(const Gradient& g, int x, int y);
Point normalize(const Point& p);

#endif
