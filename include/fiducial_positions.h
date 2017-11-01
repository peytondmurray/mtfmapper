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
#ifndef FIDUCIAL_POSITIONS_H
#define FIDUCIAL_POSITIONS_H

#include "common_types.h"

class Fiducial {
  public:
    Fiducial(double img_x=0, double img_y=0, double real_x=0, double real_y=0, int code=0, int quadrant=0) 
    : icoords(img_x, img_y), rcoords(real_x, real_y), code(code), quadrant(quadrant) {}

    Point2d icoords;
    Point2d rcoords;
    int code;
    int quadrant;
};

const int n_fiducials = 22;

const Fiducial main_fiducials[n_fiducials] = {
    {0, 0, 30.000000, 0.000000, 0, 0},
    {0, 0, -30.000000, 0.000000, 0, 0},
    {0, 0, 52.130517, 36.502181, 1, 2},
    {0, 0, 36.502181, -52.130517, 2, 3},
    {0, 0, -52.130517, -36.502181, 3, 0},
    {0, 0, -36.502181, 52.130517, 4, 1}, 
    {0, 0, 81.915204, -57.357644, 5, 3},
    {0, 0, -58.791585, -83.963085, 6, 0},
    {0, 0, -86.010965, 60.225526, 7, 1}, 
    {0, 0, 61.659467, 88.058845, 8, 2},
    {0, 0, 112.763114, 41.042417, 9, 2},
    {0, 0, 41.897468, -115.112346, 10, 3},
    {0, 0, -117.461578, -42.752518, 11, 0},
    {0, 0, -43.607568, 119.810809, 12, 1}, 
    {0, 0, 44.462619, 122.160041, 13, 2},
    {0, 0, 124.509272, -45.317669, 14, 3},
    {0, 0, -46.172719, -126.858504, 15, 0},
    {0, 0, -129.207735, 47.027770, 16, 1}, 
    {0, 0, 98.994949, 98.994949, 17, 2},
    {0, 0, 100.762716, -100.762716, 18, 3},
    {0, 0, -102.530483, -102.530483, 19, 0},
    {0, 0, -104.298250, 104.298250, 20, 1}
};

#endif
