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


#include "include/undistort.h"

// find closest point in radmap lookup table, then apply
// quadratic interpolation to refine the value
cv::Point2d Undistort::transform_point(double c, double r) {
    double px = (c + offset.x - centre.x);
    double py = (r + offset.y - centre.y);
    
    double rad = 2*sqrt(px*px + py*py); // convert to half-pixel spacing to match how radmap was built
    int rad_f = (int)rad;
    
    if (rad_f > (int)radmap.size() - 3) {
        if (radmap.size() == 0) { // this should only happen if we call transform_point before unmap, which should not happen
            logger.error("Warning: radmap not initialized in transform_point. Calling build_radmap, but fix the code!\n");
            build_radmap();
        } else {
            logger.error("Warning: radmap range exceeded in transform_point. Clamping (rad-%lf, rad_f=%d)\n", rad, rad_f);
            rad_f = (int)radmap.size() - 3;
        }
    }
    
    double w = rad - rad_f;
    // quadratic interpolation through half-pixel spaced points
    double fc = radmap[rad_f];
    double fb = -1.5*radmap[rad_f] + 2*radmap[rad_f+1] - 0.5*radmap[rad_f+2];
    double fa =  0.5*radmap[rad_f] - radmap[rad_f+1] + 0.5*radmap[rad_f+2];
    double rad_d = 2*(fa*w*w + fb*w + fc);
    
    // avoid divide-by-zero
    if (rad < 1e-8) {
        return centre;
    }
    return Point2d(px, py) * (rad_d/rad) + centre - offset;
}

void Undistort::build_radmap(void) {
    radmap.clear();
    double maxrad = 1.1*sqrt(centre.x*centre.x + centre.y*centre.y); // TODO: add better support for non-centered CoD
    for (int rad=0; rad <= 2*(maxrad + 2); rad++) { // half-pixel spacing
        Point2d tp = slow_transform_point(centre.x + 0.5*rad, centre.y);
        double rd = norm(tp - Point2d(centre.x, centre.y));
        radmap.push_back(rd);
    }
}