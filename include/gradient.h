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
#ifndef GRADIENT_H
#define GRADIENT_H

#include <cmath>
#include "include/common_types.h"

class Gradient {
public:
    Gradient(const cv::Mat& in_img, bool snapshot = false, bool thin = false);
    virtual ~Gradient(void);

    inline const float* grad_x(void) const {
        return _gradient_x;
    }

    inline const float* grad_y(void) const {
        return _gradient_y;
    }
    
    inline float grad_x(int x, int y) const {
        return _gradient_x[y*_width + x];
    }

    inline float grad_y(int x, int y) const {
        return _gradient_y[y*_width + x];
    }

    inline const float* grad_magnitude(void) const {
        return _gradient_m;
    }
    
    inline float grad_magnitude(int x, int y) const {
        return _gradient_m[y*_width + x];
    }

    inline int width(void) const {
        return _width;
    }

    inline int height(void) const {
        return _height;
    }

    static void _blur(const cv::Mat& image, float **smoothedim);
private:

    void _compute_gradients(const float* smoothed_im, int rows, int cols, float** grad_x, float** grad_y);

protected:
    int _width;
    int _height;

    float*      _gradient_x;
    float*      _gradient_y;
    float*      _gradient_m;
    float*      _smoothed;

    bool    _thin;
};

#endif // GRADIENT_H


