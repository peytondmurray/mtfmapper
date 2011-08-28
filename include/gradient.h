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


