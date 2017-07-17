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


#include "include/undistort_rectilinear.h"
#include <complex>

typedef std::complex<double> cplex;

Undistort_rectilinear::Undistort_rectilinear(const cv::Rect& r, const vector<double>& in_coeffs) : Undistort(r) {
    coeffs = vector<double>(2, 0);
    for (size_t i=0; i < 2; i++) coeffs[i] = in_coeffs[i];
    
    radius_norm = sqrt(centre.x*centre.x + centre.y*centre.y);
}


static bool laguerre(vector<cplex>& a, cplex& x, int& its) {
    const int MR=8;
    const int MT=10;
    const int MAXIT=MT*MR;
    
    const double EPS=std::numeric_limits<double>::epsilon();
    static const double frac[MR+1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    
    int m = a.size() - 1;
    for (int iter=1; iter <= MAXIT; iter++) {
        its = iter;
        cplex b = a[m];
        double err = std::abs(b);
        cplex d(0,0);
        cplex f(0,0);
        double abx = std::abs(x);
        for (int j=m-1; j >= 0; j--) {
            f = x*f + d;
            d = x*d + b;
            b = x*b + a[j];
            err = std::abs(b) + abx*err;
        }
        err *= EPS;
        if (std::abs(b) <= err) return true; // on the root
        cplex g = d/b;
        cplex g2 = g*g;
        cplex h = g2 - 2.0*f/b;
        cplex sq = std::sqrt(double(m-1)*(double(m)*h-g2));
        cplex gp = g + sq;
        cplex gm = g - sq;
        double abp = std::abs(gp);
        double abm = std::abs(gm);
        gp = (abp < abm) ? gm : gp;
        cplex dx = std::max(abp, abm) > 0.0 ? double(m)/gp : std::polar(1+abx, double(iter));
        cplex x1 = x - dx;
        if (x == x1) return true;
        if (iter % MT != 0) {
            x = x1;
        } else {
            x -= frac[iter/MT]*dx;
        }
    }
    // no convergence ...
    return false;
}

static void lroots(vector<cplex>& a, vector<cplex>& roots) {
    const double EPS = 1e-14;
    int its;
    
    int m = a.size() - 1;
    vector<cplex> ad(a);
    for (int j=m-1; j >= 0; j--) {
        cplex x(0,0);
        vector<cplex> ad_v(j+2);
        for (int jj=0; jj < j+2; jj++) {
            ad_v[jj] = ad[jj];
        }
        laguerre(ad_v, x, its);
        if (fabs(x.imag()) <= 2.0*EPS*fabs(x.real())) {
            x = cplex(x.real(), 0.0);
        }
        roots[j] = x;
        cplex b = ad[j+1];
        for (int jj=j; jj >= 0; jj--) { // deflate the poly
            cplex c = ad[jj];
            ad[jj] = b;
            b = x*b + c;
        }
    }
}

static double laguerre_smallest_positive_root(double ru, double k1, double k2) {
    vector<cplex> a = {ru, -1.0, ru*k1, 0.0, ru*k2};
    vector<cplex> roots(4);
    lroots(a, roots);
    
    double minroot = 1e30;
    for (size_t i=0; i < roots.size(); i++) {
        if (roots[i].imag() == 0 && roots[i].real() >= 0 && roots[i].real() < minroot) {
            minroot = roots[i].real();
        }
    }
    
    // try to refine the root
    cplex root(minroot, 0.0);
    int its;
    laguerre(a, root, its);
    
    return root.real();
}

static bool lagsolve(double ru, double k1, double k2, double& root) {
    // we are looking for the roots of
    // P(r) = ru*k2*root^4 + ru*k1*root^2 - root + ru = 0
    // we would prefer the smallest positive root
    
    if (fabs(k1) < 1e-8 && fabs(k2) < 1e-8) {
        root = ru;
    } else {
        if (fabs(k2) < 1e-8) { // we only have a quadratic
            // a == 1
            double b = -1.0/(k1*ru);
            double c = 1.0/k1;
            double q = -0.5*(b + std::copysign(sqrt(b*b - 4*c), b));
            // force negative roots to become very large
            double r1 = q;
            double r2 = c/q;
            if (r1 < 0) {
                root = r2;
            } else {
                if (r2 < 0) {
                    root = r1;
                } else {
                    root = std::min(r1, r2);
                }
            }
        } else {
            root = laguerre_smallest_positive_root(ru, k1, k2);    
        }
    }
    return true;
}

cv::Point2d Undistort_rectilinear::slow_transform_point(double col, double row) {
    double px = (col + offset.x - centre.x);
    double py = (row + offset.y - centre.y);
    
    double ru = sqrt((px)*(px) + (py)*(py)) / radius_norm; 
    
    if (ru == 0) {
        px = col;
        py = row;
    } else {
        double rd=0;
        bool s = lagsolve(ru, coeffs[0], coeffs[1], rd); // not the fastest method, but it should be robust
        
        if (s) {
          px = px*rd/ru + centre.x - offset.x;
          py = py*rd/ru + centre.y - offset.y;
        } else {
          px = col;
          py = row;
        }
    }
    
    return cv::Point2d(px, py);
}    

cv::Point2d Undistort_rectilinear::inverse_transform_point(double col, double row) { // explicit inverse model
    double px = (col + offset.x - centre.x);
    double py = (row + offset.y - centre.y);

    double rd = sqrt((px)*(px) + (py)*(py)) / radius_norm; 
    double r2 = rd*rd;
    double ru = 1 + (coeffs[0] + coeffs[1]*r2)*r2;
    
    px = px/ru + centre.x - offset.x;
    py = py/ru + centre.y - offset.y;
    
    return cv::Point2d(px, py);
}

// note: second parameter is the raw Bayer image, which must also be padded out
cv::Mat Undistort_rectilinear::unmap(const cv::Mat& in_src, cv::Mat& rawimg) {
    
    const double buffer = 0.025;
    Point2d pi = inverse_transform_point(centre.x - (max_val.x+buffer*std::max(in_src.cols, in_src.rows)), centre.y);
    int pad_left = pi.x < 0 ? ceil(-pi.x) : 0;
    pi = inverse_transform_point(centre.x, centre.y - (max_val.y+buffer*std::max(in_src.cols, in_src.rows)));
    int pad_top = pi.y < 0 ? ceil(-pi.y) : 0;
    
    #if 0
    // TODO: we can add some box blur here to prevent the target squares from 
    // developing ragged edges
    // ideally, we should use the magnification norm(map_x(r,c) - c, map_y(r,c) - r) to
    // adjust the size of the box blur locally. It might be worthwhile to construct an
    // integral image to speed this up a bit
    cv::Mat bimg;
    cv::blur(in_src, bimg, cv::Size(5,5));
    cv::remap(bimg, timg, map_x, map_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    #endif
    
    return unmap_base(in_src, rawimg, pad_left, pad_top);
}
    
    
