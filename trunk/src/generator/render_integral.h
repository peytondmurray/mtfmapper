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
#ifndef RENDER_INTEGRAL_H
#define RENDER_INTEGRAL_H

#include "include/common_types.h"

#include <cv.h>
#include <math.h>
#include "render.h"


//==============================================================================
class Render_rectangle_integral : public Render_rectangle {
  public:
    Render_rectangle_integral(double cx, double cy, double width, double height, double angle, double in_sigma=6.0) : 
        Render_rectangle(cx, cy, width, height, angle, in_sigma) {
        
        for (int k=0; k < 4; k++) {
            yvals.push_back(bases[k][1]);
        }
        sort(yvals.begin(), yvals.end());
        
        left_base = vector<cv::Vec2d>(3);
        left_dir  = vector<cv::Vec2d>(3);
        right_base = vector<cv::Vec2d>(3);
        right_dir  = vector<cv::Vec2d>(3);
        
        for (int k=0; k < 3; k++) {
            find_intersection_vectors((yvals[k]+yvals[k+1])/2.0, 
                left_base[k], left_dir[k], 
                right_base[k], right_dir[k]
            );
        }
    }
    
    double evaluate(int x, int y, double object_value, double background_value) const {
   
        double accumulator = 0;
        
        const double eps = 1e-7;
        if ( fabs(yvals[0] - yvals[1]) > 1e-6 ) {
            
            if (yvals[0] <= y && y <=yvals[1]) { //
                accumulator += adaptive_simpson(yvals[0], y, eps, 30, x, y);
                accumulator += adaptive_simpson(y, yvals[1], eps, 30, x, y);
            } else {
                accumulator += adaptive_simpson(yvals[0], yvals[1], eps, 30, x, y);
            }
        }
        if ( fabs(yvals[1] - yvals[2]) > 1e-6 ) {
            if (yvals[1] <= y && y <=yvals[2]) { //
                accumulator += adaptive_simpson(yvals[1], y, eps, 30, x, y);
                accumulator += adaptive_simpson(y, yvals[2], eps, 30, x, y);
            } else {
                accumulator += adaptive_simpson(yvals[1], yvals[2], eps, 30, x, y);
            }
        }
        if ( fabs(yvals[2] - yvals[3]) > 1e-6 ) {
            if (yvals[2] <= y && y <=yvals[3]) { //
                accumulator += adaptive_simpson(yvals[2], y, eps, 30, x, y);
                accumulator += adaptive_simpson(y, yvals[3], eps, 30, x, y);
            } else {
                accumulator += adaptive_simpson(yvals[2], yvals[3], eps, 30, x, y);
            }
        }
        accumulator /= sqrt(2.0);
        
        return background_value*(1 - accumulator) + accumulator*object_value;
    }
      
  private:
  
    vector<cv::Vec2d> left_base;
    vector<cv::Vec2d> right_base;
    vector<cv::Vec2d> left_dir;
    vector<cv::Vec2d> right_dir;
    vector<double>    yvals;
  
    inline double gaussian_integral(double x) const {
        if (x >= 0) return (1 + erf(x))*0.5;
        if (x < 0) return (1 - erf(-x))*0.5;
        return 0;
    }
    
    inline void find_intersections(double y, double& left, double& right) const {
        int yidx = 0;
        if (y > yvals[1]) yidx++;
        if (y > yvals[2]) yidx++;
        
        double lt = (y - left_base[yidx][1])/left_dir[yidx][1];
        left = (lt*left_dir[yidx] + left_base[yidx])[0];
        
        double rt = (y - right_base[yidx][1])/right_dir[yidx][1];
        right = (rt*right_dir[yidx] + right_base[yidx])[0];
    }
    
    inline void find_intersection_vectors(double y, 
        cv::Vec2d& left_base, cv::Vec2d& left_dir,
        cv::Vec2d& right_base, cv::Vec2d& right_dir
        ) const {
        double left = 1e50;
        double right = -1e50;
        for (int k=0; k < 4; k++) {
            int k1 = (k + 1) % 4;
            cv::Vec2d dir = bases[k] - bases[k1];
            double ydelta = y - bases[k1][1];
            double t = ydelta / (dir[1]);
            if (t >= 0 && t <= 1) {
                cv::Vec2d p = dir*t + bases[k1];
                if (p[0] < left) { left = p[0];  left_base = bases[k1];  left_dir = dir; }
                if (p[0] > right){ right = p[0]; right_base = bases[k1]; right_dir = dir; }
            }
        }
    }
    
    inline double func(double y, double px, double py) const {
        double fy = exp(-(y - py)*(y - py)/(2*sigma*sigma)) / (sigma*sqrt(M_PI));
        
        double left;
        double right;
        
        find_intersections(y, left, right);
        right = (right - px) / (sqrt(2.0)*sigma);
        left  = (left - px) / (sqrt(2.0)*sigma);
        
        return fy * (gaussian_integral(right) - gaussian_integral(left));
    }
    
    double adaptive_simpson_aux(double a, double b, 
        double epsilon, double S, double fa, 
        double fb, double fc, int bottom, 
        double px, double py, int depth) const {
        
        double c = (a+b)/2;
        double h = b - a;
        double d = (a+c)/2;
        double e = (c+b)/2;
        double fd = func(d, px, py);
        double fe = func(e, px, py);
        double Sleft = (h/12)*(fa + 4*fd + fc);
        double Sright = (h/12)*(fc + 4*fe + fb);
        double S2 = Sleft + Sright;
        
        if (bottom <= 0 || (depth >= 4 && fabs(S2 - S) <= 15*epsilon) ) {
            return S2 + (S2 - S)/15;
        }
        return adaptive_simpson_aux(a, c, epsilon/2, Sleft, fa, fc, fd, bottom - 1, px, py, depth + 1) +
               adaptive_simpson_aux(c, b, epsilon/2, Sright, fc, fb, fe, bottom - 1, px, py, depth + 1);
    }
    
    double adaptive_simpson(double a, double b,
        double epsilon, int max_depth,
        double px, double py) const {
        
        double c = (a+b)/2;
        double h = b - a;
        double fa = func(a, px, py);
        double fb = func(b, px, py);
        double fc = func(c, px, py);
        double S = (h/6)*(fa + 4*fc + fb);
        
        return adaptive_simpson_aux(a, b, epsilon, S, fa, fb, fc, max_depth, px, py, 1);
    }          
      
};

#endif // RENDER_H