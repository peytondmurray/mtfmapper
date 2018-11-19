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
#ifndef RENDER_H
#define RENDER_H

#include "include/common_types.h"

#include <cv.h>
const double transparent = -1.0;

#include "normal_sampler.h"

//==============================================================================
class Render_target {
  public:
      Render_target(void) {}
      virtual double evaluate(int x, int y, double oject_value, double background_value) const = 0;
};


//==============================================================================
class Render_rectangle : public Render_target {
  public:
    Render_rectangle(double tlx, double tly, double width, double height, double angle, double in_sigma=6.0) : sigma(in_sigma) {
        bases[0] = cv::Vec2d(tlx, tly);
        bases[1] = cv::Vec2d(tlx + cos(angle)*width, tly - sin(angle)*width);
        bases[2] = cv::Vec2d(tlx + cos(angle)*width + sin(angle)*height, tly - sin(angle)*width + cos(angle)*height);
        bases[3] = cv::Vec2d(tlx + sin(angle)*height, tly + cos(angle)*height);
              
        normals[0] = (bases[2] - bases[1]); 
        normals[1] = (bases[3] - bases[2]); 
        normals[2] = (bases[0] - bases[3]); 
        normals[3] = (bases[1] - bases[0]); 
        for (size_t i=0; i < 4; i++) {
            double n = norm(normals[i]);
            normals[i] = normals[i]*(1.0/n);
        }
              
        hs = 22; // seems like enough samples for up to sigma=6, at least
              
        int nsamples = hs*2 + 1;
        pos_x   = cv::Mat_<double>(nsamples, nsamples);
        pos_y   = cv::Mat_<double>(nsamples, nsamples);
              
        normal_sampler sampler;
        for (int ss_x=-hs; ss_x <= hs; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs; ss_y++) {
                  
                double ex = 0;
                double ey = 0;
                    
                sampler.rnorm2d(ex, ey, sigma);
                    
                pos_x(ss_y+hs, ss_x+hs) = ex;
                pos_y(ss_y+hs, ss_x+hs) = ey;
            }
        } // supersamples
    }
    
    double evaluate(int x, int y, double object_value, double background_value) const {
   
        double accumulator = 0;
        for (int ss_x=-hs; ss_x <= hs; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs; ss_y++) {
            
                double ex = pos_x(ss_y+hs, ss_x+hs);
                double ey = pos_y(ss_y+hs, ss_x+hs);
                
                if ( is_inside(ex + x, ey + y) ) {
                    accumulator += object_value;
                } else {
                    accumulator += background_value;
                }
            }
        } // supersamples
         
        double value = accumulator / ((2*hs+1)*(2*hs+1));
        return value;
    }
      
  private:
    inline bool is_inside(double x, double y) const {
        bool inside = true;
        for (int i=0; i < 4 && inside; i++) {
            cv::Vec2d dir(x - bases[i][0], y - bases[i][1]);
            double dot = dir.dot(normals[i]);
            if (dot < 0) {
                inside = false;
            }
        }
        return inside;
    }      
      
    double sigma;
    int    hs;
    cv::Vec2d normals[4];
    cv::Vec2d bases[4];
    cv::Mat_<double> pos_x;
    cv::Mat_<double> pos_y;
};

#endif // RENDER_H
