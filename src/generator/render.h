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
      virtual double evaluate(double x, double y, double oject_value, double background_value) const = 0;
};


//==============================================================================
class Render_rectangle : public Render_target {
  public:
    typedef enum {
        GAUSSIAN,
        GAUSSIAN_SAMPLED,
        AIRY,
        AIRY_PLUS_BOX,
        AIRY_PLUS_4DOT_OLPF
    } Render_type;
  
    Render_rectangle(double cx, double cy, double width, double height, double angle, 
        double in_sigma=6.0, double minor_sigma=6.0, double theta=0, bool init=true) : sigma(in_sigma), cx(cx), cy(cy) {
        bases[0] = cv::Vec2d(width/2, height/2);
        bases[1] = cv::Vec2d(-width/2, height/2);
        bases[2] = cv::Vec2d(-width/2, -height/2);
        bases[3] = cv::Vec2d(width/2, -height/2);
        
        // rotate corners
        for (size_t i=0; i < 4; i++) {
            double ox = bases[i][0];
            double oy = bases[i][1];
            bases[i][0] = cos(angle)*ox - sin(angle)*oy + cx;
            bases[i][1] = sin(angle)*ox + cos(angle)*oy + cy;
        }
              
        normals[0] = (bases[2] - bases[1]); 
        normals[1] = (bases[3] - bases[2]); 
        normals[2] = (bases[0] - bases[3]); 
        normals[3] = (bases[1] - bases[0]); 
        for (size_t i=0; i < 4; i++) {
            double n = norm(normals[i]);
            normals[i] = normals[i]*(1.0/n);
        }
        
        if (init) {
              
            hs = 22; // seems like enough samples for up to sigma=6, at least
            int nsamples = SQR(hs*2 + 1);
            pos_x   = vector<double>(nsamples);
            pos_y   = vector<double>(nsamples);
        
            normal_sampler sampler;
            for (int sidx=0; sidx < nsamples; sidx++) {
                sampler.rnorm2d(pos_x[sidx], pos_y[sidx], sigma, minor_sigma, theta);
            } 
        }
    }
    
    virtual ~Render_rectangle(void) {
    }
    
    double evaluate(double x, double y, double object_value, double background_value) const {
   
        double accumulator = 0;
        for (size_t sidx=0; sidx < pos_x.size(); sidx++) {
                if ( is_inside(pos_x[sidx] + x, pos_y[sidx] + y) ) {
                    accumulator += object_value;
                } else {
                    accumulator += background_value;
                }
        } // supersamples
         
        double value = accumulator / ((2*hs+1)*(2*hs+1));
        return value;
    }
      
  public:
    inline bool is_inside(double x, double y) const {
        bool inside = true;
        for (int i=0; i < 4 && inside; i++) {
            double dot = normals[i][0]*(x - bases[i][0]) + 
                         normals[i][1]*(y - bases[i][1]);
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
    vector<double> pos_x;
    vector<double> pos_y;
    double cx;
    double cy;
};

#endif // RENDER_H
