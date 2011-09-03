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
      virtual double evaluate(int x, int y, double background) const = 0;
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
        weights = cv::Mat_<double>(nsamples, nsamples);
        pos_x   = cv::Mat_<double>(nsamples, nsamples);
        pos_y   = cv::Mat_<double>(nsamples, nsamples);
              
        normal_sampler sampler;
        weightsum = 0;
        for (int ss_x=-hs; ss_x <= hs; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs; ss_y++) {
                  
                double ex = 0;
                double ey = 0;
                    
                sampler.rnorm2d(ex, ey, 2*sigma);
                    
                double weight = exp(-(ex*ex + ey*ey)/(2*sigma*sigma)) / (2*M_PI*sigma*sigma);
                weights(ss_y+hs, ss_x+hs) = weight; 
                pos_x(ss_y+hs, ss_x+hs) = ex;
                pos_y(ss_y+hs, ss_x+hs) = ey;
                weightsum += weight;
            }
        } // supersamples
    }
    
    double evaluate(int x, int y, double background) const {
   
        double accumulator = 0;
        for (int ss_x=-hs; ss_x <= hs; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs; ss_y++) {
            
                double ex = pos_x(ss_y+hs, ss_x+hs);
                double ey = pos_y(ss_y+hs, ss_x+hs);
                double weight = weights(ss_y+hs, ss_x+hs);
                
                if ( is_inside(ex + x, ey + y) ) {
                    accumulator += background * weight;
                } else {
                    accumulator += (1.0 - background) * weight;
                }
            }
        } // supersamples
         
        double value = accumulator / weightsum;
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
    double weightsum;
    int    hs;
    cv::Vec2d normals[4];
    cv::Vec2d bases[4];
    cv::Mat_<double> weights;
    cv::Mat_<double> pos_x;
    cv::Mat_<double> pos_y;
};

#endif // RENDER_H
