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
#ifndef RENDER_IMPORTANCE_SAMPLING_H
#define RENDER_IMPORTANCE_SAMPLING_H

#include "include/common_types.h"

#include <cv.h>
#include "normal_sampler.h"
#include "render.h"

#include "render_poly.h"

#include <algorithm>

#ifndef SQR
  #define SQR(x) ((x)*(x))
#endif

class Airy_sampler {
  public:
    Airy_sampler(double lambda, double pitch, double N) 
      : lambda(lambda), pitch(pitch), N(N),
        nsamples(7*100000), x(nsamples, 0), y(nsamples, 0) {
        
        // first point is replicated to allow interpolation later
        x[0] = 0;
        y[0] = 0;
        step = 7/double(nsamples);
        for (size_t n=1; n < nsamples; n++) {
            x[n] = sqrt(((n-1)/double(nsamples)))*7;
            y[n] = y[n-1] + 4*jinc(x[n])*jinc(x[n]) * step;
            
        }
        for (size_t n=0; n < nsamples; n++) {
            y[n] /= y[nsamples-1];
        }
    }
    
    double rairy2d(double& sx, double& sy, normal_sampler& sampler) {
        double px;
        double py;
        sampler.runif2d(px, py, 0.5, 0.5);
        px += 0.5; // undo the scaling of runif2d
        py += 0.5;
        size_t xidx = upper_bound(y.begin(), y.end(), px) - y.begin();
        double xp;
        // interpolate tp prevent stratification
        if (xidx < nsamples - 1) {
            double alpha = (px - y[xidx])/(y[xidx+1] - y[xidx]);
            xp = x[xidx]*(1-alpha) + x[xidx+1]*alpha;
        } else {
            xp = x[xidx];
        }
        double rad = (xp) * (lambda/pitch)*N;
        sx = rad*cos(py*2*M_PI);
        sy = rad*sin(py*2*M_PI);
        if (xidx == 0) return 0;
        return (y[xidx] - y[xidx-1])/step;
    }
    
    inline double jinc(double x) {
        if (fabs(x) == 0) {
            return 0.5;
        }
        return j1(x*M_PI)/(x*M_PI);
    }
    
    double lambda;
    double pitch;
    double N;
    
    size_t nsamples;
    vector<double> x;
    vector<double> y;
    double step;
};

//==============================================================================
class Render_rectangle_is : public Render_rectangle {
  public:
    Render_rectangle_is(double cx, double cy, double width, double height, double angle, 
        double in_sigma=6.0, double minor_sigma=6.0, double theta=0) 
        : Render_rectangle(cx, cy, width, height, angle, in_sigma, minor_sigma, theta),
          poly(cx, cy, width, height, angle, in_sigma)  {
        
        hs = 8; // seems like enough samples for up to sigma=6, at least
        
        int nsamples = hs*2 + 1;
        pos_x   = cv::Mat_<double>(nsamples, nsamples);
        pos_y   = cv::Mat_<double>(nsamples, nsamples);
        weights = cv::Mat_<double>(nsamples, nsamples);
        
        aperture = 2.8;
        lambda = 0.55;
        pitch = 4.73;
        
        
        normal_sampler sampler;
        Airy_sampler airy(lambda, pitch, aperture);
        
        for (int ss_x=-hs; ss_x <= hs; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs; ss_y++) {
                  
                double ex = 5*ss_x+hs/double(hs);
                double ey = 5*ss_y+hs/double(hs);
                
                double sample_prob = 1.0;
                
                sample_prob = airy.rairy2d(ex, ey, sampler);
                //sampler.runif2d(ex, ey, 8, 8);
                double rad = sqrt(ex*ex + ey*ey);
                
                rad /= (lambda/pitch)*aperture;
                double sinc_weight = 4*jinc(rad)*jinc(rad);
                
                pos_x(ss_y+hs, ss_x+hs) = ex;
                pos_y(ss_y+hs, ss_x+hs) = ey;
                weights(ss_y+hs, ss_x+hs) = sinc_weight/sample_prob;
                //fprintf(stderr, "%le %le %le %le\n", ex, ey, weights(ss_y+hs, ss_x+hs), sinc_weight);
                //fprintf(stderr, "%le %le %le %le\n", ex, ey, gauss_prob, sinc_weight);
            }
        } // supersamples
        printf("using IS renderer : %d\n", nsamples*nsamples);
    }
    
    double evaluate(int x, int y, double object_value, double background_value) const {
    
        double accumulator = 0;
        double wsum = 0;
        int sindex = 0;
        
        double sk = 0;
        double mk = 0;
        double prev_sk = 0;
        double prev_mk = 0;
        
        Render_rectangle_poly sup(0, 0, 1, 1, 0);
        
        vector<double> sampled;
        bool done = false;
        for (int ss_x=-hs; ss_x <= hs && !done; ss_x++) {
            for (int ss_y=-hs; ss_y <= hs && !done; ss_y++) {
            
                double ex = pos_x(ss_y+hs, ss_x+hs);
                double ey = pos_y(ss_y+hs, ss_x+hs);
                
                double weight = weights(ss_y + hs, ss_x + hs);
                double sample = 0;
                #if 1
                  #if 0
                    double area = poly.evaluate_x(sup, ex + x, ey + y);
                    sample = object_value * area + 
                             background_value * (1 - area);
                    wsum += weight;
                  #else
                    const double olpf_x[4] = {-0.375, -0.375, 0.375, 0.375};
                    const double olpf_y[4] = {-0.375, 0.375, -0.375, 0.375};
                    for (int k=0; k < 4; k++) {
                        double area = poly.evaluate_x(sup, ex + x + olpf_x[k], ey + y + olpf_y[k]);
                        sample += 0.25*(object_value * area + background_value * (1 - area));
                    }
                    wsum += weight;
                  #endif
                #else
                if ( is_inside(ex + x, ey + y) ) {
                    sample = object_value;
                } else {
                    sample = background_value;
                }
                wsum += weight;
                #endif
                
                accumulator += sample*weight;
                if (sindex == 0) {
                    mk = prev_mk = sample*weight/wsum;
                } else {
                    mk = prev_mk + weight/wsum*(sample - prev_mk);
                    sk = prev_sk + weight*(sample - prev_mk)*(sample - mk);
                                
                    prev_sk = sk; 
                    prev_mk = mk;
                }
                double sdev = sqrt(sk/wsum);
                
                if (x == 72 && y == 70) { 
                    sampled.push_back(accumulator/wsum);
                    
                    if (sampled.size() % 1000 == 0 || (2*hs+1)*(2*hs+1) - sampled.size() < 10) {
                        double smean = 0;
                        for (size_t n=0; n < sampled.size(); n++) {
                            smean += sampled[n];
                        }
                        smean /= double(sampled.size());
                        double svar = 0;
                        for (size_t n=0; n < sampled.size(); n++) {
                            svar += (sampled[n] - smean)*(sampled[n] - smean);
                        }
                        svar /= double(sampled.size());
                        printf("%d %le %le (%le) : (%le / %le)\n", sindex, mk, accumulator/wsum, sdev, smean, sqrt(svar));
                    }
                }
                
                const double thresh = 0.5/(65536.0*65536.0);
                if (sindex > 1000 && sdev < thresh) {
                    done = true;
                }
                
                sindex++;
            }
        } // supersamples
         
        double value = accumulator / (wsum);
        return value;
    }
      
  protected:
    inline double sinc(double x) {
        if (fabs(x) < 1e-6) {
            return 1.0;
        }
        return sin(x*M_PI)/(x*M_PI);
    }
    
    inline double jinc(double x) {
        if (fabs(x) == 0) {
            return 0.5;
        }
        return j1(x*M_PI)/(x*M_PI);
    }
  
    cv::Mat_<double> weights;
    
    double aperture;
    double pitch;
    double lambda;
    
    Render_rectangle_poly poly;
    vector<Render_rectangle_poly> supports;
};

#endif // RENDER_H
