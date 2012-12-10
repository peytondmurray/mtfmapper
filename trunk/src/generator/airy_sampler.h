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
#ifndef AIRY_SAMPLER_H
#define AIRY_SAMPLER_H

#include "include/common_types.h"

class Airy_sampler {
  public:
    Airy_sampler(double lambda, double pitch, double N) 
      : lambda(lambda), pitch(pitch), N(N),
        nsamples(3*7000), x(nsamples, 0), y(nsamples, 0) {
        
        // first point is replicated to allow interpolation later
        x[0] = 0;
        y[0] = 0;
        for (size_t n=1; n < nsamples; n++) {
            x[n] = ((n-1)/double(nsamples))*diam;
            y[n] = y[n-1] + 4*jinc(x[n])*jinc(x[n]) * ((M_PI*x[n]*x[n] - M_PI*x[n-1]*x[n-1])); // correct PRASA line
            // We can bias the edges a little to ensure more samples are taken there, but these samples
            // have very low weights, thus contribute vary little. Potentially they may help to ensure the early 
            // convergence tests function better ...
            //y[n] = y[n-1] + (4*jinc(x[n])*jinc(x[n])*(1-1e-5) + (1e-5)) * ((M_PI*x[n]*x[n] - M_PI*x[n-1]*x[n-1]));  
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
        py *= 2*M_PI;
        size_t xidx = upper_bound(y.begin(), y.end(), px) - y.begin();
        double xp;
        // interpolate to prevent stratification
        if (xidx < nsamples - 1) {
            double alpha = (px - y[xidx])/(y[xidx+1] - y[xidx]);
            xp = x[xidx]*(1-alpha) + x[xidx+1]*alpha;
        } else {
            xp = x[xidx];
        }
        double rad = (xp) * (lambda/pitch)*N;
        sx = rad*cos(py);
        sy = rad*sin(py);
        if (xidx == 0) return 1e-6;
        return (y[xidx] - y[xidx-1]) / ((M_PI*x[xidx]*x[xidx] - M_PI*x[xidx-1]*x[xidx-1])); 
    }
    
    inline static double jinc(double x) {
        if (fabs(x) == 0) {
            return 0.5;
        }
        return _j1(x*M_PI)/(x*M_PI);
    }
    
    static const double diam;
    
    double lambda;
    double pitch;
    double N;
    
    size_t nsamples;
    vector<double> x;
    vector<double> y;
};

const double Airy_sampler::diam = 45.0;

#endif

