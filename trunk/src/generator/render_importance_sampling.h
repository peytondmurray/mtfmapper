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
#include "airy_sampler.h"
#include "render.h"
#include "render_poly.h"

#include <algorithm>
using std::swap;

#ifndef SQR
  #define SQR(x) ((x)*(x))
#endif


//==============================================================================
class Render_rectangle_is : public Render_rectangle {
  protected:
    // we do not want this constructor to be called directly, but rather through
    // the named constructor idiom (below)
    Render_rectangle_is(double cx, double cy, double width, double height, double angle, 
        Render_type render_type=AIRY_PLUS_BOX,
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, int hs=60) 
        : Render_rectangle(cx, cy, width, height, angle, 1.0, 1.0, 0.0, false),
          aperture(in_aperture), pitch(in_pitch), lambda(in_lambda),
          poly(cx, cy, width, height, angle), sup(0, 0, 1, 1, 0), render_type(render_type), hs(hs) {

    }

  public:
    virtual ~Render_rectangle_is(void) {
    }
    
    double evaluate(double x, double y, double object_value, double background_value) const {
    
        double accumulator = 0;
        double wsum = 0;
        
        double sk = 0;
        double mk = 0;
        double prev_sk = 0;
        double prev_mk = 0;
        
        int samples_threshold = lrint(sqrt((double)nsamples)*12);
        
        // take initial batch of samples, checking for convergence along the way
        for (size_t sidx=0; sidx < size_t(samples_threshold); sidx++) {
        
            const double& weight = weights[sidx];
            double sample = sample_core(pos_x[sidx], pos_y[sidx], x, y, object_value, background_value);
            wsum += weight;
            
            accumulator += sample*weight;
            
            // calculate running mean and variance
            if (sidx == 0) {
                mk = prev_mk = sample*weight/wsum;
            } else {
                mk = prev_mk + weight/wsum*(sample - prev_mk);
                sk = prev_sk + weight*(sample - prev_mk)*(sample - mk);
                            
                prev_sk = sk; 
                prev_mk = mk;
            }
        } 
        
        const double thresh = SQR(0.5/(65536.0*65536.0));
        // if we have not converged yet, take the rest of the samples
        if (sk/wsum > thresh) {
            for (size_t sidx=samples_threshold; sidx < pos_x.size(); sidx++) {
                
                const double& weight = weights[sidx];
                double sample = sample_core(pos_x[sidx], pos_y[sidx], x, y, object_value, background_value);
                wsum += weight;
                
                accumulator += sample*weight;
            } 
        }
         
        double value = accumulator / (wsum);
        return value;
    }
    
    virtual string get_mtf_curve(void) const {
        return string("not defined");
    }
    
    virtual string get_psf_curve(void) const {
        return string("not defined");
    }
    
    virtual double get_mtf50_value(void) const {
        return 0;
    }
      
  protected:
    void initialise(void) {
    
        nsamples = SQR(hs*2 + 1);
        pos_x   = vector<double>(nsamples);
        pos_y   = vector<double>(nsamples);
        weights = vector<double>(nsamples);
        
        normal_sampler sampler;
        Airy_sampler airy(lambda, pitch, aperture);
        
        vector< pair<double, int> > radius_list;
        
        for (int sidx=0; sidx < nsamples; sidx++) {
            double& ex = pos_x[sidx];
            double& ey = pos_y[sidx];
            
            double sample_prob = airy.rairy2d(ex, ey, sampler);
            double rad = sqrt(ex*ex + ey*ey);
            
            rad /= (lambda/pitch)*aperture;
            
            // collect the indices of the outermost points
            // this ensures that the early stopping criterion works well
            if (rad > Airy_sampler::diam*0.15) {
                radius_list.push_back(make_pair<double, int>(fabs(rad - Airy_sampler::diam*0.8), sidx));
            }
            
            double jinc_weight = 4*Airy_sampler::jinc(rad)*Airy_sampler::jinc(rad);
            
            weights[sidx] = jinc_weight/sample_prob;
        } // supersamples
        
        // now swap out the outermost points with the first points
        sort(radius_list.begin(), radius_list.end());
        for (size_t i=0; i < radius_list.size(); i++) {
            swap(pos_x[i], pos_x[radius_list[i].second]);
            swap(pos_y[i], pos_y[radius_list[i].second]);
            swap(weights[i], weights[radius_list[i].second]);
        }
        
        printf("using IS renderer with %d samples per pixel\n", nsamples);
    }
  
    virtual double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const  = 0;
    
    double bisect_airy(double (*f)(double x, double s, double _p), double p = 0) const { // simple, but robust
        const int nmax = 100;
        const double tol = 1e-7;
        double s = (lambda/pitch) * aperture;
        double a = 1e-5;
        double b = min((1-1e-5) / s, 1.0);
        double fa = f(a, s, p);
        double fb = f(b, s, p);
        for (int n=0; n < nmax; n++) {
            double c = 0.5 * (a + b);
            double fc = f(c, s, p);
            if (fc == 0 || (b - a) < tol) {
                return c;
            }
            if (fc*fa >= 0) {
                a = c;
                fa = fc;
            } else {
                b = c;
                fb = fc;
            }
        }
        return a;
    }
    
    vector<double> weights;
    
    double aperture;
    double pitch;
    double lambda;
    
    Render_rectangle_poly poly;
    Render_rectangle_poly sup;
    
    int nsamples;
    
    Render_type render_type;
    int hs;
};


#include "render_is_airy.h"
#include "render_is_airybox.h"
#include "render_is_airyolpf.h"

Render_rectangle_is* build_psf(Render_rectangle::Render_type render_t, double cx, double cy, double width, double height, double angle,                  
    double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55, double olpf_split=0.375,
    int hs=0) {

    switch (render_t) {
        case Render_rectangle::AIRY: 
            return new Render_rectangle_is_airy(cx, cy, width, height, angle, 
                    in_aperture, in_pitch, in_lambda,
                    hs == 0 ? 60 : hs
                );
            break;
        case Render_rectangle::AIRY_PLUS_BOX:
            return new Render_rectangle_is_airybox(cx, cy, width, height, angle, 
                    in_aperture, in_pitch, in_lambda,
                    hs == 0 ? 40 : hs
                );
            break;
        case Render_rectangle::AIRY_PLUS_4DOT_OLPF:
            return new Render_rectangle_is_airyolpf(cx, cy, width, height, angle, 
                    in_aperture, in_pitch, in_lambda, 
                    olpf_split, hs == 0 ? 30 : hs
                );
            break;
        default:
            printf("Warning! Trying to use an unsupported PSF type with Render_rectangle_is path!\n");
            break;
    }
    return 0;
}

#endif // RENDER_H
