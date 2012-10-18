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
#ifndef RENDER_IMPORTANCE_SAMPLING_AIRYOLPF_H
#define RENDER_IMPORTANCE_SAMPLING_AIRYOLPF_H

#include "include/common_types.h"

#include "airy_sampler.h"
#include "render.h"
#include "render_poly.h"
#include "render_importance_sampling.h"

//==============================================================================
class Render_rectangle_is_airyolpf : public Render_rectangle_is {
  public:
    Render_rectangle_is_airyolpf(double cx, double cy, double width, double height, double angle, 
        double in_aperture=8, double in_pitch=4.73, double in_lambda=0.55) 
        : Render_rectangle_is(
            cx, cy, width, height, angle, 
            AIRY_PLUS_4DOT_OLPF, in_aperture, in_pitch, in_lambda,
            30 // half-samples
          )
          {
          
          initialise();
    }

    virtual ~Render_rectangle_is_airyolpf(void) {
    }
    
    virtual string get_mtf_curve(void) const {
        char buffer[1024];
        double scale = (lambda/pitch) * aperture;
        sprintf(buffer, "2.0/pi*abs( 0.875*sin(1.75*x*pi)/(1.75*x*pi) + 0.125*sin(0.25*x*pi)/(0.25*x*pi) )*(acos(x*%lg) - (x*%lg)*sqrt(1-(x*%lg)*(x*%lg)))", 
            scale, scale, scale, scale
        );
        return string(buffer);
    }
    
    virtual string get_psf_curve(void) const {
        return(string("not implemented (no simple analytical form)"));
    }
    
    virtual double get_mtf50_value(void) const {
        return bisect_airy(&airy_olpf_mtf);
    }
      
  protected:
    virtual inline double sample_core(const double& ex, const double& ey, const double& x, const double& y,
        const double& object_value, const double& background_value) const {
    
        double sample = 0;
        
        const double olpf_x[4] = {-0.375, -0.375, 0.375, 0.375};
        const double olpf_y[4] = {-0.375, 0.375, -0.375, 0.375};
        for (int k=0; k < 4; k++) {
            double area = poly.evaluate_x(sup, ex + x + olpf_x[k], ey + y + olpf_y[k]);
            sample += 0.25*(object_value * area + background_value * (1 - area));
        }
                
        return sample;
    }
    
    double static airy_olpf_mtf(double x, double s)  {
        return 2.0/M_PI * 
          fabs( 0.875*sin(1.75*x*M_PI)/(1.75*x*M_PI) + 0.125*sin(0.25*x*M_PI)/(0.25*x*M_PI) ) * 
          (acos(x*s) - (x*s)*sqrt(1-(x*s)*(x*s))) - 0.5;
    }
};

#endif // RENDER_H
