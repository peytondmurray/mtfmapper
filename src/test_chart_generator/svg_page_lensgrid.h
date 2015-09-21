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
#ifndef SVG_PAGE_LENSGRID_H
#define SVG_PAGE_LENSGRID_H

#include "svg_page.h"
#include "include/mtf50_edge_quality_rating.h"

class Svg_page_lensgrid : public Svg_page {
  public:
    Svg_page_lensgrid(const string& page_spec, const string& fname) 
    : Svg_page(page_spec, fname), centre(0.5, 0.5*sqrt(2.0)),
      off(0), scale(10.0/100) {
    
    }
    
    
    
    void render(void) {
        grid(0.04, 17, 11); 
    }
    
  protected:
    virtual iPoint project(double x, double y) {
        
        x = floor(x*width);
        y = floor(y*width);
        
        return iPoint(int(x), int(y));
    }

    double quant_angle(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    void place_trapezoid(double xpos, double ypos, double width) {
        dPoint delta = dPoint(xpos, ypos) - centre;
        
        //#define SQUARES
        
        //rotated_square(xpos, ypos, width*0.5, eang);
        double norm = sqrt(delta.x*delta.x + delta.y*delta.y);
        delta.x /= norm;
        delta.y /= norm;
        
        dPoint omx = dPoint(xpos + delta.x*width*sqrt(0.5), ypos + delta.y*width*sqrt(0.5));
        dPoint tang = dPoint(-delta.y, delta.x);
        
        dPoint tl(omx.x - 0.5*width*tang.x, omx.y - 0.5*width*tang.y);
        iPoint p = project(tl.x, tl.y);
        fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
        
        fprintf(stderr, "4\n%lf %lf\n", p.y*scale + off, p.x*scale + off);
        
        dPoint rad = tl - centre;
        norm = sqrt(rad.x*rad.x + rad.y*rad.y);
        rad.x /= norm;
        rad.y /= norm;
        #ifdef SQUARES
        rad = delta;
        #endif
        
        dPoint bl(tl.x - rad.x*width, tl.y - rad.y*width);
        p = project(bl.x, bl.y);
        fprintf(fout, "%d,%d ", p.x, p.y);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
        
        dPoint tr(omx.x + 0.5*width*tang.x, omx.y + 0.5*width*tang.y);
        
        rad = tr - centre;
        norm = sqrt(rad.x*rad.x + rad.y*rad.y);
        rad.x /= norm;
        rad.y /= norm;
        #ifdef SQUARES
        rad = delta;
        #endif
        
        dPoint br(tr.x - rad.x*width, tr.y - rad.y*width);
        p = project(br.x, br.y);
        fprintf(fout, "%d,%d ", p.x, p.y);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
        
        p = project(tr.x, tr.y);
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", p.x, p.y, style.c_str());
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
    }

    void grid(double swidth, size_t nrows, size_t ncols) {
        const int nangles = 2;
        
        for (int a=0; a < nangles; a++) {
            double angles[nangles] = {1.5/180.0*M_PI, -1.5/180.0*M_PI};
            double sgn = angles[a] > 0 ? 1 : -1;
            
            for (double rad=0.2; rad < 0.95; rad += 0.175) {
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a]),
                    (sqrt(0.5) + rad*0.5*sin(angles[a])) + sgn * swidth * 0.75,
                    swidth
                );
                
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a] + M_PI),
                    (sqrt(0.5) + rad*0.5*sin(angles[a] + M_PI)) - sgn * swidth * 0.75,
                    swidth
                );
            }
            
            for (double rad=0.2; rad < sqrt(2)*0.95; rad += 0.175) {

                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a] + M_PI/2) - sgn * swidth * 0.75,
                    (sqrt(0.5) + rad*0.5*sin(angles[a] + M_PI/2)),
                    swidth
                );
                
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a] + 3*M_PI/2) + sgn * swidth * 0.75,
                    (sqrt(0.5) + rad*0.5*sin(angles[a] + 3*M_PI/2)),
                    swidth
                );
            }
        }    
        for (double rad=0.38; rad < 1.1; rad += 0.16) {    
            double angles[4] = {45.0/180.0*M_PI, 135.0/180.0*M_PI, 225.0/180.0*M_PI, 315.0/180.0*M_PI};
            
            for (int a=0; a < 4; a++) {
                dPoint slant(-sin(angles[a]) * 0.75 * swidth, cos(angles[a]) * 0.75 * swidth);
            
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a]) - slant.x,
                    (sqrt(0.5) + rad*0.5*sin(angles[a])) - slant.y,
                    swidth
                );
                
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a]) + slant.x,
                    (sqrt(0.5) + rad*0.5*sin(angles[a])) + slant.y,
                    swidth
                );
            }
        }
        
        for (double rad=0.5; rad < 1.4; rad += 0.16) {    
            double angles[4] = {65.0/180.0*M_PI, 115.0/180.0*M_PI, 245.0/180.0*M_PI, 295.0/180.0*M_PI};
            
            for (int a=0; a < 4; a++) {
                dPoint slant(-sin(angles[a]) * 0.0 * swidth, cos(angles[a]) * 0.0 * swidth);
            
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a]) - slant.x,
                    (sqrt(0.5) + rad*0.5*sin(angles[a])) - slant.y,
                    swidth
                );
                
                /*
                place_trapezoid(
                    0.5 + rad*0.5*cos(angles[a]) + slant.x,
                    (sqrt(0.5) + rad*0.5*sin(angles[a])) + slant.y,
                    swidth
                );
                */
            }
        }
        
        iPoint p = project(0.985, sqrt(2)*0.985);
        fprintf(stderr, "3\n%lf %lf\n", p.y*scale + off, p.x*scale + off);
        p = project(0.987, sqrt(2)*0.985);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
        p = project(0.987, sqrt(2)*0.9865);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
    }
    
    dPoint centre;
    double off;
    double scale;
};

#endif
