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
#ifndef SVG_PAGE_GRID_H
#define SVG_PAGE_GRID_H

#include "svg_page.h"

class Svg_page_grid : public Svg_page {
  public:
    Svg_page_grid(const string& page_spec, const string& fname) : Svg_page(page_spec, fname) {
    
    }
    
    void render(void) {
        double bsize=0.04; // block size ....
        grid(0.1, 0.1, bsize, bsize, 11, 11);
    }
    
  protected:
    void grid(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols) {
    
        double allowed_angles[] = {-4,4,-6,6,-86,86,-94,94,-96,96,-84,84};
    
        dPoint centre(0.5, 0.5);
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
      
                double ypos = ry * 2.0 * sheight + tly;
                double xpos = rx * 2.0 * swidth  + tlx;
                
                dPoint delta = dPoint(xpos, ypos) - centre;
                
                double eang = atan2(delta.y, delta.x) + M_PI/4.0;
                
                
                size_t m_idx = 0;
                for (size_t q=1; q < 12; q++) {
                    if ( fabs(allowed_angles[q]/180.0*M_PI - eang) < fabs(allowed_angles[m_idx]/180.0*M_PI - eang)) {
                        m_idx = q;
                        
                    }
                }
                eang = (allowed_angles[m_idx] + 45)/180.0*M_PI;
                
          
                rotated_square(xpos, ypos, swidth*0.5, eang);
                fprintf(fout, "\n");
                if (ry < nrows-1) {
                    xpos += swidth;
                    ypos += sheight;
                    rotated_square(xpos, ypos, swidth*0.5, eang);
                    fprintf(fout, "\n");
                }
            }
        } 
    }
};

#endif
