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
