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
#ifndef SVG_PAGE_PERSPECTIVE_H
#define SVG_PAGE_PERSPECTIVE_H

#include "svg_page.h"

class Svg_page_perspective : public Svg_page {
  public:
    Svg_page_perspective(const string& page_spec, const string& fname) 
      : Svg_page(page_spec, fname), cop(0,0,0), pl(0,0,2000), normal(0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)), fl(50), sensor(15.6, 23.6) {
        set_viewing_parameters(pl[2], fl, 15.6);
        compute_perspective_extremes();
    }
    
    void render(void) {
        double ang=10.0/180.0*M_PI;
        double bsize=0.0214; // block size ....
        
        strip(-6*2*1.5*2*bsize, -30*bsize, 2.5*bsize, 2*bsize, 15, 5, -ang, 1);
        strip( 6*2*1.5*2*bsize, -30*bsize, 2.5*bsize, 2*bsize, 15, 5, -ang, 0);
        perspective_rectangle(-0.2, 0.5, 0.4, 0.5, 1.0/180.0*3.14159, 0);
    }
    
  protected:  
  
    dPoint project_core(double x, double y) {
        //Vec3d d(x, y, fl/4.0); // fudge factor --- this definitely still needs work
        
        Vec3d d(x/(0.5*sensor[1]), y/(0.4*sensor[1]), focal_distance);
        
        d = d * (1.0/norm(d));
        double t = ((pl - cop).ddot(normal)) / (d.dot(normal));
        Vec3d pi = cop + d * t;
        
        Vec3d cr = normal.cross(Vec3d(1,0,0));
        dPoint p;
        p.y = (pi - pl).dot(cr);
        p.x = (pi - pl).dot(Vec3d(1,0,0));
        
        return p;
    }
  
    virtual iPoint project(double x, double y) {
    
        dPoint p = project_core(x,y);
        
        p.x = floor(p.x/(2*extremes[0])*width + 0.5*width);
        p.y = floor(p.y/(2*extremes[1])*height + 0.5*height);
        
        return iPoint(int(p.x), int(p.y));
    }
    
    void set_viewing_parameters(double distance_to_target, double focal_length, double sensor_height, double angle=-45.0/180*M_PI) {
        pl[2] = distance_to_target;
        fl = focal_length;
        normal[1] = cos(angle);
        normal[2] = sin(angle);
        sensor[0] = sensor_height;
        
        // where should the focal plane be so that the image size matches the sensor size:
        const double desired_angle = 0.5 * 17.735 / 180.0 * M_PI;
        focal_distance = 0.5 * sensor_height / tan(desired_angle);
        
        printf("focal distance = %lf\n", focal_distance);
        printf("size: (%lf, %lf)\n", tan(desired_angle)*focal_distance*2, tan(0.5*26.56/180.0*M_PI)*focal_distance*2);
        
        compute_perspective_extremes();
    }
    
    void perspective_rectangle(double tlx, double tly, double width, double height, double angle, bool right) {
        if (right) {
            iPoint p = project(tlx, tly);
            fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
            p = project(tlx - cos(angle)*width, tly + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx - sin(angle)*width - cos(angle)*width, tly - cos(angle)*height + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx - sin(angle)*width, tly - cos(angle)*height);
            fprintf(fout, "%d,%d\" style=\"%s\"/>", p.x, p.y, style.c_str());
          
        } else {
            iPoint p = project(tlx, tly);
            fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
            p = project(tlx + cos(angle)*width, tly + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx + sin(angle)*width + cos(angle)*width, tly - cos(angle)*height + sin(angle)*height);
            fprintf(fout, "%d,%d ", p.x, p.y);
            p = project(tlx + sin(angle)*width, tly - cos(angle)*height);
            fprintf(fout, "%d,%d\" style=\"%s\"/>", p.x, p.y, style.c_str());
          
        }
    }
    
    void strip(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols, double ang, bool left) {
  
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
      
                if (left) {
                    double ypos = ry * 2.0 * sheight + tly;
                    double xpos = rx * 2.0 * swidth  + tlx;
              
              
                    perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 0);
                    fprintf(fout, "\n");
                    if (ry < nrows-1) {
                        xpos += swidth;
                        ypos += sheight;
                        perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 0);
                        fprintf(fout, "\n");
                    }
                } else {
                    double ypos = ry * 2.0 * sheight + tly;
                    double xpos = tlx - rx * 2.0 * swidth;
              
                    perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 1);
                    fprintf(fout, "\n");
                    if (ry < nrows-1) {
                        xpos -= swidth;
                        ypos += sheight;
                        perspective_rectangle(xpos, ypos, swidth*0.8, sheight*0.8, ang, 1);
                        fprintf(fout, "\n");
                    }
                }
      
            } // columns
        } // rows
    }
    
    void compute_perspective_extremes(void) {    
    
        dPoint p1 = project_core(0, -1);
        dPoint p2 = project_core(0, 1);
        
        extremes[1] = std::max(fabs(p1.y), fabs(p2.y)) * width / double(height);
        
        p1 = project_core(-1, 0);
        p2 = project_core( 1, 0);
        
        extremes[0] = std::max(fabs(p1.x), fabs(p2.x));
        
        printf("extremes: %le, %le\n", extremes[0], extremes[1]);
    }
    
    Vec3d cop;      // centre of projection
    Vec3d pl;       // position of focal plane
    Vec3d normal;     // normal vector of target
    double fl;          // focal length of lens
    Vec2d extremes;
    Vec2d sensor;
    double focal_distance;
};

#endif
