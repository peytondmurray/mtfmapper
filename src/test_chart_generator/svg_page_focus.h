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
#ifndef SVG_PAGE_FOCUS_H
#define SVG_PAGE_FOCUS_H

#include "svg_page.h"

#include <string>

#include "include/fiducial_positions.h"

using std::string;
using std::ostringstream;

class Svg_page_focus : public Svg_page {
  public:
    Svg_page_focus(const string& page_spec, const string& fname) 
      : Svg_page(page_spec, fname, 10000), cop(0,0,0), pl(0,0,2000), 
        normal(0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)), fl(50), 
        sensor(15.6, 23.6), outside_limits(false), width_scale(1), clipid(1) {
        
        // must call set_viewing_parameters before render
    }
    
    void render(void) {
        double ang=5.0/180.0*M_PI;
        double bsize=0.005; // block size ....
        
        //strip(-bsize, -2.8*10*bsize, 2.5*bsize, 10*bsize, 7, 1, -ang);
        //strip(-5, -150, 10, 60, 7, 1, -ang); // for direct blocks, non-perspective
        //strip(-5, -150, 45, 9, 25, 1, -ang); // for direct blocks
        strip(-5*bsize, -0.8, 45*bsize, 5*bsize, 44, 1, -ang); // for perspective blocks
        //codedcircles(5);
        
        coded_sector_circles(10);
        
        chevrons();
        
        // print out chart specification
        fprintf(fout, "\n");
        text(chart_description, 0.1, 0.97, 30);
        
        if (outside_limits) {
            printf("Warning: chart size exceeds the dimensions of the page you have specified.\n");
            printf("\tThe most likely cause of this is that the specified distance is too small,\n");
            printf("\ti.e., you are too close to the chart; increase distance with -d.\n");
        }
    }
    
    void set_viewing_parameters(double distance_to_target, double angle=-45.0/180*M_PI) {
        pl[2] = distance_to_target;
        normal[1] = cos(angle);
        normal[2] = sin(angle);
        
        double effective_angle = 2*atan(width_mm/(2*distance_to_target));
        
        width_scale = fl * tan(effective_angle/2);
        double f = 15.6/(2*tan(effective_angle/2));
        printf("maximum focal length (15.6 mm sensor height): %.2lf mm\n", f);
        printf("maximum focal length (24 mm sensor height): %.2lf mm\n", 24/(2*tan(effective_angle/2)));
        printf("distance-to-focal-length ratio: %lfx\n", pl[2]/f);
        
        ostringstream ss;
        ss << "Chart size: " << page_size << ", ";
        ss << "design viewing distance: " << pl[2] << " mm, ";
        ss << "design viewing angle: " << lrint(angle/M_PI*1800)/10.0 << " degrees, ";
        ss << "maximum focal length (sensor height = " << 15.6 << " mm): " << lrint(f*10)/10.0 << " mm, ";
        ss << "ratio at maximum focal length: " << lrint(pl[2]/f*10)/10.0;
        chart_description = ss.str();
    }
    
  protected:  
  
    dPoint project_core(double x, double y) {
        Vec3d d(x*width_scale, y*width_scale, fl);
        
        d = d * (1.0/norm(d));
        double t = ((pl - cop).ddot(normal)) / (d.dot(normal));
        Vec3d pi = cop + d * t;
        
        Vec3d cr = normal.cross(Vec3d(1,0,0)); // basis vector in plane, along y direction
        dPoint p;
        p.y = (pi - pl).dot(cr);
        p.x = (pi - pl).dot(Vec3d(1,0,0));
        
        return p;
    }
  
    virtual iPoint project(double x, double y) {
    
        dPoint p = project_core(x,y);
        
        p.x = floor(p.x*sscale + 0.5*width);
        p.y = floor(p.y*sscale + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = true;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    virtual iPoint scale(double x, double y) {
    
        dPoint p(x,y);
        
        p.x = floor(p.x*sscale + 0.5*width);
        p.y = floor(p.y*sscale + 0.5*height);
        
        if (p.x < 0 || p.x > width ||
            p.y < 0 || p.y > height) {
            
            outside_limits = true;
        }
        
        return iPoint(int(p.x), int(p.y));
    }
    
    
    void perspective_rectangle(double cx, double cy, double width, double height, double angle) {
        fprintf(fout, "  <polygon points=\"");
        for (int i=0; i < 4; i++) {
            double theta = 2*M_PI*i/double(4) + M_PI/4.0; 
            double urx = 0.5*width*cos(theta);
            double ury = 0.5*height*sin(theta);
            double rx = urx*cos(angle) - ury*sin(angle);
            double ry = urx*sin(angle) + ury*cos(angle);
            iPoint p = project(cx + rx, cy + ry);
            fprintf(fout, "%d,%d ", p.x, p.y);
        }
        fprintf(fout, "\" style=\"%s\"/>", style.c_str());
    }
    
    void strip(double tlx, double tly, double swidth, double sheight, size_t nrows, size_t ncols, double ang) {
        
        for (size_t ry=0; ry < nrows; ry++) {
            for (size_t rx=0; rx < ncols; rx++) {
      
            
                double ypos = ry * 1.4 * sheight + tly; // perspective
                //double ypos = ry * 0.8 * sheight + tly;   // direct
                double xpos = tlx - rx * 2.0 * swidth;
          
                //direct_rectangle(xpos, ypos, swidth, sheight, ang);
                perspective_rectangle(xpos, ypos, swidth, sheight, -ang);
                fprintf(fout, "\n");
      
            } // columns
        } // rows
    }
    
    void perspective_poly(double cx, double cy, double rad, int sides=64) {
            fprintf(fout, "  <polygon points=\"");
            for (int i=0; i < sides; i++) {
                double theta = 2*M_PI*i/double(sides+1);
                iPoint p = project(cx + rad*cos(theta), cy + rad*sin(theta));
                fprintf(fout, "%d,%d ", p.x, p.y);
            }
            fprintf(fout, "\" style=\"%s\"/>", style.c_str());
          
    }
    
    iPoint arc_poly_outer_hollow_bounds(double cx, double cy, double srad, double h, double xsign=1.0) {
        if (xsign > 0) {
            iPoint p = scale(cx + h, cy - srad);
            return iPoint(p.x, int(srad*sscale));
        } else {
            iPoint p = scale(cx - srad, cy - srad);
            iPoint p2 = scale(cx - h, cy + srad);
            return iPoint(p.x, p2.x - p.x);
        }
    }
    
    iPoint arc_poly_inner_hollow_bounds(double cx, double cy, double srad, double h, double xsign=1.0) {
        if (xsign > 0) {
            iPoint p = scale(cx, cy - srad);
            return iPoint(p.x, int(h*sscale));
        } else {
            iPoint p = scale(cx - h, cy - srad);
            iPoint p2 = scale(cx, cy + srad);
            return iPoint(p.x, p2.x - p.x);
        }
    }
    
    void emit_clipper(iPoint& current, iPoint range, const iPoint& ybounds) {
        if (current.x >= 0) { 
            fprintf(fout, " <rect x=\"%d\" y=\"%d\" width=\"%d\"  height=\"%d\"/>", current.x, ybounds.x, current.y, ybounds.y);
        }
        current = range;
    }
    
    void coded_poly(double cx, double cy, double rad, int code) {
    
        // black dot
        iPoint p = scale(cx, cy);
        fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"black\" />\n", p.x, p.y, int(rad*sscale));
        
        if (code == 0) return; // no further clipping
        
        fprintf(fout, "  <defs> <clipPath id=\"f%04d\"> ", clipid);
        double h = rad/1.4;
        double htheta = acos(1 - h/rad);
        
        double xlim = rad*cos(htheta);
        double srad = 0.6*rad;
        
        iPoint ybounds(int(scale(0, cy - srad).y), int(scale(0, cy + srad).y));
        iPoint current(-100, -100);
        if (code & 8) {
            current = arc_poly_outer_hollow_bounds(cx, cy, srad, xlim, -1.0);
        }
        if (code & 4) {
            iPoint range = arc_poly_inner_hollow_bounds(cx, cy, srad, xlim, -1.0);
            if (range.x == (current.x + current.y)) {
                current.y += range.y;
            } else {
                emit_clipper(current, range, ybounds);
            }
        }
        if (code & 2) {
            iPoint range = arc_poly_inner_hollow_bounds(cx, cy, srad, xlim);
            if (range.x == (current.x + current.y)) {
                current.y += range.y;
            } else {
                emit_clipper(current, range, ybounds);
            }
        }
        if (code & 1) {
            iPoint range = arc_poly_outer_hollow_bounds(cx, cy, srad, xlim);
            if (range.x == (current.x + current.y)) {
                current.y += range.y;
            } else {
                emit_clipper(current, range, ybounds);
            }
        }
        emit_clipper(current, current, ybounds);
        fprintf(fout, " </clipPath> </defs>\n");
        
        p = scale(cx, cy);
        fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"white\" clip-path=\"url(#f%04d)\" />\n", p.x, p.y, int(0.6*rad*sscale), clipid);
        clipid++;
    }
    
    void wedge(double cx, double cy, double inner_rad, double outer_rad, double xangle, double fraction, int depth=1) {
        double epsilon = 0;
        double epsilon2 = 0;
        
        if (fraction >= 0.999999) {
            iPoint p = scale(cx, cy);
            fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"black\" z-index=\"%d\"/>\n", p.x, p.y, int(outer_rad*sscale), depth);
            fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"white\" z-index=\"%d\"/>\n", p.x, p.y, int(inner_rad*sscale), depth+1);
        } else {
        
            int lf = fraction > 0.5 ? 1 : 0;
            double x = cx - outer_rad * cos((xangle-epsilon2)/180.0*M_PI);
            double y = cy - outer_rad * sin((xangle-epsilon2)/180.0*M_PI);
            iPoint p = scale(x, y);
            fprintf(fout, "  <path \n");
            fprintf(fout, "\td=\"M %d,%d ", p.x, p.y);
            x = cx - outer_rad * cos((xangle + 360*fraction+epsilon2)/180.0*M_PI);
            y = cy - outer_rad * sin((xangle + 360*fraction+epsilon2)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "A %lf %lf 0 %d,1 %d,%d ", sscale*(outer_rad+epsilon), sscale*(outer_rad+epsilon), lf, p.x, p.y);
            x = cx - inner_rad * cos((xangle + 360*fraction+epsilon2)/180.0*M_PI);
            y = cy - inner_rad * sin((xangle + 360*fraction+epsilon2)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "L %d,%d ", p.x, p.y);
            x = cx - inner_rad * cos((xangle-epsilon2)/180.0*M_PI);
            y = cy - inner_rad * sin((xangle-epsilon2)/180.0*M_PI);
            p = scale(x, y);
            fprintf(fout, "A %lf,%lf 0 %d,0 %d,%d ", sscale*inner_rad, sscale*inner_rad, lf, p.x, p.y);
            fprintf(fout, "z\" \n\tstyle=\"%s\"\n z-index=\"%d\" />\n", style.c_str(), depth);
        }
    }

    
    void sector_circle(double cx, double cy, double rad, int code) {
        iPoint p = scale(cx, cy);
        fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"black\" z-index=\"1\" />\n", p.x, p.y, int(rad*sscale));
        
        const int outer_step = 5;
        const int inner_step = 3;
        
        const int outer_vals = outer_step + 1;
        const int inner_vals = inner_step + 1;
        
        int inner_code = inner_step - code % inner_vals;
        int outer_code = outer_step - (code / inner_vals) % outer_vals;
        
        if (outer_code > 0) {
            fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"white\" z-index=\"2\" />\n", p.x, p.y, int(0.6*rad*sscale));
            wedge(cx, cy, 0.4*rad, 0.6*rad, 90, outer_code*1.0/double(outer_step), 3);
        }
        if (inner_code > 0) {
            if (outer_code == 0) {
                fprintf(fout, "  <circle cx=\"%d\" cy=\"%d\" r=\"%d\" stroke=\"none\" fill=\"white\" z-index=\"2\" />\n", p.x, p.y, int(0.4*rad*sscale));
            }
            wedge(cx, cy, 0.1*rad, 0.4*rad, 270, inner_code*1.0/double(inner_step), 5);
        }
    }
    
    void coded_sector_circles(double swidth) {
        const double hshift = 20;
        
        #if 1
        
        const double cdist = 30;
        sector_circle(cdist, 0, swidth/2.0, 0);
        fprintf(fout, "\n");
        printf("{0, 0, %lf, %lf, 0, 0},\n", cdist, 0.0);
        sector_circle(-cdist, 0, swidth/2.0, 0);
        fprintf(fout, "\n");
        printf("{0, 0, %lf, %lf, 0, 0},\n", -cdist, 0.0);
        
        const double theta = 35/180.0*M_PI;
        const double radius = 90*sqrt(0.5);
        int code = 2;
        for (double phi=0; phi > -2*M_PI; phi -= M_PI/2.0) {
            double cx = radius*cos(phi + theta);
            double cy = radius*sin(phi + theta);
            
            sector_circle(cx, cy, swidth/2.0, code);
            fprintf(fout, "\n");
            
            int quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
            if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                quad = 1 - quad;
            }
            
            printf("{0, 0, %lf, %lf, %d, %d}%c\n", cx, cy, code, quad, quad == 3 ? ' ' : ','); // don't actually need quadrants for these four
            
            quad++;
            code += 2;
        }
        
        vector<double> base_rad{100, 120, 130, 140};
        vector<double> offset_angle{-35, 20, 70, 45};
        for (size_t i=0; i < base_rad.size(); i++) {
        
            double rad=0;
            for (double phi=0; phi > -2*M_PI; phi -= M_PI/2.0) {
                double cx = (base_rad[i]+rad)*cos(phi + offset_angle[i]/180.0*M_PI);
                double cy = (base_rad[i]+rad)*sin(phi + offset_angle[i]/180.0*M_PI);
                
                sector_circle(cx, cy, swidth/2.0, code);
                
                int quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
                if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                    quad = 1 - quad;
                }
                
                printf("{0, 0, %lf, %lf, %d, %d}%c\n", cx, cy, code, quad, quad == 3 ? ' ' : ','); // don't actually need quadrants for these four
                
                quad++;
                rad += 2.5;
            }
            
            code += 2;
        }
        
        #else 
        // just use what we have in fiducials include ...
        for (int i=0; i < n_fiducials; i++) { // defined in "include/fiducial_positions.h"
            coded_poly(main_fiducials[i].rcoords.x, main_fiducials[i].rcoords.y, swidth, main_fiducials[i].code);
            fprintf(fout, "\n");
        }
        #endif
        
        for (int y=-80; y < 80; y++) {
            iPoint c1 = scale(swidth-1 + hshift/2, 2*y);
            iPoint c2 = scale(swidth+1 + hshift/2, 2*y+1);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
        
        const double rwidth = 2;
        for (int x=0; x < 25; x++) {
            iPoint c1 = scale( (4*x-1)*rwidth*0.5 + 2*hshift, -0.25);
            iPoint c2 = scale( (4*x+1)*rwidth*0.5 + 2*hshift,  0.25);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
        for (int x=0; x < 25; x++) {
            iPoint c1 = scale( (4*(-x)-1)*rwidth*0.5 - 2*hshift, -0.25);
            iPoint c2 = scale( (4*(-x)+1)*rwidth*0.5 - 2*hshift, 0.25);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
    }
    
    void codedcircles(double swidth) {
        const double hshift = 20;
        
        for (int rx=-1; rx <= 1; rx++) {
            coded_poly(rx*9.1*swidth + hshift, 0, swidth/2.0, 0);
        }
        
        for (int y=1; y <= 12; y++) {
            coded_poly(hshift, -1.5*y*swidth, swidth/2.0, (y-1) % 15 + 1);
        }
        
        for (int y=1; y <= 12; y++) {
            coded_poly(hshift, 1.5*y*swidth, swidth/2.0, (y+1) % 15 + 1);
        }
        
        
        #if 0
        const double theta = 35/180.0*M_PI;
        const double radius = 90*sqrt(0.5);
        for (double phi=0; phi > -2*M_PI; phi -= M_PI/2.0) {
            double cx = radius*cos(phi + theta);
            double cy = radius*sin(phi + theta);
            
            coded_poly(cx, cy, swidth, 13);
            fprintf(fout, "\n");
            
            // A: 52.130516 -36.502182
            // B: -36.502182 -52.130516
            // C: -52.130516 36.502182
            // D: 36.502182 52.130516
            
            int quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
            if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                quad = 1 - quad;
            }
            
            printf("{0, 0, %lf, %lf, 13, %d},\n", cx, cy, quad);
            
            cx = radius*cos(phi - theta);
            cy = radius*sin(phi - theta);
            
            quad = ((cx < 0) ? 1 : 0) | ((cy > 0) ? 2 : 0);
            if (quad <= 1) { // ugly, but swap quadrants 0 and 1
                quad = 1 - quad;
            }
            
            coded_poly(cx, cy, swidth, 14);
            fprintf(fout, "\n");
            
            printf("{0, 0, %lf, %lf, 14, %d}%c\n", cx, cy, quad, quad == 3 ? ' ' : ',');
            
            quad++;
        }
        #else 
        // just use what we have in fiducials include ...
        for (int i=0; i < n_fiducials; i++) { // defined in "include/fiducial_positions.h"
            coded_poly(main_fiducials[i].rcoords.x, main_fiducials[i].rcoords.y, swidth, main_fiducials[i].code);
            fprintf(fout, "\n");
        }
        #endif
        
        /*
        // corner indicators
        for (int rx=-1; rx <= 1; rx+=2) {
            coded_poly(rx*27*swidth, 130, swidth, 14);
            fprintf(fout, "\n");
        }
        
        for (int rx=-1; rx <= 1; rx+=2) {
            coded_poly(rx*27*swidth, -150, swidth, 13);
            fprintf(fout, "\n");
        }
        */
        
        for (int y=-80; y < 80; y++) {
            iPoint c1 = scale(swidth-1 + hshift, 2*y);
            iPoint c2 = scale(swidth+1 + hshift, 2*y+1);
            fprintf(fout, "\n");
            rect(c1.x, c1.y, c2.x - c1.x, c2.y - c1.y);
        }
    }
    
    void direct_rectangle(double cx, double cy, double width, double height, double angle) {
            fprintf(fout, "  <polygon points=\"");
            for (int i=0; i < 4; i++) {
                double theta = 2*M_PI*i/double(4) + M_PI/4.0; 
                double urx = 0.5*width*cos(theta);
                double ury = 0.5*height*sin(theta);
                double rx = urx*cos(angle) - ury*sin(angle);
                double ry = urx*sin(angle) + ury*cos(angle);
                iPoint p = scale(cx + rx, cy + ry);
                fprintf(fout, "%d,%d ", p.x, p.y);
            }
            fprintf(fout, "\" style=\"%s\"/>", style.c_str());
          
    }
    
    void chevrons(void) {
        const int fine_rows = 320;
        const double hshift = 20;
        for (int row=0; row < fine_rows; row++) {
            direct_rectangle(-7 + hshift, row*1 - fine_rows/2, 3.5, 0.5, 45.0/180*M_PI);
            fprintf(fout, "\n");
        }
        for (int row=0; row < fine_rows/4; row++) {
            for (int col=0; col < 3; col ++) {
                direct_rectangle(col - 7 + hshift, row*4 - fine_rows/2, 3.5, 0.5, -45.0/180*M_PI);
                fprintf(fout, "\n");
            }
        }
    }
    
    Vec3d cop;      // centre of projection
    Vec3d pl;       // position of focal plane
    Vec3d normal;     // normal vector of target
    double fl;          // focal length of lens
    Vec2d sensor;
    bool outside_limits;
    double width_scale;
    std::string chart_description;
    int clipid;
};

#endif
