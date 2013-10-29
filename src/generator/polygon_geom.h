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
#ifndef POLYGON_GEOM_H
#define POLYGON_GEOM_H


#include <cv.h>
using namespace cv;

//#include "render.h"



//==============================================================================
class Polygon_geom  {
  public:
    static const int max_verts_per_poly = 70; // max number of point in intersected polygon
  
    Polygon_geom(double cx, double cy, double width, double height, double angle, int nverts=4) 
    : cx(cx), cy(cy), nvertices(nverts) {
      //Render_polygon(cx, cy, width, height, angle, 6.0, 6.0, 0.0, false, nverts),

        construct_regular_polygon(width, height, angle);
        own_area = compute_area();
        printf("built polygon with %d vertices, area=%lf\n\n", nverts, own_area);
    }
    
    virtual ~Polygon_geom(void) {
    }

    void construct_regular_polygon(double width, double height, double angle) {
        printf("nvertices=%d\n", nvertices);
        if (false /*nvertices == 4*/) {
            printf("rendering a plain square\n");
            // default is a square
            bases[0] = cv::Vec2d(width/2, height/2);
            bases[1] = cv::Vec2d(-width/2, height/2);
            bases[2] = cv::Vec2d(-width/2, -height/2);
            bases[3] = cv::Vec2d(width/2, -height/2);
            
            // rotate corners
            for (size_t i=0; i < 4; i++) {
                double ox = bases[i][0];
                double oy = bases[i][1];
                bases[i][0] = cos(angle)*ox - sin(angle)*oy + cx;
                bases[i][1] = sin(angle)*ox + cos(angle)*oy + cy;
            }
                  
            normals[0] = (bases[2] - bases[1]); 
            normals[1] = (bases[3] - bases[2]); 
            normals[2] = (bases[0] - bases[3]); 
            normals[3] = (bases[1] - bases[0]); 
            for (size_t i=0; i < 4; i++) {
                double n = norm(normals[i]);
                normals[i] = normals[i]*(1.0/n);
            }
            
            for (int i=0; i < nvertices; i++) {
                printf("** %lf %lf (%lf %lf)\n", bases[i][0], bases[i][1], normals[i][0], normals[i][1]);
            }
        } else {
            printf("rendering a polygon with %d sides\n", nvertices);
            assert(nvertices >= 3);
            for (int i=0; i < nvertices; i++) {
                double phi = i*M_PI*2/double(nvertices);
                bases[i][0] = width/2*cos(angle+phi - M_PI/4.0) + cx;
                bases[i][1] = height/2*sin(angle+phi - M_PI/4.0) + cy;
                normals[i][0] = cos(angle+phi);
                normals[i][1] = -sin(angle+phi);
                printf("%lf %lf (%lf %lf)\n", bases[i][0], bases[i][1], normals[i][0], normals[i][1]);
            }
            printf("horizontal width: %lf\n", norm(bases[1] - bases[0]));
        }
    }

    inline bool is_inside(double x, double y) const {
        bool inside = true;
        for (int i=0; i < nvertices && inside; i++) {
            double dot = normals[i][0]*(x - bases[i][0]) + 
                         normals[i][1]*(y - bases[i][1]);
            if (dot < 0) {
                inside = false;
            }
        }
        return inside;
    }
    
    double evaluate_x(const Polygon_geom& b, double xoffset = 0, double yoffset = 0)  const {
    
        double points_x[max_verts_per_poly];
        double points_y[max_verts_per_poly];
        int points_len = 0;
        
        intersect(points_x, points_y, points_len, b, xoffset, yoffset);
        return compute_area(points_x, points_y, points_len);
    }
    
    Vec2d t_intersect(const Vec2d& v1, const Vec2d& d1, 
                   const Vec2d& v2, const Vec2d& d2) const {
                   
        double denom = (d2[1]*d1[0] - d2[0]*d1[1]);
        
        if (fabs(denom) < 1e-11) {
           printf("denom zero?\n");
           printf("1: (%lf, %lf), (%lf, %lf)\n", v1[0], v1[1], d1[0], d1[1]);
           printf("2: (%lf, %lf), (%lf, %lf)\n", v2[0], v2[1], d2[0], d2[1]);
           return Vec2d(0.0, 0.0);
        }
        
        double u = (d2[0]*(v1[1] - v2[1]) - d2[1]*(v1[0] - v2[0])) / denom;
        
        Vec2d pi = v1 + u*d1;
                   
        return pi;               
    }
    
    inline bool t_intersect(double& pix, double& piy, 
                     const double& v1x, const double& v1y,
                     const double& d1x, const double& d1y,
                     const double& v2x, const double& v2y,
                     const double& d2x, const double& d2y) const {
                   
        double denom = (d2y*d1x - d2x*d1y);
        
        if (fabs(denom) < 1e-11) {
           printf("denom zero?\n");
           return false;
        }
        
        double u = (d2x*(v1y - v2y) - d2y*(v1x - v2x)) / denom;
        
        pix = v1x + u*d1x;
        piy = v1y + u*d1y;
                   
        return true;               
    }
    
    double compute_area(double* points_x, double* points_y, int points_len) const {
        double A = 0;
        for (int i=0; i < points_len; i++) {
            int ni = (i+1) % points_len;
            A += points_x[i]*points_y[ni] - points_x[ni]*points_y[i];
        }
        return 0.5 * fabs(A);
    }
    
    double compute_area(void) const {
        double A = 0;
        for (int i=0; i < nvertices; i++) {
            int ni = (i+1) % nvertices;
            A += bases[i][0]*bases[ni][1] - bases[ni][0]*bases[i][1];
        }
        return 0.5 * fabs(A);
    }
    
    void intersect(double* points_x, double* points_y, int& points_len, 
        const Polygon_geom& b, double xoffset = 0, double yoffset = 0) const {
        
        for (int i=0; i < b.nvertices; i++) { 
            points_x[i] = b.bases[i][0] + xoffset;
            points_y[i] = b.bases[i][1] + yoffset;
            points_len++;
        }
        
        for (int e=0; e < nvertices; e++) {
            intersect_core(points_x, points_y, points_len, e, nvertices);
        }
    }

    // this looks like the Sutherland-Hodgman algorithm
    // the clipping polygon must be convex (req. by SH algo)
    // will produce overlapping edges if a concave point exists
    // outside of the clipping polygon (this wiil definitely happen
    // if concave object polygon is presented). The solution is to
    // decompose all input polygons into convex components ...
    void intersect_core(double* inpoints_x, double* inpoints_y, int& in_len, int e, int nedges) const {
        int ne = (e + 1) % nedges;
        
        double Px = bases[e][0];
        double Py = bases[e][1];
        
        double Dx = bases[ne][0] - Px;
        double Dy = bases[ne][1] - Py;
         
        double Nx = -Dy;
        double Ny = Dx;
        
        double outpoints_x[max_verts_per_poly];
        double outpoints_y[max_verts_per_poly];
        int out_idx = 0;
        
        double Sx = inpoints_x[in_len - 1];
        double Sy = inpoints_y[in_len - 1];
        
        double dSx = Sx - Px;
        double dSy = Sy - Py;
        
        for (size_t i=0; i < size_t(in_len); i++) {
            const double& Ex = inpoints_x[i];
            const double& Ey = inpoints_y[i];
            
            double dEx = Ex - Px;
            double dEy = Ey - Py;
            
            if ( (Nx*dEx + Ny*dEy) >= 0) {
                if ( (Nx*dSx + Ny*dSy) < 0) {
                    t_intersect(
                        outpoints_x[out_idx], outpoints_y[out_idx],
                        Sx, Sy,
                        Ex - Sx, Ey - Sy,
                        Px, Py,
                        Dx, Dy
                    );
                    out_idx++;
                }
                outpoints_x[out_idx] = Ex;
                outpoints_y[out_idx] = Ey;
                out_idx++;
            } else {  
                if ( (Nx*dSx + Ny*dSy) >= 0) {
                    t_intersect(
                        outpoints_x[out_idx], outpoints_y[out_idx],
                        Sx, Sy,
                        Ex - Sx, Ey - Sy,
                        Px, Py,
                        Dx, Dy
                    );
                    out_idx++;
                }
            }
            Sx = Ex;
            Sy = Ey;
            dSx = dEx;
            dSy = dEy;
        }
    
        in_len = out_idx;
        memcpy(inpoints_x, outpoints_x, sizeof(double)*out_idx);
        memcpy(inpoints_y, outpoints_y, sizeof(double)*out_idx);
    }
    
    double own_area;

    cv::Vec2d normals[max_verts_per_poly];
    cv::Vec2d bases[max_verts_per_poly];

    double cx;
    double cy;
    int nvertices;
};

#endif // RENDER_H
