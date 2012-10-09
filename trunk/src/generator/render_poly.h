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
#ifndef RENDER_POLY_H
#define RENDER_POLY_H


#include <cv.h>
using namespace cv;

#include "render.h"

//==============================================================================
class Render_rectangle_poly : public Render_rectangle {
  public:
    Render_rectangle_poly(double cx, double cy, double width, double height, double angle) 
    : Render_rectangle(cx, cy, width, height, angle, false) {
    }
    
    virtual ~Render_rectangle_poly(void) {
    }
    
    double evaluate_x(const Render_rectangle& b, double xoffset = 0, double yoffset = 0)  const {
        double rvalue = compute_area(intersect(b, xoffset, yoffset));
        return rvalue;
    }
    
    // intersection of two rays defined by origin (v1 and v2) and directions (d1 and d2)
    bool intersect(const Vec2d& v1, const Vec2d& d1, 
                   const Vec2d& v2, const Vec2d& d2, 
                   Vec2d& pi) const {
                   
        double denom = (d2[1]*d1[0] - d2[0]*d1[1]);
        
        if (fabs(denom) < 1e-11) {
           pi[0] = pi[1] = 0;
           return false;
        }
        
        double u;
        u = (d2[0]*(v1[1] - v2[1]) - d2[1]*(v1[0] - v2[0])) / denom;
        
        pi = v1 + u*d1;
                   
        return true;               
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
    
    double compute_area(const vector<Vec2d>& points) const {
        double A = 0;
        int N = (int)points.size();
        for (int i=0; i < N; i++) {
            int ni = (i+1) % N;
            A += points[i][0]*points[ni][1] - points[ni][0]*points[i][1];
        }
        return 0.5 * fabs(A);
    }
    
    vector<Vec2d> intersect(const Render_rectangle& b, double xoffset = 0, double yoffset = 0) const {
        vector<Vec2d> points;
        
        for (int i=0; i < 4; i++) {
            points.push_back(cv::Point_<double>(b.bases[i][0] + xoffset, b.bases[i][1] + yoffset));
        }
        
        for (int e=0; e < 4; e++) {
            points = intersect_core(points, e);
        }
        
        return points;    
    }
    
    vector<Vec2d> intersect_core(const vector<Vec2d>& inpoints, int e) const {
        int ne = (e + 1) % 4;
        
        Vec2d P = bases[e];
        Vec2d D = bases[ne] - bases[e];
        
        Vec2d N = Vec2d(-D[1], D[0]);
        
        vector<Vec2d> outpoints;
        
        Vec2d S = inpoints[inpoints.size() - 1];
        Vec2d dS = S - P;
        for (size_t i=0; i < inpoints.size(); i++) {
            const Vec2d& E = inpoints[i];
            Vec2d dE = E - P;
            
            if (N.dot(dE) >= 0) {
                if (N.dot(dS) < 0) {
                    outpoints.push_back( t_intersect(S, E-S, P, D) );
                }
                outpoints.push_back(E);
            } else {  
                if (N.dot(dS) >= 0) {
                    outpoints.push_back( t_intersect(S, E-S, P, D) );
                }
            }
            S = E;
            dS = dE;
        }
    
        return outpoints;
    }
    
      
};

#endif // RENDER_H
