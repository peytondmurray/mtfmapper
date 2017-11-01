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
#ifndef MTF_RENDERER_EDGES_H
#define MTF_RENDERER_EDGES_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_edges : public Mtf_renderer {
  public:
    Mtf_renderer_edges(const std::string& fname, 
      const std::string& sfrname,
      bool lpmm_mode=false, double pixel_size=1.0) 
      :  ofname(fname), sfrname(sfrname),
         lpmm_mode(lpmm_mode), pixel_size(pixel_size) {
      
    }
    
    void render(const vector<Block>& blocks) {
        Point2d centr(0,0);
        for (size_t i=0; i < blocks.size(); i++) {
            centr += blocks[i].get_centroid();
        }
        centr = centr*(1.0/double(blocks.size()));
        
    
    
        FILE* fout = fopen(ofname.c_str(), "wt");
        FILE* sfrout = fopen(sfrname.c_str(), "wt");
        vector<int> corder(4);
        vector<int> eorder(4);
        for (size_t i=0; i < blocks.size(); i++) {
        
            // start with first corner that is left of 12:00
            for (size_t k=0; k < 4; k++) {
                Point2d dir = blocks[i].get_corner(k) - blocks[i].get_centroid();
                if (dir.x <= 0 && dir.y <= 0) {
                    corder[0] = k;
                }
                if (dir.x <= 0 && dir.y > 0) {
                    corder[1] = k;
                }
                if (dir.x > 0 && dir.y > 0) {
                    corder[2] = k;
                }
                if (dir.x > 0 && dir.y <= 0) {
                    corder[3] = k;
                }
            }
            
            // now find the first edge centroid going clockwise
            for (size_t k=0; k < 4; k++) {
                Point2d cdir = blocks[i].get_corner(corder[k]) - blocks[i].get_centroid();
                double maxcross = -1;
                for (size_t j=0; j < 4; j++) {
                    Point2d edir = blocks[i].get_edge_centroid(j) - blocks[i].get_centroid();
                    double n = cdir.x*edir.y - cdir.y*edir.x;
                    double dot = cdir.x*edir.x + cdir.y*edir.y;
                    if (n > maxcross && dot >= 0) {
                        eorder[k] = j;
                        maxcross = n;
                    }
                }
            }
            
            for (size_t k=0; k < 4; k++) {
                int j = corder[k];
                int l = eorder[k];
                double val = blocks[i].get_mtf50_value(l);
                Point2d ec = blocks[i].get_edge_centroid(l);
                Point2d cr = blocks[i].get_corner(j);
                fprintf(fout, "%d %lf %lf %lf %lf %lf\n",
                    int(i),
                    ec.x, ec.y,
                    lpmm_mode ? val*pixel_size : val,
                    cr.x, cr.y
                );
                
                fprintf(sfrout, "%d %lf %lf ", int(i), ec.x, ec.y);
                
                double edge_angle = atan2(-blocks[i].get_normal(l).x, blocks[i].get_normal(l).y);
                fprintf(sfrout, "%lf ", angle_reduce(edge_angle));
                
                Point2d cent = blocks[i].get_edge_centroid(l);

                Point2d dir = cent - centr;
                dir = dir * (1.0/norm(dir));

                Point2d norm = blocks[i].get_normal(l);
                double delta = dir.x*norm.x + dir.y*norm.y;
                
                fprintf(sfrout, "%lf ", acos(fabs(delta))/M_PI*180.0);
                
                const vector<double>& sfr = blocks[i].get_sfr(l);
                for (size_t j=0; j < sfr.size(); j++) {
                    fprintf(sfrout, "%lf ", sfr[j]);
                }
                fprintf(sfrout, "\n");
            }
        }    
        fclose(fout);
    }
    
  private:
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    string ofname;
    string sfrname;
    bool filter;
    double angle;
    bool    lpmm_mode;
    double  pixel_size;
};

#endif
