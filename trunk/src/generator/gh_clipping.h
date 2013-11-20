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
#ifndef GH_CLIPPING_H
#define GH_CLIPPING_H


#include <cv.h>
using namespace cv;

#include <vector>
using std::vector;

#include "geom.h"

class Polygon_geom;

//==============================================================================
namespace GH_clipping {

    struct gh_vertex {
        double x;
        double y;
        int next;       // we use integers here to avoid multiple dynamic allocations
        int prev;       // and the corresponding pain of freeing up the linked list again
        bool isect;
        int  en;        // intersections: (0 = exit, 1 = entry), vertices(0 = outside, 1 = inside, 2 = on)
        double alpha;
        int neighbour;
        int next_poly;  // we also use a global vector for storing the polygons ... but there are only two?
    };

                      
    int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly);
    void gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size);
    void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index);
    void gh_phase_three(vector<gh_vertex>& verts, int vs, int first_isect_index, vector<Polygon_geom>& polys);
    
};

#endif // GH_CLIPPING_H
