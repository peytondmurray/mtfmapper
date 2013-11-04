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
#ifndef MULTIPOLYGON_GEOM_H
#define MULTIPOLYGON_GEOM_H


#include <cv.h>
using namespace cv;

#include <vector>
using std::vector;

#include <algorithm>
using std::sort;

#include "geom.h"
#include "polygon_geom.h"

//==============================================================================
class Multipolygon_geom : public Geometry {
  public:

    Multipolygon_geom(double cx, double cy, const string& fname) 
    : Geometry(cx, cy, 0), bounding_box(0,0,1,1,0,4), total_vertices(0) {

        // we should subtract the centroid from the input polys, so that
        // cx,cy can point to the true centroid?

        FILE* fin = fopen(fname.c_str(), "rt");
        if (!fin) {
            fprintf(stderr, "Error. Could not open polygon geometry file %s. Aborting\n", fname.c_str());
            exit(1);
        }

        // TODO: maybe one day this can read SVG?
        while (!feof(fin)) {
            int nverts = 0;
            int nread = fscanf(fin, "%d", &nverts);
            if (nread == 1) {
                vector<cv::Vec2d> verts(nverts);
                for (int i=0; i < nverts; i++) {
                    nread = fscanf(fin, "%lf %lf", &verts[i][0], &verts[i][1]);
                    verts[i][0] += cx; // TODO: hack to centre the polygon on (cx,cy) ?
                    verts[i][1] += cy;
                }
                total_vertices += nverts;
                parts.push_back(Polygon_geom(verts));
            }
        }

	compute_bounding_box();

        own_area = 1;
    }
    
    Multipolygon_geom(void) {
        own_area = 1;
    }
    
    virtual ~Multipolygon_geom(void) {
    }

    void compute_bounding_box(void) {
        // TODO: we can be smarter here by computing eigenvectors
        // and aligning the box with the object ...

        // compute a bounding box
        bounding_box.bases[0][1] = 1e12;
        bounding_box.bases[1][0] = -1e12;
        bounding_box.bases[2][1] = -1e12;
        bounding_box.bases[3][0] = 1e12;
        for (size_t p=0; p < parts.size(); p++) {
            for (int i=0; i < parts[p].nvertices; i++) {
                // bounding_box.bases: 0=top, 1=right, 2=bottom, 3=left
                bounding_box.bases[0][1] = std::min(parts[p].bases[i][1], bounding_box.bases[0][1]);
                bounding_box.bases[1][0] = std::max(parts[p].bases[i][0], bounding_box.bases[1][0]);
                bounding_box.bases[2][1] = std::max(parts[p].bases[i][1], bounding_box.bases[2][1]);
                bounding_box.bases[3][0] = std::min(parts[p].bases[i][0], bounding_box.bases[3][0]);
            }
        }
        
        bounding_box.bases[0][0] = bounding_box.bases[3][0];
        bounding_box.bases[1][1] = bounding_box.bases[0][1];
        bounding_box.bases[2][0] = bounding_box.bases[1][0];
        bounding_box.bases[3][1] = bounding_box.bases[2][1];

        bounding_box.normals[0][0] = 0;
        bounding_box.normals[0][1] = -1;
        bounding_box.normals[1][0] = 1;
        bounding_box.normals[1][1] = 0;
        bounding_box.normals[2][0] = 0;
        bounding_box.normals[2][1] = 1;
        bounding_box.normals[3][0] = -1;
        bounding_box.normals[3][1] = 0;
    }

    bool is_inside(double x, double y) const {
        bool inside = false;

        for (size_t p=0; p < parts.size(); p++) {
            inside |= parts[p].is_inside(x, y);
        }

        return inside;
    }
    
    double intersection_area(const Geometry& b, double xoffset = 0, double yoffset = 0)  const {

        double bounds_area = 1;
        if (total_vertices > 6) {
            bounds_area = b.intersection_area(bounding_box, xoffset, yoffset);
        }

        double total_area = 0;
        if (bounds_area > 1e-11) {
            for (size_t p=0; p < parts.size(); p++) {
                total_area += b.intersection_area(parts[p], xoffset, yoffset);
            }
        }
        
        return total_area;
    }
    

    double cx;
    double cy;

    Polygon_geom bounding_box;

    vector<Polygon_geom> parts;
    int total_vertices;
};

#endif // MULTIPOLYGON_GEOM_H
