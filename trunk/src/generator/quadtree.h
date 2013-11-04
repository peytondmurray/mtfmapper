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
#ifndef QUADTREE_H
#define QUADTREE_H

#include "geom.h"
#include "multipolygon_geom.h"

//==============================================================================
class Quadtree : public Multipolygon_geom {
  public:
     // constructor for root of tree
    Quadtree(double cx, double cy, const string& fname)
     : Multipolygon_geom(cx, cy, fname),  
       q_tl(0), q_tr(0), q_bl(0), q_br(0), leaf(false) {
    
        partition_polygons(0);
    }

    // constructor for children      
    //Quadtree(const Polygon_geom& bounds, int level = 0) {
    //    partition_polygons(level+1);
    //}
    
    Quadtree(const Polygon_geom& bounds) 
    : q_tl(0), q_tr(0), q_bl(0), q_br(0), leaf(false) {
        
        for (size_t i=0; i < bounds.bases.size(); i++) {
            fprintf(stderr, "%lf %lf\n", bounds.bases[i][0], bounds.bases[i][1]);
        }
        fprintf(stderr, "%lf %lf\n\n", bounds.bases[0][0], bounds.bases[0][1]);
        
        bounding_box = bounds;
        total_vertices = 0;
    }
    
    void add_poly(const Polygon_geom& p) {
        total_vertices += p.bases.size();
        parts.push_back(p);
    }

    virtual double intersection_area(const Geometry& b, double xoffset = 0, double yoffset = 0) const {
        // first, check this level's bounding box
        
        double bounds_area = 1;
        bounds_area = b.intersection_area(bounding_box, xoffset, yoffset);
        
        if (bounds_area < 1e-11) {
            return bounds_area;
        }
          
        // if this is a leaf, check parts, else pass on to
        // children
        double area = 0;
        for (size_t p=0; p < parts.size(); p++) {
            area += b.intersection_area(parts[p], xoffset, yoffset);
        }
        if (q_tl) area += q_tl->intersection_area(b, xoffset, yoffset);
        if (q_tr) area += q_tr->intersection_area(b, xoffset, yoffset);
        if (q_bl) area += q_bl->intersection_area(b, xoffset, yoffset);
        if (q_br) area += q_br->intersection_area(b, xoffset, yoffset);
        return area;
    }
    
    virtual bool is_inside(double x, double y) const {
        printf("\n\nnot defined\n\n");
        return x + y; // just to keep the warnings down
    }
    
    void partition_polygons(int level = 0) {
        printf("entering partition at level %d with %d parts\n", level, (int)parts.size());
        
        vector<double> all_x;
        vector<double> all_y;

        if (true || level == 0) {
            for (size_t p=0; p < parts.size(); p++) {
                for (size_t v=0; v < parts[p].bases.size(); v++) {
                    all_x.push_back(parts[p].bases[v][0]);
                    all_y.push_back(parts[p].bases[v][1]);
                }
            }
        } else {
            for (size_t v=0; v < 4; v++) {
                all_x.push_back(bounding_box.bases[v][0]);
                all_y.push_back(bounding_box.bases[v][1]);
            }
        }

        sort(all_x.begin(), all_x.end());
        sort(all_y.begin(), all_y.end());

        printf("\nmidx=%lf, midy=%lf\n", all_x[all_x.size()/2], all_y[all_y.size()/2]);

        // somewhat verbose, but manually build the four quadrants

        Polygon_geom tl;
        Polygon_geom tr;
        Polygon_geom bl;
        Polygon_geom br;

        const double eps = 0; //1e-8;

        size_t last = all_x.size() - 1;

        cv::Vec2d v_min(all_x[0], all_y[0]);
        //cv::Vec2d v_mid(all_x[last/2] + 1e-6, all_y[last/2] + 1e-6);
        cv::Vec2d v_mid((all_x[last]+all_x[0])/2.0, (all_y[last]+all_y[0])/2.0 );
        cv::Vec2d v_max(all_x[last], all_y[last]);

        tl.bases[3] = cv::Vec2d(v_min[0], v_min[1]);
        tl.bases[2] = cv::Vec2d(v_min[0], v_mid[1]);
        tl.bases[1] = cv::Vec2d(v_mid[0], v_mid[1]);
        tl.bases[0] = cv::Vec2d(v_mid[0], v_min[1]);
        tl.rebuild();

        tr.bases[3] = cv::Vec2d(v_mid[0]+eps, v_min[1]);
        tr.bases[2] = cv::Vec2d(v_mid[0]+eps, v_mid[1]);
        tr.bases[1] = cv::Vec2d(v_max[0]+eps, v_mid[1]);
        tr.bases[0] = cv::Vec2d(v_max[0]+eps, v_min[1]);
        tr.rebuild();

        bl.bases[3] = cv::Vec2d(v_min[0], v_mid[1]+eps);
        bl.bases[2] = cv::Vec2d(v_min[0], v_max[1]+eps);
        bl.bases[1] = cv::Vec2d(v_mid[0], v_max[1]+eps);
        bl.bases[0] = cv::Vec2d(v_mid[0], v_mid[1]+eps);
        bl.rebuild();

        br.bases[3] = cv::Vec2d(v_mid[0]+eps, v_mid[1]+eps);
        br.bases[2] = cv::Vec2d(v_mid[0]+eps, v_max[1]+eps);
        br.bases[1] = cv::Vec2d(v_max[0]+eps, v_max[1]+eps);
        br.bases[0] = cv::Vec2d(v_max[0]+eps, v_mid[1]+eps);
        br.rebuild();

        for (size_t p=0; p < parts.size(); p++) {
            Polygon_geom np;
            
            if (tl.intersect(parts[p], np)) {
                if (!q_tl) {
                    q_tl = new Quadtree(tl);
                }
                q_tl->add_poly(np);
                printf("poly %d went to tl: %d\n", p, level);
            }
                
            if (tr.intersect(parts[p], np)) {
                if (!q_tr) {
                    q_tr = new Quadtree(tr);
                }
                q_tr->add_poly(np);
                printf("poly %d went to tr: %d\n", p, level);
            }

            if (bl.intersect(parts[p], np)) {
                if (!q_bl) {
                    q_bl = new Quadtree(bl);
                }
                q_bl->add_poly(np);
                printf("poly %d went to bl: %d\n", p, level);
            }

            if (br.intersect(parts[p], np)) {
                if (!q_br) {
                    q_br = new Quadtree(br);
                }
                q_br->add_poly(np);
                printf("poly %d went to br: %d\n", p, level);
            }
        }
        
        printf("about to recurse into child QT nodes\n");
        bool valid_children = false;
        const int max_verts_per_quad = 30;
        const int maxdepth = 6;
        if (q_tl && q_tl->total_vertices > max_verts_per_quad && level < maxdepth) {
            q_tl->partition_polygons(level+1);
            valid_children = true;
        }
        
        if (q_tr && q_tr->total_vertices > max_verts_per_quad && level < maxdepth) {
            q_tr->partition_polygons(level+1);
            valid_children = true;
        }
        
        if (q_bl && q_bl->total_vertices > max_verts_per_quad && level < maxdepth) {
            q_bl->partition_polygons(level+1);
            valid_children = true;
        }
        
        if (q_br && q_br->total_vertices > max_verts_per_quad && level < maxdepth) {
            q_br->partition_polygons(level+1);
            valid_children = true;
        }

        parts.clear();
    }

    double cx;
    double cy;
    double own_area;
      
    Quadtree* q_tl;
    Quadtree* q_tr;
    Quadtree* q_bl;
    Quadtree* q_br;
    
    bool leaf;
};

#endif // QUADTREE_H

