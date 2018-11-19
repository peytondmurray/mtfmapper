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
#ifndef MRECTANGLE_H
#define MRECTANGLE_H

#include "common_types.h"
#include "peak_detector.h"
#include "edge_record.h"

const double rect_il_thresh = 3.75;
const double adjust = 0.15;

class Mrectangle {
  public:
  
    Mrectangle(void) { }
    
    // build a rectangular buffer of "width" around midpoint of side k
    Mrectangle(const Mrectangle& b, size_t k, double width) 
      : thetas(4, .0), centroids(4, Point2d(0.0,0.0)), valid(true), 
        corners(4, Point2d(0.0,0.0)), edges(4, Point2d(0.0,0.0)), 
        normals(4, Point2d(0.0,0.0)), boundary_length(0) {
        
        assert(k < b.centroids.size());
        assert(b.corners.size() == 4);
        assert(b.corner_map[k].size() == 2);

        const Point2d& n  = b.normals[k];
        Point2d c1 = b.corners[b.corner_map[k][0]];
        Point2d c2 = b.corners[b.corner_map[k][1]];
        
        Point2d diff1 = Point2d(c1.x - c2.x, c1.y - c2.y);
        double len = sqrt(diff1.ddot(diff1));
        Point2d pn(-n.y, n.x);
        double e_adj = min(len*adjust, 4.0)/len;
        e_adj = max(e_adj, 0.08); // TODO: this must probably be larger on wider PSFs
        Point2d delta(-pn.x*len*e_adj, -pn.y*len*e_adj);
        
        if (ndiff(avg(c1,c2), c1).ddot(pn) > 0) {
            delta.x = -delta.x;
            delta.y = -delta.y;
        }
        
        c1.x += delta.x;
        c1.y += delta.y;
        c2.x -= delta.x;
        c2.y -= delta.y;

        corners[0] = Point2d(c1.x + width*n.x, c1.y + width*n.y);
        corners[1] = Point2d(c1.x - width*n.x, c1.y - width*n.y);
        corners[2] = Point2d(c2.x + width*n.x, c2.y + width*n.y);
        corners[3] = Point2d(c2.x - width*n.x, c2.y - width*n.y);

        normals[3] = ndiff(corners[2], corners[0]);
        normals[2] = ndiff(corners[1], corners[3]);
        normals[0] = ndiff(corners[0], corners[1]);
        normals[1] = ndiff(corners[3], corners[2]);
        
        centroids[0] = avg(corners[2], corners[0]);
        centroids[1] = avg(corners[1], corners[3]);
        centroids[2] = avg(corners[0], corners[1]);
        centroids[3] = avg(corners[2], corners[3]);
        
        tl.x = 1e50;
        br.x = -1e50;
        tl.y = 1e50;
        br.y = -1e50;
        for (size_t k=0; k < 4; k++) {
            Point2d& c = corners[k];
            if (c.x < tl.x) tl.x = c.x;
            if (c.x > br.x) br.x = c.x;
            if (c.y < tl.y) tl.y = c.y;
            if (c.y > br.y) br.y = c.y;
        }
        
        tl.x = floor(tl.x);
        br.x = ceil(br.x);
        tl.y = floor(tl.y);
        br.y = ceil(br.y);
    }
    
    // reposition a rectangle using new estimates of the centroids and normals
    Mrectangle(const Mrectangle& b, const vector<Edge_record>& edge_records) 
      : thetas(4, 0.0), centroids(4, Point2d(0.0,0.0)), valid(true), 
        corners(4, Point2d(0.0,0.0)), edges(4, Point2d(0.0,0.0)), 
        normals(4, Point2d(0.0,0.0)), boundary_length(0) {
        
        Point2d sq_centre(0,0);
        for (int k=0; k < 4; k++) {
            sq_centre += edge_records[k].centroid;
        }
        sq_centre *= 0.25;
    
        for (int k=0; k < 4; k++) {
            centroids[k] = edge_records[k].centroid;
            normals[k] = Point2d(cos(edge_records[k].angle), sin(edge_records[k].angle)); 
            Point2d delta = centroids[k] - sq_centre;
            delta *= 1.0/(norm(delta));
            double dot = normals[k].x*delta.x + normals[k].y*delta.y;
            if (dot < 0) {
                normals[k] = -normals[k];
            }
            edges[k].x = normals[k].y;
            edges[k].y = normals[k].x;
            thetas[k] = atan2(normals[k].y, normals[k].x);
        }
        
        boundary_length = b.boundary_length;
        corner_map = vector< vector<int> >(4);
        
        int corner_idx = 0;
        for (size_t k1=0; k1 < 3 && corner_idx < 4; k1++) {
            for (size_t k2=k1+1; k2 < 4 && corner_idx < 4; k2++) {
                // check if edges are approximately normal
                double dot = normals[k1].x*normals[k2].x + normals[k1].y*normals[k2].y;
                if ( fabs(fabs(dot)) < M_PI/6.0 ) {
                    // edges are approximately normal, find intersection
                    Point2d isect(0.0,0.0);
                    intersect(centroids[k1], normals[k1], centroids[k2], normals[k2], isect);
                    corner_map[k1].push_back(corner_idx);
                    corner_map[k2].push_back(corner_idx);
                    corners[corner_idx++] = isect;
                }
            }
        }
        
        tl.x = 1e50;
        br.x = -1e50;
        tl.y = 1e50;
        br.y = -1e50;
        for (size_t k=0; k < 4; k++) {
            Point2d& c = corners[k];
            if (c.x < tl.x) tl.x = c.x;
            if (c.x > br.x) br.x = c.x;
            if (c.y < tl.y) tl.y = c.y;
            if (c.y > br.y) br.y = c.y;
        }
        
        tl.x = floor(tl.x);
        br.x = ceil(br.x);
        tl.y = floor(tl.y);
        br.y = ceil(br.y);
    }
    
    bool corners_ok(void) const {
        bool ok = true;
        if (corner_map.size() != 4) {
            return false;
        }
        for (int k=0; k < 4; k++) {
            ok &= corner_map[k].size() == 2;
        }
        return ok;
    }
  
    Point2d ndiff(const Point2d& a, const Point2d& b) {
        Point2d diff(a.x - b.x, a.y - b.y);
        double l = sqrt(diff.ddot(diff));
        diff.x /= l;
        diff.y /= l;
        return diff;
    }
    
    Point2d avg(const Point2d& a, const Point2d& b) {
        return Point2d( (a.x+b.x)/2.0, (a.y+b.y)/2.0 );
    }
  
    Mrectangle(const vector<double>& in_thetas, const vector<double>& data_thetas, 
        const vector<Point2d>& points, const Gradient& g, double thresh=5.0/180.0*M_PI) 
      : thetas(in_thetas), centroids(4, Point2d(0.0,0.0)), valid(false), 
        corners(4, Point2d(0.0,0.0)), edges(4, Point2d(0.0,0.0)), 
        normals(4, Point2d(0.0,0.0)), corner_map(4) {
      
        if (thetas.size() == 4) {
    
            vector<double> c_weight(4, 0);
            for (size_t i=0; i < points.size(); i++) { 
                for (size_t k=0; k < 4; k++) {
                    if (Peak_detector::angular_diff(data_thetas[i], thetas[k]) < thresh) {
                        double w = g.grad_magnitude(points[i].x, points[i].y);
                        centroids[k].x += points[i].x * w;
                        centroids[k].y += points[i].y * w;
                        c_weight[k] += w;
                    } 
                }
            }
            for (size_t k=0; k < 4; k++) {
                centroids[k].x /= c_weight[k];
                centroids[k].y /= c_weight[k];
            }
            
            for (size_t k=0; k < 4; k++) {
                edges[k].x = sin(thetas[k]);
                edges[k].y = cos(thetas[k]);
                normals[k].x = cos(thetas[k]);
                normals[k].y = sin(thetas[k]);
            }
            
            // now we have the normal to each edge, as well as a point on each line
            // TODO: sanity checks to see if this is a rectangle
            
            // check if most points are closer than dthresh pixels from any line
            size_t outlier_count = 0;
            for (size_t i=0; i < points.size(); i++) {
                double min_dist = 1e50;
                for (size_t k=0; k < 4; k++) {
                    Point2d v(centroids[k].x - points[i].x, centroids[k].y - points[i].y);
                    double dist = fabs(v.ddot(normals[k]));
                    min_dist = min(min_dist, dist);
                }
                if (min_dist > rect_il_thresh) {
                    outlier_count++;
                }
            }
            size_t corner_idx = 0;
            if (outlier_count > points.size()/8) {
                valid = false;
            } else {
                valid = true;
                boundary_length = points.size();
                
                for (size_t k1=0; k1 < 3 && corner_idx < 4; k1++) {
                    for (size_t k2=k1+1; k2 < 4 && corner_idx < 4; k2++) {
                        // check if edges are approximately normal
                        
                        double dot = acos(edges[k1].ddot(edges[k2]));
                        if ( fabs(fabs(dot) - M_PI/2.0) < M_PI/5.1429 ) { // about 35 degrees
                            // edges are approximately normal, find intersection
                            Point2d isect(0.0,0.0);
                            intersect(centroids[k1], normals[k1], centroids[k2], normals[k2], isect);
                            
                            corner_map[k1].push_back(corner_idx);
                            corner_map[k2].push_back(corner_idx);
                            corners[corner_idx++] = isect;
                        } 
                    }
                }
            }
            if (corner_idx < 4) {
                valid = false;
            }
        } else {
            valid = false;
        }
    }
    
    bool intersect(const Point2d& p1, const Point2d& d1, const Point2d& p2, 
        const Point2d& d2, Point2d& isect) {
        
        double dk1 = p1.ddot(d1);
        double dk2 = p2.ddot(d2);
        
        double a1 = d1.x;
        double a2 = d2.x;
        double b1 = d1.y;
        double b2 = d2.y;
        
        double det = (a1*b2 - a2*b1);
        
        if ( fabs(det) < 1e-12 ) {
            // lines are actually parallel. this is impossible?
            printf("Warning: determinant near-zero: %s, %d\n", __FILE__, __LINE__);
            return false;
        }
        
        isect.x = (b2*dk1 - b1*dk2)/det;
        isect.y = (a1*dk2 - a2*dk1)/det;
        return true;
    }
    
    bool is_inside(const Point2d& p) const {
        // a point is inside the rectangle if it falls along the negative direction of each normal
        
        for (size_t k=0; k < 4; k++) {
            Point2d ldir(p.x - centroids[k].x, p.y - centroids[k].y);
            double dot = ldir.ddot(normals[k]);
            if (dot > 0) {
                return false;
            }
        }
        return true;
    }
    
    Point2d get_centroid(size_t i) const {
        return centroids[i];
    }
    
    void print(void) const {
        for (int k=0; k < 4; k++) {
            printf("c(%lf, %lf), e(%lf, %lf), n(%lf, %lf), t(%lf), cr(%lf, %lf), map(%d, %d)\n",
                centroids[k].x, centroids[k].y,
                edges[k].x, edges[k].y,
                normals[k].x, normals[k].x,
                thetas[k],
                corners[k].x, corners[k].y,
                corner_map[k][0], corner_map[k][1]
            );
            
        }
    }
    
    vector<double> thetas;
    vector<Point2d>  centroids;  
    bool           valid;
    vector<Point2d>  corners;
    vector<Point2d>  edges;
    vector<Point2d>  normals;
    vector< vector<int> > corner_map;
    Point2d          tl;
    Point2d          br;
    size_t         boundary_length;
};

#endif
