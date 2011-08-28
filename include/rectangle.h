#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "common_types.h"
#include "peak_detector.h"

const double rect_il_thresh = 2.0;
const double adjust = 0.15;

class Rectangle {
  public:
  
    Rectangle(void) { }
    
    // build a rectangular buffer of "width" around midpoint of side k
    Rectangle(const Rectangle& b, size_t k, double width) 
      : thetas(4, .0), centroids(4, Point(0.0,0.0)), valid(true), 
        corners(4, Point(0.0,0.0)), edges(4, Point(0.0,0.0)), 
        normals(4, Point(0.0,0.0)) {
        
        assert(k < b.centroids.size());
        assert(b.corners.size() == 4);
        assert(b.corner_map[k].size() == 2);
        
        
        const Point& n  = b.normals[k];
        Point c1 = b.corners[b.corner_map[k][0]];
        Point c2 = b.corners[b.corner_map[k][1]];
        
        Point diff1 = Point(c1.x - c2.x, c1.y - c2.y);
        double len = sqrt(diff1.ddot(diff1));
        Point pn(-n.y, n.x);
        Point delta(-pn.x*len*adjust, -pn.y*len*adjust);
        
        if (ndiff(avg(c1,c2), c1).ddot(pn) > 0) {
            delta.x = -delta.x;
            delta.y = -delta.y;
        }
        
        c1.x += delta.x;
        c1.y += delta.y;
        c2.x -= delta.x;
        c2.y -= delta.y;
        
        
        corners[0] = Point(c1.x + width*n.x, c1.y + width*n.y);
        corners[1] = Point(c1.x - width*n.x, c1.y - width*n.y);
        corners[2] = Point(c2.x + width*n.x, c2.y + width*n.y);
        corners[3] = Point(c2.x - width*n.x, c2.y - width*n.y);
        
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
            Point& c = corners[k];
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
  
    Point ndiff(const Point& a, const Point& b) {
        Point diff(a.x - b.x, a.y - b.y);
        double l = sqrt(diff.ddot(diff));
        diff.x /= l;
        diff.y /= l;
        return diff;
    }
    
    Point avg(const Point& a, const Point& b) {
        return Point( (a.x+b.x)/2.0, (a.y+b.y)/2.0 );
    }
  
    Rectangle(const vector<double>& in_thetas, const vector<double>& data_thetas, 
        const vector<Point>& points, double thresh=5.0/180.0*M_PI) 
      : thetas(in_thetas), centroids(4, Point(0.0,0.0)), valid(false), 
        corners(4, Point(0.0,0.0)), edges(4, Point(0.0,0.0)), 
        normals(4, Point(0.0,0.0)), corner_map(4) {
      
        if (thetas.size() == 4) {
    
            vector<size_t> c_count(4);
            for (size_t i=0; i < points.size(); i++) { 
                for (size_t k=0; k < 4; k++) {
                    if (Peak_detector::angular_diff(data_thetas[i], thetas[k]) < thresh) {
                        c_count[k]++;
                        centroids[k].x += points[i].x;
                        centroids[k].y += points[i].y;
                    } 
                }
            }
            for (size_t k=0; k < 4; k++) {
                centroids[k].x /= double(c_count[k]);
                centroids[k].y /= double(c_count[k]);
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
                    Point v(centroids[k].x - points[i].x, centroids[k].y - points[i].y);
                    double dist = fabs(v.ddot(normals[k]));
                    min_dist = std::min(min_dist, dist);
                }
                if (min_dist > rect_il_thresh) {
                    outlier_count++;
                }
            }
            size_t corner_idx = 0;
            if (outlier_count > points.size()/10) {
                valid = false;
            } else {
                valid = true;
                boundary_length = points.size();
                
                for (size_t k1=0; k1 < 3; k1++) {
                    for (size_t k2=k1+1; k2 < 4; k2++) {
                        // check if edges are approximately normal
                        
                        double dot = acos(edges[k1].ddot(edges[k2]));
                        if ( fabs(fabs(dot) - M_PI/2.0) < M_PI/8.0 ) {
                            // edges are approximately normal, find intersection
                            Point isect(0.0,0.0);
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
    
    bool intersect(const Point& p1, const Point& d1, const Point& p2, 
        const Point& d2, Point& isect) {
        
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
    
    bool is_inside(const Point& p) const {
        // a point is inside the rectangle if it falls along the negative direction of each normal
        
        for (size_t k=0; k < 4; k++) {
            Point ldir(p.x - centroids[k].x, p.y - centroids[k].y);
            double dot = ldir.ddot(normals[k]);
            if (dot > 0) {
                return false;
            }
        }
        return true;
    }
    
    Point get_centroid(size_t i) const {
        return centroids[i];
    }
    
    vector<double> thetas;
    vector<Point>  centroids;  
    bool           valid;
    vector<Point>  corners;
    vector<Point>  edges;
    vector<Point>  normals;
    vector< vector<int> > corner_map;
    Point          tl;
    Point          br;
    size_t         boundary_length;
};

#endif
