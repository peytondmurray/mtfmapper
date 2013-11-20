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

#include <vector>
using std::vector;

#include "geom.h"
#include "gh_clipping.h"

class Multipolygon_geom;

typedef enum {
    INSIDE = 1,
    OUTSIDE = 0,
    ON = 2
} Vertex_type;

//==============================================================================
class Polygon_geom : public Geometry {
  public:
    static const int max_verts_per_poly = 200; // max number of point in intersected polygon

    friend class Multipolygon_geom;
    
    
  
    Polygon_geom(double cx=0, double cy=0, double width=1, double height=1, double angle=0, int nverts=4) 
    : Geometry(cx, cy, 0), nvertices(nverts) {

        construct_regular_polygon(width, height, angle);
        precompute_point_test();
        own_area = compute_area();
        printf("built polygon with %d vertices, area=%lf\n\n", nverts, own_area);
    }

    Polygon_geom(const vector<cv::Vec2d>& verts) {

        nvertices = verts.size();

        bases = verts;
        normals = vector<cv::Vec2d>(nvertices);

        // compute normals
        int prev = nvertices - 1;
        for (int i=0; i < nvertices; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = (prev + 1) % nvertices;
        }
        
        precompute_point_test();

        compute_bounding_box();
        
        own_area = compute_area();
        printf("built polygon with %d vertices, area=%lf\n\n", nvertices, own_area);
    }

    virtual ~Polygon_geom(void) {
    }

    void rebuild(void) {
	// compute normals
        int prev = nvertices - 1;
        for (int i=0; i < nvertices; i++) {
            cv::Vec2d d = bases[i] - bases[prev];
            normals[i][0] = -d[1] / norm(d);
            normals[i][1] = d[0] / norm(d);
            prev = (prev + 1) % nvertices;
        }
        
        precompute_point_test();

        compute_bounding_box();
        
        own_area = compute_area();
        printf("rebuilt polygon with %d vertices, area=%lf\n\n", nvertices, own_area);
    }

    void construct_regular_polygon(double width, double height, double angle) {
        bases   = vector<cv::Vec2d>(nvertices);
        normals = vector<cv::Vec2d>(nvertices);

        printf("nvertices=%d\n", nvertices);
        printf("rendering a polygon with %d sides\n", nvertices);
        assert(nvertices >= 3);
        for (int i=0; i < nvertices; i++) {
            double phi = i*M_PI*2/double(nvertices);
            bases[i][0] = width/2*cos(angle+phi - M_PI/4.0) + cx;
            bases[i][1] = height/2*sin(angle+phi - M_PI/4.0) + cy;
            normals[i][0] = cos(angle+phi);
            normals[i][1] = -sin(angle+phi);
            //printf("%lf %lf (%lf %lf)\n", bases[i][0], bases[i][1], normals[i][0], normals[i][1]);
            if (nvertices > 4) fprintf(stderr, "%lf %lf\n", bases[i][0], bases[i][1]);
        }

        compute_bounding_box();
        
    }

    void compute_bounding_box(void) {
        // compute a bounding box
        bb_bases[0][1] = 1e12;
        bb_bases[1][0] = -1e12;
        bb_bases[2][1] = -1e12;
        bb_bases[3][0] = 1e12;

        for (int i=0; i < nvertices; i++) {
            // bb_bases: 0=top, 1=right, 2=bottom, 3=left
            bb_bases[0][1] = std::min(bases[i][1], bb_bases[0][1]);
            bb_bases[1][0] = std::max(bases[i][0], bb_bases[1][0]);
            bb_bases[2][1] = std::max(bases[i][1], bb_bases[2][1]);
            bb_bases[3][0] = std::min(bases[i][0], bb_bases[3][0]);
        }
        bb_bases[0][0] = bb_bases[3][0];
        bb_bases[1][1] = bb_bases[0][1];
        bb_bases[2][0] = bb_bases[1][0];
        bb_bases[3][1] = bb_bases[2][1];

        bb_normals[0][0] = 0;
        bb_normals[0][1] = -1;
        bb_normals[1][0] = 1;
        bb_normals[1][1] = 0;
        bb_normals[2][0] = 0;
        bb_normals[2][1] = 1;
        bb_normals[3][0] = -1;
        bb_normals[3][1] = 0;

        bb_area = compute_bb_area();

        for (int k=0; k < 4; k++) {
            printf("BB: %lf %lf (%lf %lf)\n", bb_bases[k][0], bb_bases[k][1], bb_normals[k][0], bb_normals[k][1]);
        }
    }
    
    void precompute_point_test(void) {
        constant = vector<double>(nvertices);
        multiple = vector<double>(nvertices);
        int j = (int)bases.size() - 1;

        for(int i=0; i < (int)bases.size(); i++) {
            if(bases[j][1] == bases[i][1]) {
                constant[i] = bases[i][0];
                multiple[i] = 0; 
            } else {
                constant[i] = bases[i][0] - (bases[i][1]*bases[j][0]) / (bases[j][1] - bases[i][1]) + 
                        (bases[i][1]*bases[i][0]) / (bases[j][1] - bases[i][1]);

                multiple[i] = (bases[j][0] - bases[i][0]) / (bases[j][1] - bases[i][1]); 
            }
            j = i; 
        }
    }

    // NOTE: indeterminate result if point falls exactly on boundary
    // TODO: since there is no "early out" option here, this is about
    // half the speed of the older test (see below)
    bool is_inside(double x, double y) const {

        int j = bases.size() - 1;
        bool  oddNodes = false;

        for (int i=0; i < (int)bases.size(); i++) {
            if ( (bases[i][1] < y && bases[j][1] >= y) || 
                 (bases[j][1] < y && bases[i][1] >= y) ) {

                oddNodes ^= (y*multiple[i] + constant[i]) < x; 
            }
            j = i; 
        }

        return oddNodes; 
    }

    
    bool xis_inside(double x, double y) const {
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
    
    inline double det(int i, double x, double y, bool& edge) const {
        int next = (i+1)%bases.size();
        double d = (bases[i][0] - x) * (bases[next][1] - y) -
                   (bases[next][0] - x) * (bases[i][1] - y);
                   
        if (d == 0) {
            edge = true;
        } else {
            edge = false;
        }
        return d;
    }
    
    #define crossing ( (bases[i][1] < y) != (bases[next][1] < y) )
    #define right_crossing  (det(i, x, y, edge) > 0) == (bases[next][1] > bases[i][1])  
    #define modify_w ( w += ((bases[next][1] > bases[i][1]) ? 2 : 0) - 1  )
   
    inline Vertex_type classify(double x, double y) const {
        bool edge = false;
        if (x == bases[0][0] && y == bases[0][1]) {
            return ON; // on first vertex
        }
        int w = 0;
        for (int i=0; i < (int)bases.size(); i++) {
            int next = (i+1)%bases.size();
            if (bases[next][1] == y) {
                if (bases[next][0] == x) {
                    return ON; // on a vertex
                } else {
                    if ( bases[i][1] == y && ((bases[next][0] > x) == (bases[i][0] < x)) ) {
                        return ON; // on an edge
                    }
                }
            }
            if (crossing) {
                if (bases[i][0] >= x) {
                    if (bases[next][0] > x) {
                        modify_w;
                    } else {
                        if (right_crossing) {
                            modify_w;
                        }
                        if (edge) return ON; // non-horizontal edge?
                    }
                } else {
                    if (bases[next][0] > x) {
                        if (right_crossing) {
                            modify_w;
                        }
                        if (edge) return ON; // non-horizontal edge
                    }
                }
            }
        }
        return w == 0 ? OUTSIDE : INSIDE; 
    }
    
    
    bool remove_degeneracy(int pass=0) {
        // This method is rather crude, and not tested for every possible case
    
        if (bases.size() < 3) return false;
        
        // first strip out doubles
        int pe =   (int)bases.size() - 1;
        vector<cv::Vec2d> fbases;
        for (size_t e=0; e < bases.size(); e++) {
            if ( norm(bases[pe] - bases[e]) > 1e-11 ) {
                fbases.push_back(bases[e]);
            }
            pe = e;
        }
        if (fbases.size() < bases.size()) {
            double old_area = own_area;
            bases = fbases;
            normals = fbases;
            nvertices = bases.size();
            rebuild();
            printf("after squashing duplicate points, area diff=%lf\n", old_area - own_area);
        }
        
        
        pe =   (int)bases.size() - 1;
        int pem1 = (int)bases.size() - 2;
        
        cv::Vec2d pd(bases[pe] - bases[pem1]);
        
        int collinear_count = 0;
        bool some_marked = false;
        vector<bool> marked(bases.size(), false);
        for (int e=0; e < (int)bases.size(); e++) {
            cv::Vec2d cd(bases[e] - bases[pe]);
            if ( fabs(cd[0]*pd[1] - cd[1]*pd[0]) < 2e-10 ) { // area of triangle less than 1e-10
                collinear_count++;
                
                
                printf("$$ found colinear points %lf %lf %lf %lf (pass %d)\n",
                    bases[e][0], bases[e][1],
                    bases[pe][0], bases[pe][1],
                    pass
                );
                for (size_t k=0; k < bases.size(); k++) {
                    printf("%lf %lf, ", bases[k][0], bases[k][1]);
                }
                printf("\n");
                // TODO: only remove a marked point if the area remains unaffected?
                
                printf("[%lf %lf] marked for removal\n", bases[pe][0], bases[pe][1]);
                
                marked[pe] = true;
                some_marked = true;
            } 
            
            pe = e;
            pd = cd;
        }
        
        if (some_marked) {
            double old_area = own_area;
            vector<cv::Vec2d> ncoords;
            for (size_t i=0; i < bases.size(); i++) {
                if (!marked[i]) {
                    ncoords.push_back(bases[i]);
                }
            }
            bases = ncoords;
            normals = ncoords;
            nvertices = bases.size();
            rebuild();
            printf("after nuking degenerate points, area diff=%lf\n", old_area - own_area);
        }
        
        
        return collinear_count == 0;
    }
    
    double intersection_area(const Geometry& ib, double xoffset = 0, double yoffset = 0, bool skip_bounds=false)  const {
    
        double points_x[max_verts_per_poly];
        double points_y[max_verts_per_poly];
        int points_len = 0;

        // TODO: this will probably throw an exception if you try to
        // pass a multipolygon as the photosite geometry. Don't do it!
        const Polygon_geom& b = dynamic_cast<const Polygon_geom&>(ib);

        if (!skip_bounds && b.nvertices >= 6) {
            intersect(points_x, points_y, points_len, b, xoffset, yoffset, true);
            double i_bb_area = compute_area(points_x, points_y, points_len);
            
            if (fabs(i_bb_area) < 1e-11) {
                return 0;
            } else {
                if (fabs(i_bb_area - b.bb_area) < 1e-11) {
                    return b.own_area;
                } 
            }

            // partial intersection, so reset and perform 
            // intersection again using full geometry
            points_len = 0;
        }

        intersect(points_x, points_y, points_len, b, xoffset, yoffset);
        return compute_area(points_x, points_y, points_len);
    }
    
    inline bool t_intersect(double& pix, double& piy, 
                     const double& v1x, const double& v1y,
                     const double& d1x, const double& d1y,
                     const double& v2x, const double& v2y,
                     const double& d2x, const double& d2y) const {
                   
        double denom = (d2y*d1x - d2x*d1y);
        
        if (fabs(denom) < 1e-11) {
           // this happens when the lines are parallel
           // the caller handles this correctly by not
           // adding an additional intersection point
           return false;
        }
        
        double u = (d2x*(v1y - v2y) - d2y*(v1x - v2x)) / denom;
        
        pix = v1x + u*d1x;
        piy = v1y + u*d1y;
                   
        return true;               
    }

    void print(void) const {
        for (int i=0; i < nvertices; i++) {
            printf("\t%lf %lf\n", bases[i][0], bases[i][1]);
        }
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
    
    bool compute_area_sign(void) const {
        double A = 0;
        for (int i=0; i < nvertices; i++) {
            int ni = (i+1) % nvertices;
            A += bases[i][0]*bases[ni][1] - bases[ni][0]*bases[i][1];
        }
        return A < 0;
    }

    double compute_bb_area(void) const {
        double A = 0;
        for (int i=0; i < 4; i++) {
            int ni = (i+1) % 4;
            A += bb_bases[i][0]*bb_bases[ni][1] - bb_bases[ni][0]*bb_bases[i][1];
        }
        return 0.5 * fabs(A);
    }
    
    void intersect(double* points_x, double* points_y, int& points_len, 
        const Polygon_geom& b, double xoffset = 0, double yoffset = 0, bool bounding_box=false) const {

        if (bounding_box) {
            for (int i=0; i < 4; i++) { 
                points_x[i] = b.bb_bases[i][0];
                points_y[i] = b.bb_bases[i][1];
                points_len++;
            }
        } else {
            for (int i=0; i < b.nvertices; i++) { 
                points_x[i] = b.bases[i][0];
                points_y[i] = b.bases[i][1];
                points_len++;
            }
        }

        
        for (int e=0; e < nvertices; e++) {
            intersect_core(points_x, points_y, points_len, e, nvertices, xoffset, yoffset);
        }
    }

    bool intersect(const Polygon_geom& b, Polygon_geom& s) const {

        double points_x[2*max_verts_per_poly];
	double points_y[2*max_verts_per_poly];
	int points_len = 0;

        for (int i=0; i < b.nvertices; i++) { 
            points_x[i] = b.bases[i][0];
            points_y[i] = b.bases[i][1];
            points_len++;
        }

        for (int e=0; e < nvertices; e++) {
            intersect_core(points_x, points_y, points_len, e, nvertices, 0, 0);
        }

        if (points_len > 0) {
            vector<cv::Vec2d> new_vertices(points_len);
            for (int v=0; v < points_len; v++) {
                new_vertices[v][0] = points_x[v];
                new_vertices[v][1] = points_y[v];
            }

            s = Polygon_geom(new_vertices);
            s.remove_degeneracy(0);
            //s.remove_degeneracy(1); // call it again, just to make sure. maybe something interesting is happening ?
            return s.compute_area() > 1e-11;
        } else {
            return false;
        }
        return true;
    }

    // this looks like the Sutherland-Hodgman algorithm
    // the clipping polygon must be convex (req. by SH algo)
    // will produce overlapping edges if a concave point exists
    // outside of the clipping polygon. these overlapping
    // edges cause no harm, because they have zero area (which
    // seems to work out fine with compute_area())
    void intersect_core(double* inpoints_x, double* inpoints_y, int& in_len, int e, int nedges,
                        double xoffset, double yoffset) const {
        int ne = (e + 1) % nedges;
        
        double Px = bases[e][0] + xoffset;
        double Py = bases[e][1] + yoffset;
        
        double Dx = bases[ne][0] - bases[e][0];
        double Dy = bases[ne][1] - bases[e][1];
         
        double Nx = -Dy;
        double Ny = Dx;

		double Nl = sqrt(Nx*Nx + Ny*Ny);
		Nx /= Nl;
		Ny /= Nl;
        
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
    
    
    // slower polygon clipping algorithm, but this one should handle concave-concave
    // clipping, and it should also avoid creating degenerate parts
    vector<Polygon_geom> intersect_greiner_horman(const Polygon_geom& in_b) {
        vector<Polygon_geom> polys;
        
        printf("Poly 0:\n");
        for (size_t v=0; v < bases.size(); v++) {
            printf("%lf %lf\n", bases[v][0], bases[v][1]);
        }
        printf("\nPoly 1:\n");
        double in_b_cx = 0;
        double in_b_cy = 0;
        for (size_t v=0; v < in_b.bases.size(); v++) {
            in_b_cx += in_b.bases[v][0];
            in_b_cy += in_b.bases[v][1];
            printf("%lf %lf\n", in_b.bases[v][0], in_b.bases[v][1]);
        }
        
        in_b_cx /= in_b.bases.size();
        in_b_cy /= in_b.bases.size();
        
        // now perturb all points of in_b that are classified as "ON",
        // forcing them to become "OUT"
        bool perturbed = false;
        vector<cv::Vec2d> nverts;
        double delta = 0;
        for (size_t v=0; v < in_b.bases.size(); v++) {
            double px = in_b.bases[v][0];
            double py = in_b.bases[v][1];
            if (classify(px, py) == ON) {
                //printf("point on edge, perturbing\n");
                perturbed = true;
                double dx = px - in_b_cx;
                double dy = py - in_b_cy;
                
                double n = sqrt(dx*dx + dy*dy);
                
                if (n < 1e-11) {
                    printf("reverting to random perturbation\n");
                    px += 1e-10*rand()/double(RAND_MAX);
                    py += 1e-10*rand()/double(RAND_MAX);
                } else {
                    
                    dx /= n;
                    dy /= n;
                    
                    if (fabs(dx) < 1e-11) {
                        dx = dx < 0 ? -1 : 1;
                    }
                    if (fabs(dy) < 1e-11) {
                        dy = dy < 0 ? -1 : 1;
                    }
                    
                    dx *= 1e-4;
                    dy *= 1e-4;
                    
                    if (sqrt(dx*dx+dy*dy) > delta) {
                        delta = sqrt(dx*dx+dy*dy);
                    }
                    
                    if (classify(px + dx, py + dy) == INSIDE) {
                        px -= dx;
                        py -= dy;
                    } else {
                        px += dx;
                        py += dy;
                    }
                    
                }
                
                if (classify(px, py) != OUTSIDE) {
                    printf("## perturbation failed to move point to outside?\n");
                }
            }
            nverts.push_back(cv::Vec2d(px, py));
        }
        
        Polygon_geom b(nverts);
        
        vector<GH_clipping::gh_vertex> verts(nvertices * b.nvertices*2 + nvertices + b.nvertices);
        vector<Polygon_geom> poly(2);
        poly[0] = *this;
        poly[1] = b;
        
        // populate the verts vector with the two polys
        int poly1_start = GH_clipping::init_gh_list(verts, bases, 0, 1);
        int vs = GH_clipping::init_gh_list(verts, b.bases, poly1_start, -1);
        
        int vs_before_intersections = vs;
        GH_clipping::gh_phase_one(verts, vs, bases.size(), b.bases.size());
        
        if (vs == vs_before_intersections) {
        
            bool all_on = true;
            for (size_t p=0; p < bases.size(); p++) {
                int cl = b.classify(bases[p][0], bases[p][1]);
                if (cl == OUTSIDE) {
                    all_on = false;
                }
            }
            
            if (all_on) {
                // *this must be entirely within b, so return *this
                printf("S inside C: returning S\n");
                polys.push_back(*this);
                return polys;
            } else {
                // maybe b is entirely inside *this?
                bool all_in = true;
                for (size_t p=0; p < b.bases.size(); p++) {
                    int cl = classify(b.bases[p][0], b.bases[p][1]);
                    if (cl == OUTSIDE) {
                        all_in = false;
                    }
                }
                
                if (all_in) {
                    printf("C inside S: returning C\n");
                    polys.push_back(b);
                    return polys;
                }
            
                // *this is entirely outside b, so return empty list
                printf("no intersection\n");
                return polys;
            }
        }
        
        GH_clipping::gh_phase_two(verts, b, 0);
        GH_clipping::gh_phase_two(verts, *this, poly1_start);
        
        for (size_t v=0; v < vs; v++) {
            printf("%d: (%lf, %lf), n=%d, p=%d, en=%d, isect=%d, neigh=%d\n",
                v, verts[v].x, verts[v].y,
                verts[v].next, verts[v].prev,
                verts[v].en, verts[v].isect,
                verts[v].neighbour
            );
        }
        
        printf("Poly 0\n");
        int cur = 0;
        do {
            int& v = cur;
            printf("%d: (%lf, %lf), n=%d, p=%d, en=%d, isect=%d, neigh=%d\n",
                v, verts[v].x, verts[v].y,
                verts[v].next, verts[v].prev,
                verts[v].en, verts[v].isect,
                verts[v].neighbour
            );
            cur = verts[v].next;
        } while (cur != 0);
        
        printf("Poly 1\n");
        cur = poly1_start;
        do {
            int& v = cur;
            printf("%d: (%lf, %lf), n=%d, p=%d, en=%d, isect=%d, neigh=%d\n",
                v, verts[v].x, verts[v].y,
                verts[v].next, verts[v].prev,
                verts[v].en, verts[v].isect,
                verts[v].neighbour
            );
            cur = verts[v].next;
        } while (cur != poly1_start);
        
        GH_clipping::gh_phase_three(verts, vs, vs_before_intersections, polys);
        
        if (perturbed) { // clean up the near-duplicates introduced by perturbation
            vector<Polygon_geom> npolys;
            for (size_t k=0; k < polys.size(); k++){
                
                vector<cv::Vec2d> nverts;
                
                int prev = polys[k].bases.size() - 1;
                for (int j=0; j < polys[k].bases.size(); j++) {
                
                    if (norm(polys[k].bases[j] - polys[k].bases[prev]) < 2.5*delta) {
                    } else {
                        nverts.push_back(polys[k].bases[j]);
                    }
                    prev = j;
                }
                if (nverts.size() > 2) {
                    npolys.push_back(Polygon_geom(nverts));
                }
            }    
            polys = npolys;
        }
        
        //for (size_t k=0; k < polys.size(); k++){
        //    printf("poly %d\n", (int)k);
        //    for (int j=0; j < polys[k].bases.size(); j++) {
        //        printf("\t%lf %lf\n", polys[k].bases[j][0], polys[k].bases[j][1]);
        //    }
        //}
        
        return polys;
    }

    vector<cv::Vec2d> normals;
    vector<cv::Vec2d> bases;

    cv::Vec2d bb_normals[4];
    cv::Vec2d bb_bases[4];

    vector<double> constant;
    vector<double> multiple;
    
    int nvertices;

    double bb_area;
};

#endif // RENDER_H
