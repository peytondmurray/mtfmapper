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

class Multipolygon_geom;

struct gh_vertex {
    double x;
    double y;
    int next;       // we use integers here to avoid multiple dynamic allocations
    int prev;       // and the corresponding pain of freeing up the linked list again
    bool isect;
    bool en;
    double alpha;
    int neighbour;
    int next_poly;  // we also use a global vector for storing the polygons ... but there are only two?
};

//==============================================================================
class Polygon_geom : public Geometry {
  public:
    static const int max_verts_per_poly = 200; // max number of point in intersected polygon

    friend class Multipolygon_geom;
  
    Polygon_geom(double cx=0, double cy=0, double width=1, double height=1, double angle=0, int nverts=4) 
    : Geometry(cx, cy, 0), nvertices(nverts) {

        construct_regular_polygon(width, height, angle);
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

    bool is_inside(double x, double y) const {
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
    
                      
    inline bool t_intersect(const gh_vertex& s0, const gh_vertex& s1,
                            const gh_vertex& c0, const gh_vertex& c1,
                            gh_vertex& is, gh_vertex& ic) const {                      
                      
        double ds_x = s1.x - s0.x;
        double ds_y = s1.y - s0.y;
        double dc_x = c1.x - c0.x;
        double dc_y = c1.y - c0.y;
                                          
        double denom = (dc_y*ds_x - dc_x*ds_y);
                                                          
        if (fabs(denom) < 1e-11) {
            printf("denom zero? lines are parallel?\n");
            //printf("1: (%lf, %lf), (%lf, %lf)\n", v1[0], v1[1], d1[0], d1[1]);
            //printf("2: (%lf, %lf), (%lf, %lf)\n", v2[0], v2[1], d2[0], d2[1]);
            return false;
        }
                                                                                                                              
        is.alpha = (dc_x*(s0.y - c0.y) - dc_y*(s0.x - c0.x)) / denom;
        ic.alpha = -(ds_x*(c0.y - s0.y) - ds_y*(c0.x - s0.x)) / denom;
        
        printf("S: (%lf, %lf), d=(%lf, %lf)\n", s0.x, s0.y, ds_x, ds_y);
        printf("C: (%lf, %lf), d=(%lf, %lf)\n", c0.x, c0.y, dc_x, dc_y);
        
        printf("is.alpha=%lf, ic.alpha=%lf\n", is.alpha, ic.alpha);
        
        // should we check that the alphas are in the range [0,1] ?
        
        if (is.alpha < 0 || is.alpha > 1 ||
            ic.alpha < 0 || ic.alpha > 1) {
            
            printf("intersection outside of end vertices, not reporting this as an intersection\n");
            return false;
        }
        
        is.x = s0.x + is.alpha*ds_x;
        is.y = s0.y + is.alpha*ds_y;
        is.isect = true;
        
        ic.x = c0.x + ic.alpha*dc_x;
        ic.y = c0.y + ic.alpha*dc_y;
        ic.isect = true;
        
        return true;        
    }
    
    //double x;
    //double y;
    //int next;       // we use integers here to avoid multiple dynamic allocations
    //int prev;       // and the corresponding pain of freeing up the linked list again
    //bool isect;
    //bool en;
    //int neighbour;
    //int next_poly;  // we also use a global vector for storing the polygons ... but there are only two?
    
    int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly) {
        for (int v=0; v < (int)in_verts.size(); v++) {
            verts[v+vs].x = in_verts[v][0];// + (v+vs+1)/1e5;
            verts[v+vs].y = in_verts[v][1];// - (v+vs+1)/2e5;
            verts[v+vs].next = (v+1) % in_verts.size() + vs;
            verts[v+vs].prev = (v + in_verts.size() - 1) % in_verts.size() + vs;
            verts[v+vs].isect = false;  // this it not an intersection, so it is a vertex
            verts[v+vs].en = false;     // thus, en is the inside/outside flag???
            verts[v+vs].neighbour = -1; // no matching vertex in neighbour yet
            verts[v+vs].next_poly = next_poly; // no next poly yet
            verts[v+vs].alpha = 10;     // any value > 1 indicates no intersection, which is the default
            
            printf("v=%d, v+vs=%d, next=%d, prev=%d\n", v, vs+v, verts[v+vs].next, verts[v+vs].prev);
        }
        return vs + in_verts.size(); // index of next free vertex entry
    }
    
    void gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size) {
        // phase one
        for (int si=0; si < poly0_size; si++) {
            int nsi = (si + 1) % poly0_size;
            for (int ci=0; ci < poly1_size; ci++) {
                int nci = (ci + 1) % poly1_size;
                
                printf("testing %d from p0, %d from p1\n", si, ci);
                
                bool isect = t_intersect(
                    verts[si], verts[nsi], 
                    verts[ci+poly0_size], verts[nci+poly0_size], 
                    verts[vs], verts[vs+1]
                );
                
                printf("result: %d (%lf, %lf), (%lf, %lf)\n", isect, verts[vs].x, verts[vs].y, verts[vs+1].x, verts[vs+1].y);
                    
                if (isect) { 
                    verts[vs].neighbour = vs+1;
                    verts[vs].next_poly = 1;
                    verts[vs+1].neighbour = vs;
                    verts[vs+1].next_poly = -1;
                    
                    // insert si into chain for poly 0
                    
                    // find the node that points to nsi
                    int current = si;
                    while (verts[current].next != nsi ) { //&& // this test may be redundant
                           //verts[vs].alpha > verts[current].alpha) {
                        current = verts[current].next;
                    }
                    // current.next == nsi
                    
                    verts[vs].next = nsi;
                    verts[vs].prev = current;
                    verts[nsi].prev = vs;
                    verts[current].next = vs;
                    printf("setting prev of %d to %d\n", nsi, vs);
                    
                    // same for ci into poly 1
                    
                    // find the node that points to nci
                    current = ci + poly0_size;
                    while (verts[current].next != (nci+poly0_size) ) { //&& // this test may be redundant
                           //verts[vs+1].alpha > verts[current].alpha) { // why is this not working ?
                        current = verts[current].next;
                    }
                    // current.next == nci+poly0_size
                    
                    verts[vs+1].next = nci + poly0_size;
                    verts[vs+1].prev = current;
                    verts[nci + poly0_size].prev = vs + 1;
                    verts[current].next = vs + 1;
                    printf("**setting prev of %d to %d\n", nci+poly0_size, vs+1);
                    
                    vs += 2;
                }
            }
        }
    }
    
    void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index) {
        bool status = !b.is_inside(verts[first_vert_index].x, verts[first_vert_index].y);
        
        printf("first vertex inside other poly: %d\n", status);
        
        // inside == true -> status == exit (true?)
        int current = first_vert_index;
        do {
            printf("current vert is %d\n", current);
            if (verts[current].isect) {
                printf("intersection (%lf, %lf): en=%d (0=exit, 1=entry)\n", verts[current].x, verts[current].y, status);
                verts[current].en = status;
                status = !status;
            }
            
            current = verts[current].next;
        } while (current != first_vert_index); // let us hope the chain is not broken ...
    }
    
    void gh_phase_three(vector<gh_vertex>& verts, int vs, int first_isect_index, vector<Polygon_geom>& polys) {
        // any vertex with a next_poly value of 1 is part of the subject polygon
        // we can set the next field to -2 when we have processed a node
        bool done = false;
        while (!done) {
            vector<cv::Vec2d> vertices;
            // scan through list to find inprocessed vertices belonging to subject polygon
            done = true;
            for (int si=first_isect_index; si < vs; si++) {
                if (verts[si].next_poly == 1 && verts[si].isect) {
                    printf("found unprocessed intersection %d\n", si);
                    done = false;
                    
                    // in the GH paper, the current vertex is added here
                    // but this produces duplication (first==last), so we just skip it here ...
                    int current = si;
                    do {
                        verts[current].next_poly = -2;
                        if (verts[current].en) {
                            do {
                                current = verts[current].next;
                                vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                                verts[current].next_poly = -2;
                            } while (!verts[current].isect);
                        } else {
                            do {
                                current = verts[current].prev;
                                vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                                verts[current].next_poly = -2;
                            } while(!verts[current].isect);
                        }
                        current = verts[current].neighbour; // swap to the other poly
                    } while (current != si); // TODO: will this work if the last intersection was in clip poly?
                    if (vertices.size() > 0) {
                        printf("adding a poly with %d vertices\n", (int)vertices.size());
                        for (size_t k=0; k < vertices.size(); k++) {
                            printf("\t%lf %lf\n", vertices[k][0], vertices[k][1]);
                        }
                        polys.push_back(Polygon_geom(vertices));
                        vertices.clear();
                    }
                }
            }
            
        }
    }
    
    // slower polygon clipping algorithm, but this one should handle concave-concave
    // clipping, and it should also avoid creating degenerate parts
    vector<Polygon_geom> intersect_greiner_horman(const Polygon_geom& b) {
        vector<Polygon_geom> polys;
    
        vector<gh_vertex> verts(nvertices * b.nvertices*2 + nvertices + b.nvertices);
        vector<Polygon_geom> poly(2);
        poly[0] = *this;
        poly[1] = b;
        
        // populate the verts vector with the two polys
        int poly1_start = init_gh_list(verts, bases, 0, 1);
        printf("next!\n");
        int vs = init_gh_list(verts, b.bases, poly1_start, -1);
        
        printf("vs is now = %d\n", vs);
        for (int i=0; i < vs; i++) {
            printf("i=%d: (%lf %lf), p=%d, n=%d\n", i, verts[i].x, verts[i].y, verts[i].prev, verts[i].next);
        }
        
        int vs_before_intersections = vs;
        gh_phase_one(verts, vs, bases.size(), b.bases.size());
        
        printf("after phase 1, vs is %d\n", vs);
        
        for (int i=0; i < vs; i++) {
            printf("i=%d: (%lf %lf), p=%d, n=%d, np=%d\n", i, verts[i].x, verts[i].y, verts[i].prev, verts[i].next, verts[i].next_poly);
        }
        
        if (vs == vs_before_intersections) {
            // either *this is entirely inside b, or the other way round
            // check a single vertex of *this to decide
            printf("** No intersections found ...\n");
            
            if (b.is_inside(bases[0][0], bases[0][1])) {
                // *this must be entirely within b, so return *this
                polys.push_back(*this);
                return polys;
            } else {
                // *this is entirely outside b, so return empty list
                return polys;
            }
        }
        
        
        // phase 2 (of original gh algo)
        printf("phase 2, poly0\n");
        gh_phase_two(verts, b, 0);
        printf("phase 2, poly1\n");
        gh_phase_two(verts, *this, poly1_start);
        
        // phase 3 ....
        gh_phase_three(verts, vs, vs_before_intersections, polys);
        
        return polys;
    }

    vector<cv::Vec2d> normals;
    vector<cv::Vec2d> bases;

    cv::Vec2d bb_normals[4];
    cv::Vec2d bb_bases[4];
    
    int nvertices;

    double bb_area;
};

#endif // RENDER_H
