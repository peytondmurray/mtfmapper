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

#include "gh_clipping.h"
#include "polygon_geom.h"

namespace GH_clipping {

static inline bool t_intersect(const gh_vertex& s0, const gh_vertex& s1,
                        const gh_vertex& c0, const gh_vertex& c1,
                        gh_vertex& is, gh_vertex& ic) {
                  
    double ds_x = s1.x - s0.x;
    double ds_y = s1.y - s0.y;
    double dc_x = c1.x - c0.x;
    double dc_y = c1.y - c0.y;
                                      
    double denom = (dc_y*ds_x - dc_x*ds_y);
                                                      
    if (fabs(denom) < 1e-11) {
        return false;
    }
                                                                                                                          
    is.alpha = (dc_x*(s0.y - c0.y) - dc_y*(s0.x - c0.x)) / denom;
    ic.alpha = -(ds_x*(c0.y - s0.y) - ds_y*(c0.x - s0.x)) / denom;
    
    
    // should we check that the alphas are in the range [0,1] ?
    if (is.alpha <= 0 || is.alpha > 1 ||  
        ic.alpha <= 0 || ic.alpha > 1) {
    
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


static void insert_vertex(vector<gh_vertex>& verts, int v1, int next, int vs) {
    int current = v1;
        
    if (verts[current].next != next) {
        do { 
            if (verts[vs].alpha > verts[verts[current].next].alpha) {
                current = verts[current].next;
            } else {
                break;
            }
        } while (verts[vs].alpha > verts[verts[current].next].alpha);
    }
    
    verts[vs].next = verts[current].next;
    verts[vs].prev = current;
    verts[verts[current].next].prev = vs;
    verts[current].next = vs;
    verts[vs].en = true;
}

static inline double dist(const gh_vertex& a, const gh_vertex& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}


int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly) {
    for (int v=0; v < (int)in_verts.size(); v++) {
        verts[v+vs].x = in_verts[v][0];
        verts[v+vs].y = in_verts[v][1];
        verts[v+vs].next = (v+1) % in_verts.size() + vs;
        verts[v+vs].prev = (v + in_verts.size() - 1) % in_verts.size() + vs;
        verts[v+vs].isect = false;  // this it not an intersection, so it is a vertex
        verts[v+vs].en = false;     // thus, en is the inside/outside flag???
        verts[v+vs].neighbour = -1; // no matching vertex in neighbour yet
        verts[v+vs].next_poly = next_poly; // no next poly yet
        verts[v+vs].alpha = 10;     // any value > 1 indicates no intersection, which is the default
    }
    return vs + in_verts.size(); // index of next free vertex entry
}



void gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size) {
    // phase one
    for (int si=0; si < poly0_size; si++) {
        int nsi = (si + 1) % poly0_size;
        for (int ci=0; ci < poly1_size; ci++) {
            int nci = (ci + 1) % poly1_size;
            
            bool isect = t_intersect(
                verts[si], verts[nsi], 
                verts[ci+poly0_size], verts[nci+poly0_size], 
                verts[vs], verts[vs+1]
            );
            
            if (isect) { 
                verts[vs].neighbour = vs+1;
                verts[vs].next_poly = 1;
                verts[vs+1].neighbour = vs;
                verts[vs+1].next_poly = -7;
                
                
                if (isect) {
                  if (dist(verts[vs], verts[si]) < 1e-11 ) {
                    printf("##new vertex similar to existing point in poly 0 : vs=%d\n", vs);
                    insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                    verts[vs+1].neighbour = si;
                    verts[si].neighbour = vs+1;
                    verts[si].next_poly = -10;
                    verts[vs+1].isect = 1;
                    verts[si].isect = 1;
                    verts[vs].isect = 0;
                    vs += 2;
                  } else {
                      if (dist(verts[vs], verts[ci+poly0_size]) < 1e-11) {
                          printf("##new vertex similar to existing point in poly 1\n");
                          insert_vertex(verts, si, nsi, vs);
                          verts[vs].neighbour = ci+poly0_size;
                          verts[ci+poly0_size].neighbour = vs;
                          verts[ci+poly0_size].next_poly = -10;
                          verts[vs].isect = 1;
                          verts[ci+poly0_size].isect = 1;
                          verts[vs+1].isect = 0;
                          vs += 2;
                      } else {
                          insert_vertex(verts, si, nsi, vs);
                          insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                          vs += 2;
                      }
                  }
                  
                  
                }
                
            }
        }
    }
}


void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index) {
    bool status = b.classify(verts[first_vert_index].x, verts[first_vert_index].y) == INSIDE ? false : true;
    
    
    // inside == true -> status == exit (true?)
    int current = first_vert_index;
    do {
        if (verts[current].next_poly == -10) {
            printf("\tflipping status of vertex %d\n", current);
            verts[current].en = !status;
            status = !status; // ??
        } 
        if (verts[current].isect) {
            printf("setting en status of vert %d\n", current);
            
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
                printf("visiting vertex %d\n", si);
                verts[si].next_poly = -2; // make sure
                done = false;
                
                // in the GH paper, the current vertex is added here
                // but this produces duplication (first==last), so we just skip it here ...
                int current = si;
                do {
                    verts[current].next_poly = -2;
                    if (verts[current].en) {
                        do {
                            current = verts[current].next;
                            printf("en: visiting %d\n", current);
                            vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                            verts[current].next_poly = -2;
                        } while (!verts[current].isect);
                    } else {
                        do {
                            current = verts[current].prev;
                            printf("ex: visiting %d\n", current);
                            vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                            verts[current].next_poly = -2;
                        } while(!verts[current].isect);
                    }
                    current = verts[current].neighbour; // swap to the other poly
                    printf("switched to neighbour %d\n", current);
                } while (current != si); 
                if (vertices.size() > 0) {
                    vector<cv::Vec2d> nverts;
                    int prev = vertices.size() - 1;
                    for (size_t vv=0; vv < vertices.size(); vv++) {
                        if (norm(vertices[prev], vertices[vv]) > 1e-11) {
                            nverts.push_back(vertices[vv]);
                        }
                        prev = vv;
                    }
                    polys.push_back(Polygon_geom(nverts));
                    vertices.clear();
                }
            }
        }
        
    }
}

}; // namespace GH_clipping

