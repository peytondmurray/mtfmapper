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

static inline bool on_edge(const gh_vertex& s0, const gh_vertex& s1,
                        const gh_vertex& p0) {
                  
    double ds_x = s1.x - s0.x;
    double ds_y = s1.y - s0.y;
    double dc_x = p0.x - s0.x;
    double dc_y = p0.y - s0.y;
                                      
    double denom = (dc_y*ds_x - dc_x*ds_y);
                                                      
    if (fabs(denom) < 1e-11) {
        // vectors are parallel, so check if p0 falls between s0 and s1
        double dot = (ds_x*dc_x + ds_y*dc_y)/sqrt(ds_x*ds_x + ds_y*ds_y);
        if (dot < 0) return false;
        return dot <= 1;
    }
    
    return false;        
}


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
    if (is.alpha < 0 || is.alpha >= 1 ||  
        ic.alpha < 0 || ic.alpha >= 1) {
    
        //printf("intersection outside of end vertices, not reporting this as an intersection\n");
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

static double dist(gh_vertex& a, gh_vertex& b) {
    return sqrt( (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)  );
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

static int check_degeneracy(vector<gh_vertex>& verts, int v1, int v2, bool skip_next=false) {
    /*
    printf("checking (%lf, %lf) -> (%lf, %lf) / (%lf, %lf) -> (%lf, %lf)\n",
        verts[v1].x, verts[v1].y,
        verts[verts[v1].next].x, verts[verts[v1].next].y,
        verts[v2].x, verts[v2].y,
        verts[verts[v2].next].x, verts[verts[v2].next].y
    );
    */
    bool isect = true;
    if (dist(verts[v1], verts[v2]) < 1e-11) {
        printf("##1 identical vertex found in first poly: %lf %lf (%d, %d)\n",verts[v1].x, verts[v1].y, v1, v2);
        verts[v1].neighbour = v2;
        verts[v2].neighbour = v1;
        verts[v1].isect = true;
        verts[v2].isect = true;
        verts[v1].en = true;
        verts[v2].en = true;
        isect = false; // TODO: should this return? orig code falls through (check others below too)
    }
    if (!skip_next && dist(verts[verts[v1].next], verts[v2]) < 1e-11) {
        int& next = verts[v1].next;
        //printf("##2 identical vertex found in first poly: %lf %lf\n",verts[v1].x, verts[v1].y);
        verts[next].neighbour = v2;
        verts[v2].neighbour = next;
        verts[next].isect = true;
        verts[v2].isect = true;
        verts[next].en = true;
        verts[v2].en = true;
        isect = false;
    }
    
    if (!skip_next && dist(verts[v1], verts[verts[v2].next]) < 1e-11) {
        int& next = verts[v2].next;
        //printf("##3 identical vertex found in first poly: %lf %lf\n",verts[next].x, verts[next].y);
        verts[v1].neighbour = next;
        verts[next].neighbour =v1;
        verts[v1].isect = true;
        verts[next].isect = true;
        verts[v1].en = true;
        verts[next].en = true;
        isect = false;
    }
    if (!skip_next && dist(verts[verts[v1].next], verts[verts[v2].next]) < 1e-11) {
        int& next = verts[v1].next;
        int& nnext = verts[v2].next;
        //printf("##4 identical vertex found in first poly: %lf %lf\n",verts[v1].x, verts[v1].y);
        verts[next].neighbour = nnext;
        verts[nnext].neighbour = next;
        verts[next].isect = true;
        verts[nnext].isect = true;
        verts[next].en = true;
        verts[nnext].en = true;
        isect = false;
    }
    
    if (isect) { // not one of the vertex matching cases above
        if (on_edge(verts[v1], verts[verts[v1].next], verts[v2]) ) {
            printf("##$vertex %d lies on edge from %d to %d\n", v2, v1, verts[v1].next);
        }
    }
    
    return isect;
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
        
        if (next_poly != 1) { // while processing the second polygon
            for (int v0=0; v0 < vs; v0++) {
               
                bool isect = check_degeneracy(verts, v0, v+vs, true);
                if (!isect) {
                    printf("matched second poly vertex with %d\n", v0);
                    // since toggling the isect flag seems to break one of the test cases
                    // the fault must be in phase two ... ?
                }
                
            }
        }
    }
    return vs + in_verts.size(); // index of next free vertex entry
}



void gh_phase_one(vector<gh_vertex>& verts, int& vs, int poly0_size, int poly1_size) {
    // phase one
    for (int si=0; si < poly0_size; si++) {
        int nsi = (si + 1) % poly0_size;
        for (int ci=0; ci < poly1_size; ci++) {
            int nci = (ci + 1) % poly1_size;
            
            //printf("testing %d from p0, %d from p1\n", si, ci);
            
            bool isect = t_intersect(
                verts[si], verts[nsi], 
                verts[ci+poly0_size], verts[nci+poly0_size], 
                verts[vs], verts[vs+1]
            );
            
            if (isect && 
                //(verts[vs].alpha < 1e-11 || verts[vs+1].alpha < 1e-11 ||
                // (1-verts[vs].alpha) < 1e-11 || (1-verts[vs+1].alpha) < 1e-11) ) { // TODO: could modify intersection and treat alpha=1 here
                (dist(verts[si], verts[vs]) < 1e-11 || dist(verts[ci+poly0_size], verts[vs]) < 1e-11 ||
                 dist(verts[verts[si].next], verts[vs]) < 1e-11 || dist(verts[verts[ci+poly0_size].next], verts[vs]) < 1e-11
                 ) ) {
            
                printf("## new vertex identical to subject start (alpha1=%lf, alpha2=%lf) : %lf %lf\n", 
                    verts[vs].alpha, verts[vs+1].alpha, verts[vs].x, verts[vs].y);
                
                // if we simply skip adding this intersection, we appear
                // to be fine (except for complete overlap)
                
                isect = check_degeneracy(verts, si, ci+poly0_size); 
                if (!isect) {
                    //printf("we found vertices in the S and C polygon that are identical, so we fudged the neighbour pointers accordingly\n");
                } else {      
                    //printf("intersection is close to one of the starting points, but the starting points are not the same...\n");
                    
                    // we should insert a single vertex into the relevant polygon
                    // and link it up with the matching neighbour point
                    
                    if (dist(verts[si], verts[vs]) < 1e-11 || dist(verts[verts[si].next], verts[vs]) < 1e-11) {
                        // the matching vertex was in the subject polygon
                        //printf("inserting vertex %d between %d and %d in poly1\n", vs+1, ci+poly0_size, nci+poly0_size);
                        insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                        
                        
                        verts[vs+1].neighbour = (dist(verts[si], verts[vs]) < 1e-11) ? si : verts[si].next;
                        verts[vs+1].next_poly = -1;
                        verts[verts[vs+1].neighbour].neighbour = vs+1;
                        verts[verts[vs+1].neighbour].isect = 1;
                        
                        verts[vs].next_poly = -5;
                        verts[vs].isect = 0;
                        
                        vs += 2; // this wastes slot vs, but it is not linked to anything, so no harm done
                        
                    } else {
                        // the matching vertex was in the clip polygon
                        //printf("inserting vertex %d between %d and %d in poly0\n", vs, si, nsi);
                        insert_vertex(verts, si, nsi, vs);
                        
                        verts[vs].neighbour = (dist(verts[ci+poly0_size], verts[vs]) < 1e-11) ? ci+poly0_size : verts[ci+poly0_size].next;
                        verts[vs].next_poly = 1;
                        verts[verts[vs].neighbour].neighbour = vs;
                        verts[verts[vs].neighbour].isect = 1;
                        
                        verts[vs+1].next_poly = -5;
                        verts[vs+1].isect = 0;
                        
                        vs += 2; // this wastes slot vs, but it is not linked to anything, so no harm done
                    }
                    
                    isect = false; // 
                }
                
                // 140, 90 lies on edge, but is not an intersection. why not?
                
            }
            
            
            //printf("result: %d (%lf, %lf), (%lf, %lf)\n", isect, verts[vs].x, verts[vs].y, verts[vs+1].x, verts[vs+1].y);
                
            if (isect) { 
                verts[vs].neighbour = vs+1;
                verts[vs].next_poly = 1;
                verts[vs+1].neighbour = vs;
                verts[vs+1].next_poly = -7;
                
                insert_vertex(verts, si, nsi, vs);
                insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                
                vs += 2;
            }
        }
    }
}

void gh_phase_two_a(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index, int nverts) {
    // If a vertex is "ON", we force it to inside?
    for (int i=first_vert_index; i < (first_vert_index+nverts); i++) {
        int status = b.classify(verts[i].x, verts[i].y);
        //verts[i].en = status == ON ? 1 : status;
        if (status == 2) {
            printf("## vertex is ON\n");
        }
        verts[i].en = status;
    }
    
}

static void remove_intersection(vector<gh_vertex>& verts, int current) {
    
    printf("removing intersection %d\n", current);
    
    // we do not actually remove the vertex
    // but merely turn it into non-intersection vertex
    
    int neigh = verts[current].neighbour;
    
    verts[current].en = verts[verts[current].prev].en;  // copy predecessor's en flag
    verts[neigh].en   = verts[verts[neigh].prev].en;    // same for neighbour
    
    verts[current].neighbour = -1;
    verts[neigh].neighbour = -1;
    
    //verts[verts[current].prev].next = verts[current].next; 
    //verts[verts[current].next].prev = verts[current].prev;
    
    // verts[current].next_poly = -3; // not really needed since we nuke the isect flag
    
    verts[current].isect = false;
    verts[neigh].isect = false;
}

void update_en_flag(vector<gh_vertex>& verts, int current) {
    printf("intersection (%lf, %lf): en=%d (0=exit, 1=entry)\n", verts[current].x, verts[current].y, verts[current].en);
    
    int& next = verts[current].next;
    int& prev = verts[current].prev;
    int& neigh = verts[current].neighbour;
    int& nprev = verts[neigh].prev;
    int& nnext = verts[neigh].next;
    
    printf("[%d]: prev(%d): isect=%d, en=%d, next(%d): isect=%d, en=%d\n",
        current,  
        prev,
        verts[prev].isect, verts[prev].en,
        next,
        verts[next].isect, verts[next].en
    );
    
    // Intersections are counted as "ON", but individual vertices could also be "ON" ??
    int sbits = 0;
    sbits += verts[prev].isect ? ON : verts[prev].en;
    sbits += (verts[next].isect ? ON : verts[next].en) * 3;
    
    int nsbits = 0;
    nsbits += verts[nprev].isect ? ON : verts[nprev].en;
    nsbits += (verts[nnext].isect ? ON : verts[nnext].en) * 3;
    
    printf("\tsbits=%d, nsbits=%d\n", sbits, nsbits);
    
    bool polarity = false;
    switch(sbits) {
    case 0: // out/out
        polarity = true;
    case 4: // in/in
        
        switch(nsbits){
        case 8: // on/on
            remove_intersection(verts, current);
            verts[neigh].en = polarity;
            break;
        case 0: // out/out
        case 4: // in/in
            remove_intersection(verts, current);
            break;
        default:
            // one of these cases should lead to removal of the intersection
            // how do we identify this case? // first check if vertex==ON fixes this
            if (verts[nprev].en && !verts[nnext].en) {
                verts[current].en = true;
            } else {
                verts[current].en = false;
            }
            // TODO: fix this mess
            /*
            printf("here: %d %d\n", verts[current].next_poly, verts[neigh].next_poly);
            if (verts[nprev].en && !verts[nnext].en) {
                verts[current].en = (verts[current].next_poly == 1 &&  verts[neigh].next_poly == -1) ? false : true;
            } else {
                printf("here too\n");
                verts[current].en = (verts[current].next_poly == 1 &&  verts[neigh].next_poly == -1) ? true : false;
            }
            */
        };
        break;
        
    case 8: // on/on  
    
        switch(nsbits) {
        case 0: // out/out
            remove_intersection(verts, current);
            verts[current].en = true;
            verts[neigh].en = false;
            break;
        case 4: // in/in
            remove_intersection(verts, current);
            verts[current].en = false;
            verts[neigh].en = true;
            break;
        case 8: // on/on
            remove_intersection(verts, current);
            verts[current].en = true;
            verts[neigh].en = true;
            break;
        case 2: // on/out
            verts[current].en = false;
            break;
        case 5: // on/in
            verts[current].en = true;
            break;
        case 6: // out/on
            printf("here!\n");
            verts[current].en = false;
            break;
        case 3: // out/in
            verts[current].en = false;
            break;
        case 7: // in/on
            verts[current].en = true;
            break;
        case 1: // in/out
            verts[current].en = true;
            break;
        };
        
        break;
    case 1: // in/out
    case 2: // on/out
    case 7: // in/on
        verts[current].en = false;
        break;
    
    case 5: // on/in
    case 3: // out/in
    case 6: // out/on
        printf("here2: nsbits=%d\n", nsbits);
        // this exact case seems to work now ...
        if (nsbits == 8) { // the neighbour is between two "ON" vertices, and is itself an intersection
            verts[neigh].en = true; // force neighbour to be "ON", which will subsequently lead to removal of this intersection?
        }
        verts[current].en = true;
        break;
    }
    
    printf("(%d): sbits=%d, nsbits=%d, outflag=%d\n", current, sbits, nsbits, verts[current].en);
}

void gh_phase_two(vector<gh_vertex>& verts, int first_vert_index = 0) {
    //printf("first vertex inside other poly: %d\n", verts[first_vert_index].en);
    
    int current = first_vert_index;
    do {
        //printf("current vert is %d, (n=%d, p=%d), isect=%d\n", current, verts[current].next, verts[current].prev, verts[current].isect);
        if (verts[current].isect) {
            
            update_en_flag(verts, current);
            
            // somehow this case fell through. The en flags are not the problem
            // which means that remove_intersection() was not called, or failed somehow ...
            //verts[3].en = 1; // fudge it, see if that helps
            //verts[4].en = 1;
            
            if (verts[current].isect) { // if we did not remove current
                update_en_flag(verts, verts[current].neighbour);
                printf("*after: (%d) (%lf, %lf): en=%d isect=%d\n", current, verts[current].x, verts[current].y, verts[current].en, verts[current].isect);
                printf("**after: (%d) (%lf, %lf): en=%d isect=%d\n", current, verts[verts[current].neighbour].x, verts[verts[current].neighbour].y, verts[verts[current].neighbour].en, verts[verts[current].neighbour].isect);
                if (verts[current].en == verts[verts[current].neighbour].en) { // this might be the problem
                    printf("calling remove on vertex %d\n", current);
                    remove_intersection(verts, current);
                }
            }
            
            //printf("*after: intersection(%d) (%lf, %lf): en=%d (0=exit, 1=entry)\n", current, verts[current].x, verts[current].y, verts[current].en);
        }
        
        current = verts[current].next;
    } while (current != first_vert_index); // let us hope the chain is not broken ...
    
    
}


void gh_phase_two_old(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index) {
    bool status = b.classify(verts[first_vert_index].x, verts[first_vert_index].y) == INSIDE ? false : true;
    
    printf("first vertex status (i.e., first isect en flag): %d\n", status);
    
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
                            printf("en: visiting vertex %d, en=%d, isect=%d\n", current, verts[current].en, verts[current].isect);
                            vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                            verts[current].next_poly = -2;
                        } while (!verts[current].isect);
                    } else {
                        do {
                            current = verts[current].prev;
                            printf("ex: visiting vertex %d, en=%d, isect=%d\n", current, verts[current].en, verts[current].isect);
                            vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
                            verts[current].next_poly = -2;
                        } while(!verts[current].isect);
                    }
                    printf("swapping from %d to neighbour %d (neighbour.en=%d)\n", current, verts[current].neighbour, verts[verts[current].neighbour].en);
                    current = verts[current].neighbour; // swap to the other poly
                } while (current != si); // TODO: will this work if the last intersection was in clip poly?
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

