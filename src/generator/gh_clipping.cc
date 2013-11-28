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
    
    
    if (fabs(is.alpha) < 1e-11 && fabs(ic.alpha - 1) < 1e-11) { // classical vertex intersection
        return false;
    }
    
    if (fabs(ic.alpha) < 1e-11 && fabs(is.alpha - 1) < 1e-11) { // same, other way around
        return false;
    }
    
    // should we check that the alphas are in the range [0,1] ?
    if (is.alpha < 0 || is.alpha > 1 ||  
        ic.alpha < 0 || ic.alpha > 1) {
    
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
    verts[vs].flag = NONE;
    verts[vs].couple = -1;
    verts[vs].cross_change = false;
}

static inline double dist(const gh_vertex& a, const gh_vertex& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

int check_degeneracy(vector<gh_vertex>& verts, int vs, int si, int ci, int nsi, int nci) {
    int& cur = vs;
    int& si_prev = verts[si].prev;
    int& si_next = nsi; 
    
    int& ci_prev = verts[ci].prev;
    int& ci_next = nci;
    
    if (dist(verts[vs], verts[si_prev]) < 1e-11) {
        printf("new vert %d identical to prev S vertex %d\n", vs, si_prev);
        return si_prev;
    }
    
    if (dist(verts[vs], verts[ci_prev]) < 1e-11) {
        printf("new vert %d identical to prev C vertex %d\n", vs, ci_prev);
        return ci_prev;
    }
    
    if (dist(verts[vs], verts[si]) < 1e-11) {
        printf("new vert %d identical to current S vertex %d\n", vs, si);
        return si;
    }
    
    if (dist(verts[vs], verts[ci]) < 1e-11) {
        printf("new vert %d identical to current C vertex %d\n", vs, ci);
        return ci;
    }
    
    
    if (dist(verts[vs], verts[si_next]) < 1e-11) {
        printf("new vert %d identical to next S vertex %d\n", vs, si_next);
        return si_next;
    }
    
    if (dist(verts[vs], verts[ci_next]) < 1e-11) {
        printf("new vert %d identical to next C vertex %d\n", vs, ci_next);
        return ci_next;
    }
    
    
    for (int k=0; k < vs; k++) {
        if (dist(verts[vs], verts[k]) < 1e-11) {
            printf("new vert %d identical to arb C vertex %d\n", vs, k);
            return k;
        }
    }
    
    return -1;
}


int init_gh_list(vector<gh_vertex>& verts, const vector<cv::Vec2d>& in_verts, int vs, int next_poly) {
    for (int v=0; v < (int)in_verts.size(); v++) {
        verts[v+vs].x = in_verts[v][0];
        verts[v+vs].y = in_verts[v][1];
        verts[v+vs].next = (v+1) % in_verts.size() + vs;
        verts[v+vs].prev = (v + in_verts.size() - 1) % in_verts.size() + vs;
        verts[v+vs].isect = false;  // this it not an intersection, so it is a vertex
        verts[v+vs].flag = NONE;     // thus, en is the inside/outside flag???
        verts[v+vs].neighbour = -1; // no matching vertex in neighbour yet
        verts[v+vs].next_poly = next_poly; // no next poly yet
        verts[v+vs].alpha = 10;     // any value > 1 indicates no intersection, which is the default
        verts[v+vs].couple = -1;
        verts[v+vs].cross_change = false;
        
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
            
            printf("intersect (%d) %d->%d and %d->%d, alpha=%lf, %lf\n", isect, si, nsi, ci+poly0_size, nci+poly0_size, verts[vs].alpha, verts[vs+1].alpha);
            
            if (isect &&
                vs - 2 > (poly1_size + poly0_size) &&
                dist(verts[vs], verts[vs-2]) < 1e-11) {
                
                printf("same intersection as the one we have just added. skipping\n");
                isect = false;
            }
            
            if (isect) { 
                verts[vs].neighbour = vs+1;
                verts[vs].next_poly = 1;
                verts[vs+1].neighbour = vs;
                verts[vs+1].next_poly = -7;
                
                int deg = check_degeneracy(verts, vs, si, ci+poly0_size, nsi, nci+poly0_size);
                
                if (deg >= 0) { // intersection already exists
                    if (verts[deg].neighbour == -1) { // if vertex already has a neighbour, just skip
                        if (deg < poly0_size) { // match comes from S, so add intersection to C
                            // check if it already exists in C
                            bool found = false;
                            
                            
                            int found_idx = poly0_size;
                            for (; found_idx < (poly0_size+poly1_size); found_idx++) {
                                if (dist(verts[found_idx], verts[deg]) < 1e-11) {
                                    found = true;
                                    break;
                                }
                            }
                            
                            int dest = vs+1;
                            
                            if (found) {
                                printf("found same vertex in other poly S : %d!\n", found_idx);
                                dest = found_idx;
                            } else {
                                insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                            }
                            
                            verts[dest].neighbour = deg;
                            verts[dest].isect = 1;
                            
                            verts[deg].neighbour = dest;
                            verts[deg].isect = 1;
                            
                            verts[vs].isect = 0; // disable redundant vertex
                            verts[vs].next_poly = -20;
                            
                            if (!found) vs += 2;
                        } else { // match comes from C, so add intersection to S
                            // check if it already exists in S
                            
                            bool found = false;
                            int found_idx = 0;
                            for (; found_idx < poly0_size; found_idx++) {
                                if (dist(verts[found_idx], verts[deg]) < 1e-11) {
                                    found = true;
                                    break;
                                }
                            }
                            
                            int dest = vs;
                            
                            if (found) {
                                printf("found same vertex in other poly C : %d !\n", found_idx);
                                dest = found_idx;
                            } else {
                                insert_vertex(verts, si, nsi, vs);
                            }
                            
                            verts[deg].neighbour = dest;
                            verts[deg].isect = 1;
                            
                            verts[dest].neighbour = deg;
                            verts[dest].isect = 1;
                              
                            verts[vs+1].isect = 0; // disable redundant vertex
                            verts[vs+1].next_poly = -20;
                            
                            if (!found) vs += 2;
                            
                        }
                    } else {
                        printf("skipping ...\n");
                    }
                } else {
                    insert_vertex(verts, si, nsi, vs);
                    insert_vertex(verts, ci+poly0_size, nci+poly0_size, vs+1);
                    vs += 2;
                }
                
            }
        }
    }
}

inline static double triangle_cross(const gh_vertex& v0, const gh_vertex& v1, const gh_vertex& v2) {
    double ds_x = v1.x - v0.x;
    double ds_y = v1.y - v0.y;
    double dc_x = v2.x - v0.x;
    double dc_y = v2.y - v0.y;
                                      
    return (dc_y*ds_x - dc_x*ds_y);
}

inline static edge_status classify_edge(int stat) {
    edge_status ret = OUT;
    switch(stat) {
    case 8: // on/on
        ret = ON;
        break;
    case 4: // in/in
        ret = IN;
        break;
    case 0: // out/out
        ret = OUT;
        break;
        
    case 5: // in/on
        ret = IN; // TODO: ON?
        break;
    case 2: // out/on
        ret = OUT; // TODO: ON?
        break;
        
    case 6: // on/out
        ret = OUT; // TODO: ON?
        break;
    case 7: // on/in
        ret = IN; // TODO: ON?
        break;
        
    case 1: // out/in
    case 3: // in/out
        printf("This should not happen!\n");
        ret = ON;
        break;    
    }
    
    return ret;
}

bool edge_present(vector<gh_vertex>& verts, int other, int current) {

    int neighbour = verts[current].neighbour;

    printf("checking if %d -> %d falls on other poly (neighbour=%d)\n",
        other, current, neighbour);


    if (neighbour >= 0) { // we have a neighbour
        // take the expensive O(n) route to check if other->current is an edge in the other polygon
        int nc = verts[current].neighbour;
        do {
            if (dist(verts[nc], verts[current]) < 1e-11) {
                
                // now check if prev/next of found_idx matches "other"
                if (dist(verts[verts[nc].prev], verts[other]) < 1e-11 ||
                    dist(verts[verts[nc].next], verts[other]) < 1e-11) {
                    
                    return true;
                }
                
            }
            nc = verts[nc].next;
        } while (nc != verts[current].neighbour);
    } 
    
    return false;
}

static void set_traversal_flag(vector<gh_vertex>& verts, const Polygon_geom& b, int current) {
    
    const int& prev = verts[current].prev;
    const int& next = verts[current].next;
    
    const int& neigh = verts[current].neighbour;
    const int& nprev = verts[neigh].prev;
    const int& nnext = verts[neigh].next;
    
    
    // classify the prev->current edge
    edge_status prev_status = edge_present(verts, prev, current) ? ON : OUT;
    printf("prev_status (%d) : %d\n", current, prev_status);
    if (prev_status != ON) {
        switch (
            b.classify( 
                0.5*(verts[prev].x + verts[current].x), 
                0.5*(verts[prev].y + verts[current].y)
            ) 
        ) {
        case ON:
            printf("Warning! (prev) This should not be possible\n");
            break;
        case INSIDE:
            prev_status = IN;
            break;
        case OUTSIDE:
            prev_status = OUT;
        }
    }
    
    // classify the current->next edge
    edge_status next_status = edge_present(verts, next, current) ? ON : OUT;
    printf("next_status (%d) : %d\n", current, next_status);
    if (next_status != ON) {
        switch (
            b.classify( 
                0.5*(verts[next].x + verts[current].x), 
                0.5*(verts[next].y + verts[current].y)
            ) 
        ) {
        case ON:
            printf("Warning! (next) This should not be possible\n");
            break;
        case INSIDE:
            next_status = IN;
            break;
        case OUTSIDE:
            next_status = OUT;
        }
    }
    
    traversal_flag& trav = verts[current].flag;
    
    int combined_edge_status = prev_status * 3 + next_status;
    printf("combined_edge_status=%d\n", combined_edge_status);
    switch (combined_edge_status) {
    case 8: // on/on
        trav = NONE;
        break;
        
    case 6: // on/out
    case 5: // in/on
    case 3: // in/out
        trav = EX;
        break;
    
    case 7: // on/in
    case 2: // out/on
    case 1: // out/in
        trav = EN;
        break;
        
    case 4: // in/in
        trav = EXEN;
        break;
    
    case 0: // out/out
        trav = ENEX;
        break;
    }
    
    // compute cross_change flag. probably not needed for every vertex ...
    if (trav == EXEN || trav == ENEX) {
        // TODO: assume this function is called with C (poly 1), so S is neigh, C is self
        double tr1_cross = triangle_cross(verts[nprev], verts[prev], verts[nnext]);
        double tr2_cross = triangle_cross(verts[nprev], verts[prev], verts[next]);
        verts[current].cross_change = tr1_cross * tr2_cross < 0;
    }
}

void gh_phase_two(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index) {

    // set the traversal flags for C (poly 1)
    int current = first_vert_index;
    do {
        if (verts[current].isect) {
            printf("setting en flag of vert %d (before=%d)\n", current, verts[current].flag);
            set_traversal_flag(verts, b, current);
            printf("\tafter: flag=%d\n", verts[current].flag);
        }
        
        current = verts[current].next;
    } while (current != first_vert_index); 
    
    // now set the "couple" field
    current = first_vert_index;
    do {
        if ( verts[current].isect && (verts[current].flag == EN || verts[current].flag == EX) ) {
            int next = verts[current].next;
            bool done = false;
            
            printf("couple: current=%d\n", current);
            
            //while (next != current && !done) { // TODO: what happens when we have three colinear points on the edge of the other poly?
                if (verts[next].isect && verts[next].flag == verts[current].flag && 
                    verts[current].couple == -1 && verts[next].couple == -1) { // only couple uncoupled vertices?
                    done = true;
                    printf("coupling vertices %d and %d\n", next, current);
                    verts[current].couple = next;
                    verts[next].couple = current;
                    
                    // should we couple the same vertices in the other polygon?
                    verts[verts[current].neighbour].couple = verts[next].neighbour;
                    verts[verts[next].neighbour].couple = verts[current].neighbour;
                    printf(" ** coupling vertices %d and %d\n", verts[next].neighbour, verts[current].neighbour);
                }
              //  next = verts[next].next;
            //}
        }
        
        current = verts[current].next;
    } while (current != first_vert_index); 
}

void gh_phase_two_b(vector<gh_vertex>& verts, const Polygon_geom& b, int first_vert_index) {

    printf("### about to set traversal flags in Poly 1\n");

    // set the traversal flags for C (poly 1)
    int current = first_vert_index;
    int start = -1;
    do {
        if (verts[current].isect) {
            start = current;
            break;
        }
        current = verts[current].next;
    } while (current != first_vert_index); 
    
    bool reverse = false;
    if (start >= 0) {
        set_traversal_flag(verts, b, start);
    
        if ( (verts[start].flag == EN && verts[verts[start].neighbour].flag == EX) ||
             (verts[start].flag == EX && verts[verts[start].neighbour].flag == EN) || 
             (verts[start].flag == ENEX && verts[verts[start].neighbour].flag == EXEN) ) { // TODO: does EXEN/ENEX ever occur?
            reverse = true;
        } else {
            if (verts[start].flag == verts[verts[start].neighbour].flag) {
                reverse = false;
            } else {
                // everything else, assume revers is true, but we would still
                // like to know about it ...
                printf("flags in undefined state (%d, %d)->%d. help!\n",
                    verts[start].flag, verts[verts[start].neighbour].flag, true);
                reverse = true;
            }
        }
        
        current = start;
        do {
            if (verts[current].isect) {
                printf("### visiting %d, flag=%d, neigh=%d, nflag=%d\n",
                    current, verts[current].flag, verts[current].neighbour, verts[verts[current].neighbour].flag);
                if (reverse) {
                    traversal_flag nflag = NONE;
                    switch (verts[verts[current].neighbour].flag) {
                    case NONE:
                        nflag = NONE; // ?? // TODO: this could be a problem ...
                        break;
                    case EN:
                        nflag = EX;
                        break;
                    case EX:
                        nflag = EN;
                        break;
                    case EXEN:
                        nflag = ENEX;
                        break;
                    case ENEX:
                        nflag = EXEN;
                        break;
                    }
                    verts[current].flag = nflag;
                } else {
                    verts[current].flag = verts[verts[current].neighbour].flag;
                }
                printf("### AFTER: visiting %d, flag=%d, neigh=%d, nflag=%d\n",
                    current, verts[current].flag, verts[current].neighbour, verts[verts[current].neighbour].flag);
            }
            current = verts[current].next;
        } while (current != start);
        
    } else {
        printf("no intersections in subject poly. weird.\n");
    }
}

static int proceed(int current, traversal_edge status, vector<gh_vertex>& verts, vector<cv::Vec2d>& vertices) {
    printf("when visiting %d, status == %d\n", current, status);
    if (status == D1) {
        vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
        printf("visiting %d, direction next\n", current);
        return verts[current].next;
    }
    if (status == D2) {
        vertices.push_back(cv::Vec2d(verts[current].x, verts[current].y));
        printf("visiting %d, direction prev\n", current);
        return verts[current].prev;
    }
    // touch both intersections
    verts[current].next_poly = -2;
    verts[verts[current].neighbour].next_poly = -2; 
    
    printf("switching to neighbour %d\n", verts[current].neighbour);
    return verts[current].neighbour;
}

static traversal_edge delete_flag1(int current, traversal_edge status, vector<gh_vertex>& verts) {
    printf("delete_flag1 on vertex %d, status=%d, flag=%d\n", current, status, verts[current].flag);
    switch (verts[current].flag) {
    case ENEX:
        verts[current].flag = NONE;
        if (verts[current].cross_change) {
            return (status == D3) ? D3 : D4;
        }
        return (status == D3) ? D4 : D3;
        break;
    case EXEN:
        verts[current].flag = (status == D3) ? EN : EX;
        return (status == D3) ? D2 : D1;
        break;
    case EN:
        verts[current].flag = NONE;
        return D1;
        break;
    case EX:
        verts[current].flag = NONE;
        return D2;
        break;
    case NONE: // never reached, parent skips this case
        break;
    }
    //assert(false);
    return D1; // not usually reachable
}

static traversal_edge delete_flag2(int current, int prev, traversal_edge status, vector<gh_vertex>& verts) {
    printf("delete_flag2 on vertex %d, status=%d, flag=%d\n", current, status, verts[current].flag);
    switch (verts[current].flag) {
    case ENEX:
        verts[current].flag = (status == D1) ? EX : EN;
        if (verts[current].cross_change) {
            return (status == D1) ? D4 : D3;
        }
        return (status == D1) ? D3 : D4;
        break;
    case EXEN:
        verts[current].flag = (status == D1) ? EN : EX;
        if (verts[current].cross_change) {
            return (status == D1) ? D4 : D3;
        }
        return (status == D1) ? D3 : D4;
        break;
    case EN:
        verts[current].flag = NONE;
        if (status == D1 && 
            verts[current].couple != -1 && 
            verts[prev].couple == current) {
            
            // mark both vertices in couple to prevent them from being selected
            verts[current].next_poly = -4;
            verts[verts[current].couple].next_poly = -5; //?
            
            return D1; // keep going in the same direction when encountering a coupled vertex
        }
        return (status == D1) ? D3 : D4;
        break;
    case EX:
        verts[current].flag = NONE;
        if (status != D1 && 
            verts[current].couple != -1 &&
            verts[prev].couple == current) {
            
            // mark both vertices in couple to prevent them from being selected
            verts[current].next_poly = -4;
            verts[verts[current].couple].next_poly = -5; //?
            
            return D2;
        }
        return (status == D1) ? D3 : D4;
        break;
    case NONE: // never reached, parent skips this case
        break;
    }
    //assert(false);
    return D1; // not usually reachable
}

int select_vertex(vector<gh_vertex>& verts, int vs) {
    
    for (int si=0; si < vs; si++) {
        if (verts[si].next_poly == 1 && verts[si].isect) {
        
            // we have an unmarked intersection
            
            if (verts[si].couple != -1) {  // si is part of an unmarked couple
                if (verts[si].flag == EN && verts[verts[si].couple].flag == EN) {
                    printf("**selecting vertex (2nd in couple) %d as starting point\n", verts[si].couple);
                    return verts[si].couple;
                } else {
                    if (verts[si].flag == EX && verts[verts[si].couple].flag == EX)  {
                        printf("**selecting vertex (1st in couple) %d as starting point\n", si);
                        return si;
                    } else {
                        // we reach this state if one of the couple was accessed through the neighbour pointer
                        verts[si].next_poly = -11;
                        continue; // try next si?
                        //assert(false);
                    }
                }
            } 
            
            printf("**selecting vertex %d as starting point\n", si);
            return si;
        }
    }
    
    return -1;
}

void gh_phase_three(vector<gh_vertex>& verts, int vs, int first_isect_index, vector<Polygon_geom>& polys) {

    bool done = false;
    while (!done) {
        vector<cv::Vec2d> vertices;
        
        int current = select_vertex(verts, vs);
        done = current < 0;
        if (!done) {
            verts[current].next_poly = -10;
            
            int prev = current;
            int start = current;
            
            traversal_edge status = delete_flag1(current, D3, verts);
            current = proceed(current, status, verts, vertices);
            
            while (current != start) {
                if (verts[current].flag) {
                    if (status == D1 || status == D2) {
                        status = delete_flag2(current, prev, status, verts); // Kim does not appear to set status here ?
                    } else {
                        status = delete_flag1(current, status, verts); // or here?
                    }
                    prev = current;
                }
                current = proceed(current, status, verts, vertices);
            }
            printf("got %d vertices for next poly\n", (int)vertices.size());
            if (vertices.size() > 2) {
                // we should scan for colinearity here, just to keep the polygons simpler
                polys.push_back(Polygon_geom(vertices));
            }
            vertices.clear();
        }
        
    }
    
}

}; // namespace GH_clipping

