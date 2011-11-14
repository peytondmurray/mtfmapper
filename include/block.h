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
#ifndef BLOCK_H
#define BLOCK_H


#include "include/rectangle.h"

#include <map>
using std::map;

class Block {
  public:
	typedef enum {TOP, LEFT, RIGHT, BOTTOM} edge_position;    

    Block(void) : rect(Mrectangle()), mtf50(4,0.0), quality(4, 0.0), centroid(0,0) {
    }

    Block(const Mrectangle& in_rect) : rect(in_rect), mtf50(4,0.0), quality(4, 0.0), centroid(0,0) {
    
        size_t top_edge_idx = 0;
        size_t bot_edge_idx = 0;
        size_t left_edge_idx = 0;
        size_t right_edge_idx = 0;
     
        for (size_t i=0; i < 4; i++) {
            centroid.x += get_edge_centroid(i).x;
            centroid.y += get_edge_centroid(i).y;
        
            if (get_edge_centroid(i).y < get_edge_centroid(top_edge_idx).y) {
                top_edge_idx = i;
            }
            
            if (get_edge_centroid(i).y > get_edge_centroid(bot_edge_idx).y) {
                bot_edge_idx = i;
            }
            
            if (get_edge_centroid(i).x < get_edge_centroid(left_edge_idx).x) {
                left_edge_idx = i;
            }
            
            if (get_edge_centroid(i).x > get_edge_centroid(right_edge_idx).x) {
                right_edge_idx = i;
            }
        }
        centroid.x /= 4;
        centroid.y /= 4;
        
        edge_lut[TOP]    = top_edge_idx;
        edge_lut[BOTTOM] = bot_edge_idx;
        edge_lut[LEFT]   = left_edge_idx;
        edge_lut[RIGHT]  = right_edge_idx;
        
        Point e1(
            get_edge_centroid(edge_lut[BOTTOM]).x - get_edge_centroid(edge_lut[TOP]).x,
            get_edge_centroid(edge_lut[BOTTOM]).y - get_edge_centroid(edge_lut[TOP]).y
        );
        
        Point e2(
            get_edge_centroid(edge_lut[RIGHT]).x - get_edge_centroid(edge_lut[LEFT]).x,
            get_edge_centroid(edge_lut[RIGHT]).y - get_edge_centroid(edge_lut[LEFT]).y
        );
        
        area = sqrt(SQR(e1.x) + SQR(e1.y)) * sqrt(SQR(e2.x) + SQR(e2.y));
    }

    void set_normal(size_t edge_number, const Point& rgrad) {
        rect.normals[edge_number] = rgrad;
    }

    Point get_normal(size_t edge_number) const {
        return rect.normals[edge_number];
    }

    Point get_edge(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.edges[edge_number];
    }
    
    double get_edge_angle(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.thetas[edge_number];
    }
    
    Point get_edge_centroid(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.centroids[edge_number];
    }
    
    Point get_edge_centroid(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return rect.centroids[it->second];
    }
    
    // quality of 1.0 means good quality, 0.0 means unusably poor
    void set_mtf50_value(size_t edge_number, double mtf50_value, double quality_value) {
        assert(edge_number < 4);
        mtf50[edge_number] = mtf50_value;
        quality[edge_number] = quality_value;
    }
    
    double get_mtf50_value(size_t edge_number) const {
        assert(edge_number < 4);
        return mtf50[edge_number];
    }
    
    double get_mtf50_value(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return get_mtf50_value(it->second);
    }
    
    Point get_centroid(void) const {
        return centroid;
    }
    
    double get_area(void) const {
        return area;
    }
    
    double get_quality(size_t edge_number) const {
        return quality[edge_number];
    }
    
    double get_quality(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return quality[it->second];
    }
    
    
    Mrectangle rect;
    vector<double> mtf50;
    vector<double> quality;
    map<edge_position, size_t> edge_lut;
    Point centroid;
    double area;
};

#endif
