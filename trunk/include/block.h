#ifndef BLOCK_H
#define BLOCK_H


#include "rectangle.h"

#include <map>
using std::map;

class Block {
  public:
    typedef enum {TOP, LEFT, RIGHT, BOTTOM} edge_position;
    
    Block(const Rectangle& in_rect) : rect(in_rect), mtf50(4,0.0), centroid(0,0) {
    
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
    
    
    
    Point get_edge_centroid(size_t edge_number) const {
        assert(edge_number < 4);
        return rect.centroids[edge_number];
    }
    
    Point get_edge_centroid(edge_position ep) const {
        map<edge_position, size_t>::const_iterator it = edge_lut.find(ep);
        return rect.centroids[it->second];
    }
    
    void set_mtf50_value(size_t edge_number, double mtf50_value) {
        assert(edge_number < 4);
        mtf50[edge_number] = mtf50_value;
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
    
    
    Rectangle rect;
    vector<double> mtf50;
    map<edge_position, size_t> edge_lut;
    Point centroid;
    double area;
};

#endif
