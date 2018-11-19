#include "include/common_types.h"
#include "include/gradient.h"

Point centroid(const Pointlist& p) {
    double mx = 0;
    double my = 0;
    
    for (size_t i=0; i < p.size(); i++) {
        mx += p[i].x;
        my += p[i].y;
    }
    mx /= p.size();
    my /= p.size();
    
    return Point(mx, my);
}

Point average_dir(const Gradient& g, int x, int y) {
    double mx = 0;
    double my = 0;
    
    double wsum = 0;
    
    int offsets[9][2] = { {0,0}, {-1,0}, {1,0}, {0,-1}, {0,1},
                                 {-1,1}, {1,1}, {1,-1}, {-1,-1} };
    
    for (int k=0; k < 9; k++) {
    
        int lx = x + offsets[k][0];
        int ly = y + offsets[k][1];
    
        mx += g.grad_x(lx,ly) * g.grad_magnitude(lx,ly);
        my += g.grad_y(lx,ly) * g.grad_magnitude(lx,ly);
        wsum += g.grad_magnitude(lx,ly);
        
    }
    
    mx /= wsum;
    my /= wsum;
    return Point(-my, mx);
}

Point normalize(const Point& p) {
    Point q;
    double norm = sqrt(p.ddot(p));
    q.x = p.x / norm;
    q.y = p.y / norm;
    return q; 
}
