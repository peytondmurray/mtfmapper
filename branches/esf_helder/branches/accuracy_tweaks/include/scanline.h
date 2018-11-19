#ifndef SCANLINE_H
#define SCANLINE_H

class scanline {
public:
    scanline(int in_start=0, int in_end=0) : start(in_start), end(in_end) {}
    
    int start;
    int end;
};

#endif 
