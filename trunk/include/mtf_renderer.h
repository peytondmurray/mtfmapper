#ifndef MTF_RENDERER_H
#define MTF_RENDERER_H

#include "common_types.h"
#include "block.h"

class Mtf_renderer {
  public:
    Mtf_renderer(void) {}
    virtual ~Mtf_renderer(void) {}
    
    virtual void render(const vector<Block>& blocks) = 0;
};

#endif
