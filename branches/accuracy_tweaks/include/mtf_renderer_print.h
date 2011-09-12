#ifndef MTF_RENDERER_PRINT_H
#define MTF_RENDERER_PRINT_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_print : public Mtf_renderer {
  public:
    Mtf_renderer_print(const std::string& fname) 
      :  ofname(fname) {
      
    }
    
    void render(const vector<Block>& blocks) {
        FILE* fout = fopen(ofname.c_str(), "wt");
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                fprintf(fout, "%lf ", val);
            }
            fprintf(fout, "\n");
        }    
        fclose(fout);
    }
    
    string ofname;
};

#endif
