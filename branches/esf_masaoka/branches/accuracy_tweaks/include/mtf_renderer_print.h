#ifndef MTF_RENDERER_PRINT_H
#define MTF_RENDERER_PRINT_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_print : public Mtf_renderer {
  public:
    Mtf_renderer_print(const std::string& fname, bool filter=false, double angle=0) 
      :  ofname(fname), filter(filter), angle(angle) {
      
    }
    
    void render(const vector<Block>& blocks) {
        FILE* fout = fopen(ofname.c_str(), "wt");
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (filter) {
                    double ba = blocks[i].get_edge_angle(k);
                    double ad = acos(cos(angle)*cos(ba) + sin(angle)*sin(ba));
                    if (fabs(ad) < 5.0/180.0*M_PI || fabs(ad - M_PI) < 5.0/180.0*M_PI) {
                        fprintf(fout, "%lf ", val);
                    }
                } else {
                    fprintf(fout, "%lf ", val);
                }
            }
            fprintf(fout, "\n");
        }    
        fclose(fout);
    }
    
    string ofname;
    bool filter;
    double angle;
};

#endif
