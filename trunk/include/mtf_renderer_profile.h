#ifndef MTF_RENDERER_PROFILE_H
#define MTF_RENDERER_PROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "gaussfilter.h"

class Mtf_renderer_profile : public Mtf_renderer {
  public:
    Mtf_renderer_profile(const std::string& prof_fname, const std::string& peak_fname) 
      :  prname(prof_fname), pfname(peak_fname) {
      
    }
    
    void render(const vector<Block>& blocks) {
        size_t largest_block = 0;
        Point centroid(0,0);
        
        map<int, double> row_max;
        
        for (size_t i=0; i < blocks.size(); i++) {
            if (blocks[i].get_area() > blocks[largest_block].get_area()) {
                largest_block = i;
            }
            centroid.x += blocks[i].get_centroid().x;
            centroid.y += blocks[i].get_centroid().y;
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0) {
                    Point cent = blocks[i].get_edge_centroid(k);
                    
                    int y = lrint(cent.y);
                    
                    map<int, double>::iterator it = row_max.find(y);
                    
                    if (it == row_max.end()) {
                        row_max[y] = val;
                    } else {
                        if (val > row_max[y]) {
                            row_max[y] = val;
                        }
                    }
                }
            }
        }
        
        const size_t w = 15;
        vector<double> med_filt_mtf(row_max.size(), 0);
        size_t i = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            vector<double> medwin;
            map<int, double>::const_iterator start = it;
            map<int, double>::const_iterator end = it;
            for (size_t j=0; j < w && start != row_max.begin(); j++, start--);
            for (size_t j=0; j < w && end != row_max.end(); j++, end++);
            for (; start != end; start++) {
                medwin.push_back(start->second);
            }
            sort(medwin.begin(), medwin.end());
            med_filt_mtf[i++] = medwin[w];
        }
        
        // apply additional smoothing
        Gaussfilter::filter(med_filt_mtf);
        
        FILE* prfile = fopen(prname.c_str(), "wt");
        i=0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            fprintf(prfile, "%d %lf %lf\n", it->first, it->second, med_filt_mtf[i++]);
        }
        fclose(prfile);
        
        centroid.x -= blocks[largest_block].get_centroid().x;
        centroid.y -= blocks[largest_block].get_centroid().y;
        
        centroid.x /= blocks.size() - 1;
        centroid.y /= blocks.size() - 1;
        
        // profile peak
        size_t peak_idx = 0;
        double min_dist = 1e50;
        for (size_t k=0; k < 4; k++) {
            double dist = sqrt(
                SQR(centroid.x - blocks[largest_block].get_edge_centroid(k).x) +
                SQR(centroid.y - blocks[largest_block].get_edge_centroid(k).y)
            );
            
            if (dist < min_dist) {
                min_dist = dist;
                peak_idx = k;
            }
        }
        
        FILE* pffile = fopen(pfname.c_str(), "wt");
        fprintf(pffile, "%lf %lf %lf\n", 
            blocks[largest_block].get_edge_centroid(peak_idx).y, 
            blocks[largest_block].get_mtf50_value(peak_idx),
            blocks[largest_block].get_mtf50_value(peak_idx)*3
        );
        fclose(pffile);
        
    }
    
    std::string prname;
    std::string pfname;
};

#endif
