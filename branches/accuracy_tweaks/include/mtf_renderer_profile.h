#ifndef MTF_RENDERER_PROFILE_H
#define MTF_RENDERER_PROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"

#include <stdlib.h>

class Mtf_renderer_profile : public Mtf_renderer {
  public:
    Mtf_renderer_profile(const std::string& prof_fname, const std::string& peak_fname, const cv::Mat& img) 
      :  prname(prof_fname), pfname(peak_fname), img(img) {
      
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
        
        // apply some median filtering to remove obvious outliers
        // this unfortunately chops the peaks of the sharp double-exponential
        // curves seen on synthetic images. ah well ...
        const size_t w = 2;
        vector<double> med_filt_mtf(row_max.size(), 0);
        size_t i = 0;
        vector<Ordered_point> ordered;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            vector<double> medwin;
            map<int, double>::const_iterator start = it;
            map<int, double>::const_iterator end = it;
            for (size_t j=0; j < w && start != row_max.begin(); j++, start--);
            for (size_t j=0; j <= w && end != row_max.end(); j++, end++);
            for (; start != end; start++) {
                medwin.push_back(start->second);
            }
            sort(medwin.begin(), medwin.end());
            med_filt_mtf[i++] = medwin[medwin.size()/2];
            double ex = (it->first - img.cols/2)/double(img.cols)*max_dot;
            ordered.push_back(Ordered_point(ex, medwin[medwin.size()/2]));
        }
        
        // apply additional smoothing
        // try various LOESS filter sizes until the sign changes in the slope
        // drops below 5% (which seems to provide relatively oscillation-free curves)
        for (size_t w2=3; w2 < std::max(ordered.size()/10, size_t(6)); w2+=3) {
            for (size_t i=0; i < ordered.size() - 1; i++) {
                Point sol;
                
                size_t start = std::max(0, int(i) - int(w2));
                size_t end   = std::min(ordered.size() - 1, i + w2);
                loess_core(ordered, start, end, ordered[i].first, sol, 1);
                    
                med_filt_mtf[i] = ordered[i].first * sol.y + sol.x;
            }
            double frac_sign_change = 0;
            for (size_t i=3; i < ordered.size()-2; i++) {
                if ( (med_filt_mtf[i] - med_filt_mtf[i-2]) * (med_filt_mtf[i+2] - med_filt_mtf[i]) < 0) {
                    frac_sign_change += 1.0;
                }
            }
            frac_sign_change /= ordered.size();
            if (frac_sign_change < 0.05) {  // less than 5% of sequential slopes experience sign changes
                break;
            }
        }
        
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
        
        FILE* gpf = fopen("profile.gnuplot", "wt");
        fprintf(gpf, "set xlab \"column (pixels)\"\n");
        fprintf(gpf, "set ylab \"MTF50 (cyc/pix)\"\n");
        fprintf(gpf, "set term png size 1024, 768\n");
        fprintf(gpf, "set output \"profile_image.png\"\n");
        fprintf(gpf, "plot \"%s\" u 1:2 t \"MTF50 (c/p) raw\" w p ps 0.25, \"%s\" u 1:3 t \"MTF50 (c/p) smoothed\" w l lw 3, \"%s\" u 1:2 t \"Expected focus point\" w i lw 3\n", 
            prname.c_str(), prname.c_str(), pfname.c_str());
        fclose(gpf);
        
        char* buffer = new char[1024];
        sprintf(buffer, "gnuplot%s profile.gnuplot", EXE_SUFFIX);
        int rval = system(buffer);
        if (rval != 0) {
            printf("Failed to execute gnuplot (error code %d)\n", rval);
            printf("You can try to execute \"%s\" to render the plots manually\n", buffer);
        } else {
            printf("Gnuplot plot completed successfully. Look for profile_image.png\n");
        }
        
        delete [] buffer;
        
    }
    
    std::string prname;
    std::string pfname;
    const cv::Mat& img;
};

#endif
