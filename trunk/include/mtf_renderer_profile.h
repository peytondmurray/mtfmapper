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
#ifndef MTF_RENDERER_PROFILE_H
#define MTF_RENDERER_PROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"

#include <stdlib.h>

class Mtf_renderer_profile : public Mtf_renderer {
  public:
    Mtf_renderer_profile(const std::string& wdir, const std::string& prof_fname, const std::string& peak_fname, const cv::Mat& img) 
      :  wdir(wdir), prname(prof_fname), pfname(peak_fname), 
         img(img), gnuplot_failure(false), gnuplot_warning(true) {
      
    }
    
    void set_gnuplot_warning(bool gnuplot) {
        gnuplot_warning = gnuplot;
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
        for (size_t w2=5; w2 < std::max(ordered.size()/10, size_t(6)); w2+=3) {
            for (size_t i=0; i < ordered.size() - 1; i++) {
                Point sol;
                
                size_t start = std::max(0, int(i) - int(w2));
                size_t end   = std::min(ordered.size() - 1, i + w2);
                loess_core(ordered, start, end, ordered[i].first, sol);
                    
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
        
        FILE* prfile = fopen((wdir+prname).c_str(), "wt");
        i=0;
        double max_med_filt = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            max_med_filt = std::max(max_med_filt, med_filt_mtf[i]);
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
        
        FILE* pffile = fopen((wdir+pfname).c_str(), "wt");
        double peak_mtf50 = blocks[largest_block].get_mtf50_value(peak_idx);
        bool peak_quality_good = true;
        if (!blocks[largest_block].get_quality(peak_idx)) {
            peak_mtf50 = max_med_filt;
            peak_quality_good = false;
        }
        fprintf(pffile, "%lf %lf %lf\n", 
            blocks[largest_block].get_edge_centroid(peak_idx).y, 
            peak_mtf50,
            peak_mtf50*3
        );
        fclose(pffile);

		       
        FILE* gpf = fopen( (wdir + string("profile.gnuplot")).c_str(), "wt");
        fprintf(gpf, "set xlab \"column (pixels)\"\n");
        fprintf(gpf, "set ylab \"MTF50 (cyc/pix)\"\n");
        fprintf(gpf, "set term png size 1024, 768\n");
        fprintf(gpf, "set output \"%sprofile_image.png\"\n", wdir.c_str());
        fprintf(gpf, "plot \"%s\" u 1:2 t \"MTF50 (c/p) raw\" w p ps 0.25, \"%s\" u 1:3 t \"MTF50 (c/p) smoothed\" w l lw 3, \"%s\" u 1:2 t \"Expected focus point\" w i lc %d lw 3\n", 
            (wdir+prname).c_str(), (wdir+prname).c_str(), (wdir+pfname).c_str(), peak_quality_good ? 3 : 1);
        fclose(gpf);
        
        char* buffer = new char[1024];
        sprintf(buffer, "gnuplot%s %sprofile.gnuplot", EXE_SUFFIX, wdir.c_str());
        int rval = system(buffer);
        if (rval != 0) {
            printf("Failed to execute gnuplot (error code %d)\n", rval);
            printf("You can try to execute \"%s\" to render the plots manually\n", buffer);
            gnuplot_failure = true;
        } else {
            printf("Gnuplot plot completed successfully. Look for profile_image.png\n");
        }
        
        delete [] buffer;
        
    }
    
    bool gnuplot_failed(void) {
        return gnuplot_failure;
    }
    
    std::string wdir;
    std::string prname;
    std::string pfname;
    const cv::Mat& img;
    bool gnuplot_failure;
    bool gnuplot_warning;
};

#endif
