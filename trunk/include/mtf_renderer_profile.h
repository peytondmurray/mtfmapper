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
    Mtf_renderer_profile(const std::string& wdir, const std::string& prof_fname, 
        const std::string& peak_fname, const std::string& gnuplot_binary,
        const cv::Mat& img) 
      :  wdir(wdir), prname(prof_fname), pfname(peak_fname), 
         gnuplot_binary(gnuplot_binary), img(img), gnuplot_failure(false), 
         gnuplot_warning(true) {
      
    }
    
    
    
    void render(const vector<Block>& blocks) {
        size_t largest_block = 0;
        Point centroid(0,0);
        
        if (blocks.size() < 10) { // probably not a valid image for profiles
            return;
        }
        
        map<int, double> row_max;
        extract_row_maxima(row_max, blocks, false);

        // check distribution: a bimodal distribution probably means the image is rotated
        bool transpose = test_for_bimodal_distribution(row_max);
        if (transpose) {
            printf("bimodal distribution found, taking transpose\n");
            extract_row_maxima(row_max, blocks, true);
        }

        // compute robust maximum
        vector<double> mtf50_values;
        for (map<int,double>::const_iterator it=row_max.begin(); it != row_max.end(); it++) {
            mtf50_values.push_back(it->second);
        }
        sort(mtf50_values.begin(), mtf50_values.end());
        // pick 95% percentile
        double effective_max = mtf50_values[0.98*mtf50_values.size()];
        // adjust it upwards a little
        effective_max /= 0.9; 

        for (size_t i=0; i < blocks.size(); i++) {
            if (blocks[i].get_area() > blocks[largest_block].get_area()) {
                largest_block = i;
            }
            centroid.x += blocks[i].get_centroid().x;
            centroid.y += blocks[i].get_centroid().y;
        }
        
        if (row_max.size() == 0) {
            printf("Warning: All mtf50 values are zero; cannot generate a profile.\n");
            return;
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
        int max_med_filt_coord = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            if (med_filt_mtf[i] >= max_med_filt) {
                max_med_filt = med_filt_mtf[i];
                max_med_filt_coord = it->first;
            }
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

        double peak_shift = 0;
        double ref_edge_position = 0;
        if (transpose) {
            ref_edge_position = blocks[largest_block].get_edge_centroid(peak_idx).x;
        } else {
            ref_edge_position = blocks[largest_block].get_edge_centroid(peak_idx).y;
        }
        peak_shift = ref_edge_position - max_med_filt_coord;

        // The LOESS-filtered curve is already a good approximation of the mean MTF50
        // at any given position along the curve. We can thus estimate mean MTF50 at
        // the position of the reference edge by directly looking up the value in the 
        // smoothed curve.
        double mean_mtf50_at_ref_edge = 0;
        min_dist = 1e50;
        i = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); it++) {
            double dist = fabs(it->first - ref_edge_position);
            if (dist < min_dist) {
                min_dist = dist;
                mean_mtf50_at_ref_edge = med_filt_mtf[i++];
            }
        }

        printf("peak_shift = %lf mtf_at_peak = %lf\n", peak_shift, mean_mtf50_at_ref_edge);
        
        FILE* pffile = fopen((wdir+pfname).c_str(), "wt");
        double peak_mtf50 = blocks[largest_block].get_mtf50_value(peak_idx);
        bool peak_quality_good = true;
        if (!blocks[largest_block].get_quality(peak_idx)) {
            peak_mtf50 = max_med_filt;
            peak_quality_good = false;
        }
        fprintf(pffile, "%lf %lf %lf\n", 
            transpose ? 
                 blocks[largest_block].get_edge_centroid(peak_idx).x :
                 blocks[largest_block].get_edge_centroid(peak_idx).y , 
            peak_mtf50,
            peak_mtf50*3
        );
        fclose(pffile);

		       
        FILE* gpf = fopen( (wdir + string("profile.gnuplot")).c_str(), "wt");
        fprintf(gpf, "set xlab \"column (pixels)\"\n");
        fprintf(gpf, "set ylab \"MTF50 (cyc/pix)\"\n");
        fprintf(gpf, "set term png size 1024, 768\n");
        fprintf(gpf, "set output \"%sprofile_image.png\"\n", wdir.c_str());
        fprintf(gpf, "plot [][0:%lf]\"%s\" u 1:2 t \"MTF50 (c/p) raw\" w p ps 0.25, \"%s\" u 1:3 t \"MTF50 (c/p) smoothed\" w l lw 3, \"%s\" u 1:2 t \"Expected focus point\" w i lc %d lw 3\n", 
            effective_max,
            (wdir+prname).c_str(), (wdir+prname).c_str(), (wdir+pfname).c_str(), peak_quality_good ? 3 : 1);
        fclose(gpf);
        
        char* buffer = new char[1024];
        sprintf(buffer, "%s %sprofile.gnuplot", gnuplot_binary.c_str(), wdir.c_str());
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

    void set_gnuplot_warning(bool gnuplot) {
        gnuplot_warning = gnuplot;
    }
    
    bool gnuplot_failed(void) {
        return gnuplot_failure;
    }

  private:
    static double angular_diff(double a, double b) {
        return acos(cos(a)*cos(b) + sin(a)*sin(b));
    }

    bool test_for_bimodal_distribution(const map<int, double>& data) {
        int minval = data.begin()->first;
        int maxval = data.rbegin()->first;
        int span = maxval - minval;

        int nbins = 60;

        // bin coarser
        vector<double> counts(nbins, 0);
        for (map<int, double>::const_iterator it = data.begin(); it != data.end(); it++) {
            int idx = (it->first - minval)*nbins/(span);
            idx = std::max(0, idx);
            idx = std::min(nbins-1, idx);
            counts[idx] = std::max(it->second, counts[idx]);
        }

        // find two breakpoints in density; we are looking for a -_- shape, which
        // indicates that the chart has been rotated
        double min_cost = 1e50;
        int min_start=0;
        int min_end=0;
        for (int b_start=10; b_start < 2*nbins/3; b_start += 2) {
            for (int b_end=std::max(nbins/2, b_start+4); b_end < nbins-10; b_end += 2) {
                double cost = 
                    IQR(counts,0, b_start) + 
                    10*IQR(counts, b_start, b_end) + 
                    IQR(counts, b_end, nbins);
                if (cost < min_cost) {
                    min_cost = cost;
                    min_start = b_start;
                    min_end = b_end;
                }
            }
        }
        double first_segment = median(counts, 0, min_start);
        double middle_segment = median(counts, min_start, min_end);
        double last_segment = median(counts, min_end, nbins);

        return middle_segment <= first_segment && middle_segment <= last_segment;
    }

    double IQR(const vector<double>& counts, int start, int end) {
        vector<double> sorted;
        start = std::max(0, start);
        end = std::min(counts.size()-1, size_t(end));
        for (size_t i=start; i < size_t(end); i++) {
            sorted.push_back(counts[i]);
        }
        sort(sorted.begin(), sorted.end());
        return sorted[lrint(sorted.size()*0.75)] - sorted[lrint(sorted.size()*0.25)];
    }

    double median(const vector<double>& counts, int start, int end) {
        vector<double> sorted;
        start = std::max(0, start);
        end = std::min(counts.size()-1, size_t(end));
        for (size_t i=start; i < size_t(end); i++) {
            sorted.push_back(counts[i]);
        }
        sort(sorted.begin(), sorted.end());
        return sorted[lrint(sorted.size()*0.5)];
    }

    void extract_row_maxima(map<int,double>& row_max, const vector<Block>& blocks, bool transpose=false) {
        row_max.clear();
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                double angle = blocks[i].get_edge_angle(k);
                double mindiff = 0;

                if (transpose) {
                    mindiff = std::min(angular_diff(angle, 0)/M_PI*180, angular_diff(angle, M_PI)/M_PI*180);
                } else {
                    mindiff = std::min(angular_diff(angle, M_PI/2)/M_PI*180, angular_diff(angle, -M_PI/2)/M_PI*180);
                }

                if (val > 0 && blocks[i].get_quality(k) >= 0.5 && mindiff < 30) {
                    Point cent = blocks[i].get_edge_centroid(k);
                    
                    int y = 0;

                    if (transpose) {
                        y = lrint(cent.x);
                        //fprintf(stderr, "%lf %lf %lf\n", cent.y, cent.x, val);
                    } else {
                        y = lrint(cent.y);
                        //fprintf(stderr, "%lf %lf %lf\n", cent.x, cent.y, val);
                    }
                    
                    map<int, double>::iterator it = row_max.find(y);
                    
                    if (it == row_max.end()) {
                        row_max[y] = val;
                    } else {
                        if (val >= row_max[y]) {
                            row_max[y] = val;
                        }
                    }
                }
            }
        }
    }
    
    std::string wdir;
    std::string prname;
    std::string pfname;
    std::string gnuplot_binary;
    const cv::Mat& img;
    bool gnuplot_failure;
    bool gnuplot_warning;
};

#endif
