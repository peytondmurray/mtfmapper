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
#ifndef MTF_RENDERER_MFPROFILE_H
#define MTF_RENDERER_MFPROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"

#include <stdlib.h>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

/*
#ifndef SQR 
#define SQR(x) ((x)*(x))
#endif
*/

#include "focus_surface.h"

class Mtf_renderer_mfprofile : public Mtf_renderer {
  public:
    Mtf_renderer_mfprofile(
        const Point& zero, const Point& transverse, const Point& longitudinal,
        double chart_scale,
        const std::string& wdir, const std::string& prof_fname, 
        const std::string& gnuplot_binary,
        const cv::Mat& img, const vector<cv::Point3d>& distance_scale,
        bool lpmm_mode=false, double pixel_size=1.0) 
      :  zero(zero), transverse(transverse), longitudinal(longitudinal),
         wdir(wdir), prname(prof_fname),
         gnuplot_binary(gnuplot_binary), img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         chart_scale(chart_scale), distance_scale(distance_scale),
         gnuplot_failure(false), gnuplot_warning(true) {
      
    }
    
    
    
    void render(const vector<Block>& blocks) {
        Point centroid(0,0);
        
        if (blocks.size() < 10) { // probably not a valid image for profiles
            return;
        }
        
        
        
        vector< pair<Point, double> > points;
        double max_trans = 0;
        vector<double> mtf50;
        //FILE* rf = fopen("rawpoints.txt", "wt");
        for (size_t i=0; i < blocks.size(); i++) {
            const double angle_thresh = 20.0;
            
            for (size_t k=0; k < 4; k++) {
                Point ec = blocks[i].get_edge_centroid(k);
                Point norm = blocks[i].get_normal(k);
                double delta = transverse.x*norm.x + transverse.y*norm.y;
                
                if (acos(fabs(delta))/M_PI*180.0 < angle_thresh) { // edge perp to tangent
                    //fprintf(rf, "%lf %lf %lf\n", ec.x, ec.y, blocks[i].get_mtf50_value(k));
                    points.push_back(make_pair(ec, blocks[i].get_mtf50_value(k)));
                    if (fabs(ec.y - zero.y) > max_trans) {
                        max_trans = fabs(ec.y - zero.y);
                    }
                    mtf50.push_back(blocks[i].get_mtf50_value(k));
                } 
            }
        }
        //fclose(rf);
        
        sort(mtf50.begin(), mtf50.end());
        double p5 = mtf50[0.05*mtf50.size()];
        double p95 = mtf50[0.95*mtf50.size()];
        
        FILE* ap = fopen("axis_projected.txt", "wt");
        vector< pair<Point, double> > ap_points; // axis projected points
        for (size_t i=0; i < points.size(); i++) {
            Point d = points[i].first - zero;
            Point coord(
                d.x*longitudinal.x + d.y*longitudinal.y,
                d.x*transverse.x + d.y*transverse.y
            );
            double val = (points[i].second - p5)/(p95-p5);
            
            // TODO: we can clip the values here to [0,1]
            
            ap_points.push_back(make_pair(coord, val));
            fprintf(ap, "%lf %lf %lf\n", coord.x, coord.y, val);
        }
        fclose(ap);
        
        Focus_surface pf(ap_points, 3, 2, chart_scale, distance_scale);

        /*
        FILE* prfile = fopen((wdir+prname).c_str(), "wt");
        i=0;
        double max_med_filt = 0;
        int max_med_filt_coord = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); ++it) {
            if (med_filt_mtf[i] >= max_med_filt) {
                max_med_filt = med_filt_mtf[i];
                max_med_filt_coord = it->first;
            }
            fprintf(prfile, "%lf %lf %lf\n", 
                it->first/pixel_size, 
                it->second*pixel_size, 
                med_filt_mtf[i++]*pixel_size
            );
        }
        
        
        fprintf(prfile, "\n\n%lf %lf %lf\n", 
            transpose ? 
                 blocks[largest_block].get_edge_centroid(peak_idx).x/pixel_size :
                 blocks[largest_block].get_edge_centroid(peak_idx).y/pixel_size , 
            peak_mtf50*pixel_size,
            peak_mtf50*3*pixel_size
        );
        fclose(prfile);
        */
		       
        FILE* gpf = fopen( (wdir + string("mfprofile.gnuplot")).c_str(), "wt");
        fprintf(gpf, "set xlab \"column (%s)\"\n", lpmm_mode ? "mm" : "pixels");
        fprintf(gpf, "set ylab \"MTF50 (%s)\"\n", lpmm_mode ? "line pairs per mm" : "cycles/pixel");
        fprintf(gpf, "set term png size 1024, 768\n");
        fprintf(gpf, "set output \"%sprofile_image.png\"\n", wdir.c_str());
        fprintf(gpf, "plot [][0:%lf]\"%s\" index 0 u 1:2 t \"MTF50 (%s) raw\" w p ps 0.25, \"%s\" index 0 u 1:3 t \"MTF50 (%s) smoothed\" w l lw 3, \"%s\" index 1 u 1:2 t \"Expected focus point\" w i lc %d lw 3\n", 
            1*pixel_size, // TODO: compute effective max (95%?)
            (wdir+prname).c_str(), lpmm_mode ? "lp/mm" : "c/p",
            (wdir+prname).c_str(), lpmm_mode ? "lp/mm" : "c/p",
            (wdir+prname).c_str(), 3); // supposed to be peak quality
        fclose(gpf);
        
        char* buffer = new char[1024];
        #ifdef _WIN32
        sprintf(buffer, "\"\"%s\" \"%smfprofile.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
        #else
        sprintf(buffer, "\"%s\" \"%smfprofile.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
        #endif
        int rval = system(buffer);
        if (rval != 0) {
            printf("Failed to execute gnuplot (error code %d)\n", rval);
            printf("You can try to execute [%s] to render the plots manually\n", buffer);
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

    Point zero;
    Point transverse;
    Point longitudinal;
    std::string wdir;
    std::string prname;
    std::string pfname;
    std::string gnuplot_binary;
    const cv::Mat& img;
    bool    lpmm_mode;
    double  pixel_size;
    double  chart_scale;
    const vector<cv::Point3d>& distance_scale;
    bool gnuplot_failure;
    bool gnuplot_warning;
};

#endif
