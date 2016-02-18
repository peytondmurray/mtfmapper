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
#include "mtf_profile_sample.h"

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
        const cv::Mat& img, const vector<cv::Point3d>& distance_scale,
        bool lpmm_mode=false, double pixel_size=1.0,
        int largest_block_index=-1) 
      :  zero(zero), transverse(transverse), longitudinal(longitudinal),
         wdir(wdir), prname(prof_fname),
         img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         chart_scale(chart_scale), distance_scale(distance_scale),
         largest_block_index(largest_block_index) {
      
    }
    
    
    
    void render(const vector<Block>& blocks) {
        Point centroid(0,0);
        
        if (blocks.size() < 10) { // probably not a valid image for profiles
            return;
        }
        
        
        
        vector< Mtf_profile_sample > points;
        double max_trans = 0;
        vector<double> mtf50;
        for (size_t i=0; i < blocks.size(); i++) {
            const double angle_thresh = 25.0;
            
            if (int(i) == largest_block_index) {
                continue; // ignore the largest block
            }
            
            double diam = 0;
            for (size_t k=0; k < 3; k++) {
                for (size_t j=k+1; j < 4; j++) {
                    double d = norm(blocks[i].get_edge_centroid(k) - blocks[i].get_edge_centroid(j));
                    diam = std::max(diam, d);
                }
            }
            
            double maxmtf = 0;
            for (size_t k=0; k < 4; k++) {
                Point norm = blocks[i].get_normal(k);
                double delta = longitudinal.x*norm.x + longitudinal.y*norm.y;
                
                
                
                if (acos(fabs(delta))/M_PI*180.0 < angle_thresh && // edge roughly aligned with long axis
                    blocks[i].get_mtf50_value(k) < 1) {            // and not a problematic edge
                    
                    maxmtf = std::max(blocks[i].get_mtf50_value(k), maxmtf);
                }
            }
            
            for (size_t k=0; k < 4; k++) {
                Point ec = blocks[i].get_edge_centroid(k);
                Point norm = blocks[i].get_normal(k);
                double delta = longitudinal.x*norm.x + longitudinal.y*norm.y;
                
                if (acos(fabs(delta))/M_PI*180.0 < angle_thresh &&
                    blocks[i].get_mtf50_value(k) == maxmtf ) {
                
                    for (int ti=-1; ti <= 1; ti++) {
                        Point mec = ec + ti*diam/4.0*Point(-norm.y, norm.x);
                        points.push_back(Mtf_profile_sample(mec, maxmtf, blocks[i].get_edge_angle(k)));
                    }
                    
                    if (fabs(ec.y - zero.y) > max_trans) {
                        max_trans = fabs(ec.y - zero.y);
                    }
                    mtf50.push_back(blocks[i].get_mtf50_value(k));
                } 
            }
        }
        
        if (mtf50.size() == 0) {
            printf("No usable edges found. Aborting manual focus peak evaluation.\n");
            return;
        }
        
        sort(mtf50.begin(), mtf50.end());
        double p2 = mtf50[0.02*mtf50.size()];
        double p98 = mtf50[0.98*mtf50.size()];
        
        vector< Mtf_profile_sample > ap_points; // axis projected points
        for (size_t i=0; i < points.size(); i++) {
            Point d = points[i].p - zero;
            Point coord(
                d.x*longitudinal.x + d.y*longitudinal.y,
                d.x*transverse.x + d.y*transverse.y
            );
            double val = softclamp(points[i].mtf, p2, p98);
            
            // TODO: we can clip the values here to [0,1]
            if (!(isnan(coord.x) || isnan(coord.y) || isnan(val))) {
                ap_points.push_back(Mtf_profile_sample(coord, val, points[i].angle));
            } else {
                printf("element %d has nans: %lf, %lf, %lf\n", (int)i, coord.x, coord.y, val);
            }
        }
        
        
        double wsum = 0;
        double covxx = 0;
        double covxy = 0;
        double covyy = 0;
        double cpx = 0;
        double cpy = 0;
        const double w = 1.0;
        for (size_t i=0; i < ap_points.size(); i++) {
            double temp = w + wsum;
            double delta_x = ap_points[i].p.x - cpx;
            double delta_y = ap_points[i].p.y - cpy;
            double rx = delta_x * w / temp;
            double ry = delta_y * w / temp;
            cpx += rx;
            cpy += ry;
            
            if (isnan(rx) || isnan(ry) || isnan(delta_x) || isnan(delta_y)) {
                printf("element (i=%d) : (%lf, %lf) produced nans: %lf %lf %lf %lf\n",
                    (int)i, ap_points[i].p.x, ap_points[i].p.y, 
                    rx, ry, delta_x, delta_y
                );
            }
                
            covxx += wsum * delta_x * rx;
            covyy += wsum * delta_y * ry;
            covxy += wsum * delta_x * ry;
                
            wsum = temp;
        }
        covxx /= wsum;
        covxy /= wsum;
        covyy /= wsum;
        covxx = sqrt(covxx);
        covyy = sqrt(covyy);
        vector< Mtf_profile_sample > ap_keep;
        for (size_t i=0; i < ap_points.size(); i++) {
            if (fabs(ap_points[i].p.x - cpx) > covxx*2 ||
                fabs(ap_points[i].p.y - cpy) > covyy*1.5) {
            } else {
                ap_keep.push_back(ap_points[i]);
            }
        }
        
        Focus_surface pf(ap_keep, 3, 2, chart_scale, distance_scale);
        
        cv::Mat channel(img.rows, img.cols, CV_8UC1);
        double imin;
        double imax;
        cv::minMaxLoc(img, &imin, &imax);
        img.convertTo(channel, CV_8U, 255.0/(imax - imin), 0);
        
        vector<cv::Mat> channels;
        channels.push_back(channel);
        channels.push_back(channel);
        channels.push_back(channel);
        cv::Mat merged;
        merge(channels, merged);
        int initial_rows = merged.rows;
        merged.resize(merged.rows + 100);
        
        for (size_t i=1; i < points.size(); i+=3) {
            int baseline = 0;
            char buffer[1024];
            sprintf(buffer, "%03d", (int)lrint(points[i].mtf*1000));
            cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, 0.5, 1, &baseline);
            cv::Point to(-ts.width/2,  ts.height/2);
            to.x += points[i].p.x;
            to.y += points[i].p.y;
            cv::putText(merged, buffer, to, 
                cv::FONT_HERSHEY_SIMPLEX, 0.5, 
                CV_RGB(20, 20, 20), 2.5, CV_AA
            );
            cv::putText(merged, buffer, to, 
                cv::FONT_HERSHEY_SIMPLEX, 0.5, 
                CV_RGB(50, 255, 255), 1, CV_AA
            );
        }
        
        draw_curve(merged, pf.ridge_peaks, cv::Scalar(30, 30, 255), 2, true);
        draw_curve(merged, pf.ridge, cv::Scalar(30, 255, 30), 3);
        draw_curve(merged, pf.ridge_p05, cv::Scalar(100, 100, 200), 1);
        draw_curve(merged, pf.ridge_p95, cv::Scalar(100, 100, 200), 1);
        
       
        
        rectangle(merged, Point(0, initial_rows), Point(merged.cols, merged.rows), cv::Scalar::all(255), CV_FILLED);
        
        int font = cv::FONT_HERSHEY_DUPLEX; 
        char tbuffer[1024];
        sprintf(tbuffer, "Focus peak at depth %.1lf mm [%.1lf,%.1lf] relative to chart origin", pf.focus_peak, pf.focus_peak_p05, pf.focus_peak_p95);
        cv::putText(merged, tbuffer, Point(50, initial_rows + (merged.rows-initial_rows)/2), font, 1, cv::Scalar::all(0), 1, CV_AA);
        
        
        imwrite(wdir + prname, merged);
        
    }

  private:
  
    void draw_curve(cv::Mat& image, const vector<Point>& data, cv::Scalar col, double width, bool points=false) {
        int prevx = 0;
        int prevy = 0;
        for (size_t i=0; i < data.size(); i++) {
            double px = data[i].x*longitudinal.x + data[i].y*longitudinal.y + zero.x;
            double py = data[i].x*transverse.x + data[i].y*transverse.y + zero.y;
            
            int ix = lrint(px);
            int iy = lrint(py);
            if (ix >= 0 && ix < image.cols &&
                iy >= 0 && iy < image.rows && i > 0) {
                
                if (points) {
                    cv::line(image, Point(ix, iy), Point(ix, iy), col, width, CV_AA);
                } else {
                    cv::line(image, Point(prevx, prevy), Point(ix, iy), col, width, CV_AA);
                }
            }
            
            prevx = ix;
            prevy = iy;
        }
    }
  
    double softclamp(double x, double lower, double upper, double p=0.98) {
        double s = (x - lower) / (upper - lower);
        if (s > p) {
            return 1.0/(1.0 + exp(-3.89182*s));
        }
        return s < 0 ? 0 : s;
    }

    Point zero;
    Point transverse;
    Point longitudinal;
    std::string wdir;
    std::string prname;
    const cv::Mat& img;
    bool    lpmm_mode;
    double  pixel_size;
    double  chart_scale;
    const vector<cv::Point3d>& distance_scale;
    int largest_block_index;
};

#endif
