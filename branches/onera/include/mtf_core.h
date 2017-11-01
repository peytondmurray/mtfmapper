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
#ifndef MTF_CORE_H
#define MTF_CORE_H

#include "include/component_labelling.h"
#include "include/gradient.h"
#include "include/block.h"
#include "include/rectangle.h"
#include "include/ellipse.h"
#include "include/edge_record.h"
#include "include/loess_fit.h"
#include "include/afft.h"
#include "include/mtf_profile_sample.h"
#include "include/bayer.h"

#include <algorithm>
#include <map>
using std::map;

typedef vector<Block> block_vector;

// global constants for ESF-fourier MTF method
// TODO: these can be dynamic parameters, with some effort
const double max_dot = 28;
const int SAMPLES_PER_PIXEL = 32;
const size_t FFT_SIZE = int(16)*SAMPLES_PER_PIXEL;
const int NYQUIST_FREQ = FFT_SIZE/16;
const double max_edge_length = 200;

class Mtf_core {
  public:
    

    Mtf_core(const Component_labeller& in_cl, const Gradient& in_g, 
             const cv::Mat& in_img, const cv::Mat& in_bayer_img, std::string bayer_subset)
      : cl(in_cl), g(in_g), img(in_img), bayer_img(in_bayer_img), absolute_sfr(false),
        snap_to(false), snap_to_angle(0), sfr_smoothing(true),
        sliding(false), samples_per_edge(0), roi_width(max_dot) {

        bayer = Bayer::from_string(bayer_subset);
        printf("bayer subset is %d\n", bayer);
      
        for (Boundarylist::const_iterator it=cl.get_boundaries().begin(); it != cl.get_boundaries().end(); ++it) {
            valid_obj.push_back(it->first);
        }
        
        cv::Mat temp;
        in_img.convertTo(temp, CV_8U, 256.0/16384.0);
        cv::cvtColor(temp, od_img, CV_GRAY2RGB);
    }
    
    ~Mtf_core(void) {
        cv::imwrite(string("detections.png"), od_img);
    }
    
    size_t num_objects(void) {
        return valid_obj.size();
    }
    
    void search_borders(const Point2d& cent, int label);
    bool extract_rectangle(const Point2d& cent, int label, Mrectangle& rect);
    double compute_mtf(const Point2d& in_cent, const map<int, scanline>& scanset, 
                       Edge_record& er,
                       double& poor, 
                       vector<double>& sfr, vector<double>& esf, vector<int> skiplist=vector<int>());
    
    vector<Block>& get_blocks(void) {
        // make a copy into an STL container if necessary
        if (detected_blocks.size() == 0) {
            for (map<int,Block>::const_iterator it=shared_blocks_map.begin();
                 it != shared_blocks_map.end(); ++it) {
                
                bool allzero = true;
                for (int k=0; k < 4 && allzero; k++) {
                    if (fabs(it->second.get_mtf50_value(k)) > 1e-6) {
                        allzero = false;
                    }
                }
                
                if (it->second.valid && !allzero) { 
                    detected_blocks.push_back(it->second);
                }
            }
        }
        return detected_blocks;
    }
    
    vector<Mtf_profile_sample>& get_samples(void) {
        return samples;
    }
    
    void set_absolute_sfr(bool val) {
        absolute_sfr = val;
    }
    
    void set_sfr_smoothing(bool val) {
        sfr_smoothing = val;
    }
    
    void set_sliding(bool val) {
        sliding = val;
    }
    
    void set_samples_per_edge(int s) {
        samples_per_edge = s;
    }
    
    void set_snap_angle(double angle) {
        snap_to = true;
        snap_to_angle = angle;
    }
    
    void set_roi_width(double width) {
        roi_width = width;
    }
    
    void process_roi_image(cv::Mat& roi_img) {
        const int n_roi=4;
        const string names[4] = {"Vh", "Hl", "Vl", "Hh"};
        
        FILE* f_mtf = fopen("mtf50.txt", "wt");
        FILE* f_sfr = fopen("sfr.csv", "wt");
        FILE* f_stats = fopen("stats.csv", "wt");
        
        FILE* fesf = fopen("e1raw.txt", "wt");
        
        fprintf(f_stats, "edge,MTF50,centroid_x,centroid_y,edge_angle,relative_orientation\n");
        
        vector< vector<double> > all_sfr;
        
        for (int roi=0; roi < n_roi; roi++) {
        
            int npixels = 0;
            map<int, scanline> scanset;
            Edge_record er;
            for (int row=0; row < roi_img.rows; row++) {
                for (int col=0; col < roi_img.cols; col++) {
                    int value = roi_img.at<uint8_t>(row, col);
                    if (value != 255 && (value & (1 << roi)) != 0) {
                        npixels++;
                        
                        {    
                            er.add_point(col, row, fabs(g.grad_x(col, row)), fabs(g.grad_y(col, row)));
                            
                            if (roi == 0) fprintf(fesf, "%d %d %d\n", col, row, img.at<uint16_t>(row, col));
                        }
                        
                        if (scanset.find(row) == scanset.end()) {
                            scanset[row] = scanline(col,col);
                        }
                        scanset[row].start = std::min(col, scanset[row].start);
                        scanset[row].end   = std::max(col, scanset[row].end);
                        
                        
                    } 
                }
            }
            
            if (npixels > 1) {
                
                er.reduce();
                Point2d cent(er.centroid);
                printf("ER reduce grad estimate: %lf\n", er.angle/M_PI*180);
                
                Point2d mean_grad(cos(er.angle), sin(er.angle));
                Point2d edge_direction(-sin(er.angle), cos(er.angle));
                vector< vector<double> > binned_gradient(max_dot*4+1);
                for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                    int y = it->first;
                    for (int x=it->second.start; x <= it->second.end; ++x) {
                        Point2d d((x) - cent.x, (y) - cent.y);
                        double dot = d.ddot(mean_grad); 
                        
                        int idot = lrint(dot*2) + max_dot;
                        if (idot >= 0 && idot < (int)binned_gradient.size()) {
                            binned_gradient[idot].push_back(g.grad_magnitude(x, y));
                        }
                    }
                }
                for (size_t k=0; k < binned_gradient.size(); k++) {
                    auto& b = binned_gradient[k];
                    sort(b.begin(), b.end());
                }
                vector<int> skiplist;
                for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                    int y = it->first;
                    for (int x=it->second.start; x <= it->second.end; ++x) {
                    
                        Point2d d((x) - cent.x, (y) - cent.y);
                        double dot = d.ddot(mean_grad); 
                        
                        if (fabs(dot) <= 2) {
                            // do nothing
                        } else {
                            int idot = lrint(dot*2) + max_dot;
                            if (idot >= 0 && idot < (int)binned_gradient.size() && binned_gradient[idot].size() > 3) {
                                
                                
                                double upper = binned_gradient[idot][(int)floor(binned_gradient[idot].size()*0.95)];
                                if (binned_gradient[idot].back() - binned_gradient[idot][binned_gradient[idot].size()-2] < 0.001) {
                                    continue;
                                }
                                
                                //printf("upper=%lf, grad=%lf, dot=%lf, idot=%d, #bins=%lu: ", upper, g.grad_magnitude(x, y), dot, idot, binned_gradient[idot].size());
                                
                                if (g.grad_magnitude(x, y) >= upper) {
                                    skiplist.push_back(y*img.cols + x);
                                    // add 8-connected neighbours
                                    // may be too much
                                    /*
                                    for (int dy=-1; dy <= 1; dy++) {
                                        for (int dx=-1; dx <= 1; dx++) {
                                            skiplist.push_back((y+dy)*img.cols + (x+dx));
                                        }
                                    }
                                    */
                                    
                                    // add 4-connected neighbours
                                    // somewhat better ?
                                    /*
                                    skiplist.push_back((y+1)*img.cols + x);
                                    skiplist.push_back((y-1)*img.cols + x);
                                    skiplist.push_back(y*img.cols + x - 1);
                                    skiplist.push_back(y*img.cols + x + 1);
                                    */
                                    
                                    //printf("skip\n");
                                } else {
                                    //printf("keep\n");
                                    
                                }
                            } 
                        }
                    }
                }
                sort(skiplist.begin(), skiplist.end());
                
                // draw the actual pixels used
                for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                    int y = it->first;
                    for (int x=it->second.start; x <= it->second.end; ++x) {
                        int idx = y*img.cols + x;
                        
                        if (!binary_search(skiplist.begin(), skiplist.end(), idx)) {
                            cv::Vec3b& color = od_img.at<cv::Vec3b>(y, x);
                            color[0] = 255;
                            color[1] = 255;
                            color[2] = 0;
                        }
                    }
                }
                
                // iterate the edge record to both centre the edge, as well as refine the edge orientation
                cent = er.centroid;
                Point2d normal(cos(er.angle), sin(er.angle)); 
                scanset.clear();
                er.clear();
                for (int row=0; row < roi_img.rows; row++) {
                    for (int col=0; col < roi_img.cols; col++) {
                        int value = roi_img.at<uint8_t>(row, col);
                        if (value != 255 && (value & (1 << roi)) != 0) {
                        
                            Point2d p(col, row);
                            Point2d d = p - cent;
                            double dot = d.x*normal.x + d.y*normal.y;
                            
                            if (fabs(dot) < 12) {
                                int idx = row*img.cols + col;
                                if (!binary_search(skiplist.begin(), skiplist.end(), idx)) {
                                    er.add_point(col, row, fabs(g.grad_x(col, row)), fabs(g.grad_y(col, row)));
                                }
                            }
                            
                            if (scanset.find(row) == scanset.end()) {
                                scanset[row] = scanline(col,col);
                            }
                            scanset[row].start = std::min(col, scanset[row].start);
                            scanset[row].end   = std::max(col, scanset[row].end);
                        } 
                    }
                }
                er.reduce();
                printf("updated: ER reduce grad estimate: %lf, centroid (%lf,%lf) -> (%lf, %lf)\n", 
                    er.angle/M_PI*180, cent.x, cent.y, er.centroid.x, er.centroid.y
                );
                
                double quality;
                
                vector <double> sfr(NYQUIST_FREQ*2, 0);
                vector <double> esf(FFT_SIZE/2, 0);
                
                double mtf50 = compute_mtf(er.centroid, scanset, er, quality, sfr, esf, skiplist);
                printf("(%lf, %lf) -> mtf50=%lf, quality=%lf\n", cent.x, cent.y, mtf50, quality);
                
                all_sfr.push_back(sfr);
                
                fprintf(f_stats, "%s,%lf,%.2lf,%.2lf,%.3lf,%.3lf\n", 
                    names[roi].c_str(), mtf50,
                    cent.x, cent.y,
                    er.angle/M_PI*180.0,
                    angle_reduce(er.angle)
                );
                
            } else {
                printf("Warning: no pixels found for ROI %d\n", roi);
            }
        }
        
        printf("got %ld ROIs\n", all_sfr.size());
        
        for (size_t roi=0; roi < all_sfr.size(); roi++) {
            fprintf(f_sfr, ",%s%c", names[roi].c_str(), roi == (all_sfr.size()-1) ? '\n' : ',');
        }
        
        for (size_t row=0; row < all_sfr[0].size(); row++) {
            double freq = double(row)/64.0;
            for (size_t roi=0; roi < all_sfr.size(); roi++) {
                fprintf(f_sfr, "%lf,%lf%c", freq, all_sfr[roi][row], roi == (all_sfr.size()-1) ? '\n' : ',');
            }
        }
        
        fclose(fesf);
        fclose(f_mtf);
        fclose(f_sfr);
        fclose(f_stats);
    }
    
    const Component_labeller& cl;
    const Gradient&           g;
    const cv::Mat&            img;
    const cv::Mat&            bayer_img;
    Bayer::bayer_t bayer;
    
    AFFT<512> afft; // FFT_SIZE = 512 ??
    vector<int> valid_obj;
    
    vector<Block> detected_blocks;  
    map<int, Block> shared_blocks_map;
    vector<Point2d> solid_ellipses;
    vector<Ellipse_detector> ellipses;
    
    vector<Mtf_profile_sample> samples;
    
    cv::Mat od_img;
  private:
    bool absolute_sfr;
    bool snap_to;
    double snap_to_angle;
    bool sfr_smoothing;
    bool sliding;
    int samples_per_edge;
    double roi_width;
    
    void process_with_sliding_window(Mrectangle& rrect);
  
    void sample_at_angle(double ea, vector<Ordered_point>& local_ordered, 
        const map<int, scanline>& scanset, const Point2d& cent,
        double& edge_length, vector<int> skiplist=vector<int>()) {

        double max_along_edge = -1e50;
        double min_along_edge = 1e50;
        
        Point2d mean_grad(cos(ea), sin(ea));
        Point2d edge_direction(-sin(ea), cos(ea));
        
        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; ++x) {
                
                
                if (skiplist.size() > 0) {
                    int idx = y*img.cols + x;
                    if (binary_search(skiplist.begin(), skiplist.end(), idx)) {
                        //printf("skipping pixel (r=%d, c=%d)\n", y, x);
                        continue;
                    }
                }

                Point2d d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                double dist_along_edge = d.ddot(edge_direction);
                
                
                if (fabs(dot) < max_dot) {
                    local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                    max_along_edge = std::max(max_along_edge, dist_along_edge);
                    min_along_edge = std::min(min_along_edge, dist_along_edge);
                }
            }
        }
        
        edge_length = max_along_edge - min_along_edge;
        //printf("edge length = %lf\n\n", edge_length);
    }

    inline double bin_at_angle(double ea, const map<int, scanline>& scanset, const Point2d& cent,
        vector<double>& sum_a, vector<double>& sum_q, vector<int>& count, vector<int> skiplist=vector<int>()) {

        Point2d mean_grad(cos(ea), sin(ea));
        Point2d edge_direction(sin(ea), cos(ea));

        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; ++x) {
                
                if (skiplist.size() > 0) {
                    int idx = y*img.cols + x;
                    if (binary_search(skiplist.begin(), skiplist.end(), idx)) {
                        continue;
                    }
                }
            
                Point2d d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                if (fabs(dot) < max_dot) {
                    int idot = lrint(dot*4 + max_dot*4);
                    double data = img.at<uint16_t>(y,x)/256.0;
                    if (idot >= 0 && idot < (int)count.size()) {
                    count[idot]++;
                        double old_sum_a = sum_a[idot];
                        sum_a[idot] += (data - sum_a[idot])/double(count[idot]);
                        sum_q[idot] += (data - old_sum_a)*(data - sum_a[idot]);
                    } 
                }
            }
        }
        double varsum = 0;
        int used = 0;
        for (size_t k=0; k < count.size(); k++) {
            if (count[k] > 2) {
                used++;
                varsum += sum_q[k] / double(count[k]-1);
            }
        }
        varsum *= sum_a.size() / double(used);
        return varsum;
    }
    
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {  
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180; 
        return quad1;
    }

};

#endif
