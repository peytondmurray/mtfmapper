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
#include "include/edge_record.h"
#include "include/loess_fit.h"

#include <map>
using std::map;

#include "tbb/tbb.h"
using namespace tbb;

typedef vector<Block> block_vector;

#include <fftw3.h>

// global constants for ESF-fourier MTF method
// TODO: these can be dynamic parameters, with some effort
const double max_dot = 16;
const int SAMPLES_PER_PIXEL = 32;
const size_t FFT_SIZE = int(max_dot)*SAMPLES_PER_PIXEL;
const int NYQUIST_FREQ = FFT_SIZE/16;
const double max_edge_length = 100;

class Mtf_core {
  public:
    Mtf_core(const Component_labeller& in_cl, const Gradient& in_g, const cv::Mat& in_img)
      : cl(in_cl), g(in_g), img(in_img) {
      
      
        // set up FFTW plan
        double *fft_in;
        fftw_complex *fft_out;
        fft_in = (double*)fftw_malloc(sizeof(double)*2*(FFT_SIZE+2));
        int nc = (FFT_SIZE)  + 1;
        fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nc);
        plan_forward = fftw_plan_dft_r2c_1d(FFT_SIZE, fft_in, fft_out, FFTW_ESTIMATE);
        fftw_free(fft_out);
        fftw_free(fft_in);
        
        
        for (Boundarylist::const_iterator it=cl.get_boundaries().begin(); it != cl.get_boundaries().end(); it++) {
            valid_obj.push_back(it->first);
        }
    }
    
    ~Mtf_core(void) {
        // clean up FFTW
        fftw_destroy_plan(plan_forward);
    }
    
    size_t num_objects(void) {
        return valid_obj.size();
    }
    
    void search_borders(const Point& cent, int label);
    bool extract_rectangle(const Point& cent, int label, Mrectangle& rect);
    double compute_mtf(const Point& in_cent, const map<int, scanline>& scanset, 
                       Edge_record& er,
                       double& poor, Point& rgrad, vector<double>& sfr);
    
    vector<Block>& get_blocks(void) {
        // make a copy into an STL container if necessary
        if (detected_blocks.size() == 0) {
            for (map<int,Block>::const_iterator it=shared_blocks_map.begin();
                 it != shared_blocks_map.end(); it++) {

                detected_blocks.push_back(it->second);
            }
        }
        return detected_blocks;
    }
    
    const Component_labeller& cl;
    const Gradient&           g;
    const cv::Mat&            img;
    
    // global plan for fourier transform
    fftw_plan plan_forward;
    vector<int> valid_obj;
    
    vector<Block> detected_blocks;  
    map<int, Block> shared_blocks_map;

  private:
    void sample_at_angle(double ea, vector<Ordered_point>& local_ordered, 
        const map<int, scanline>& scanset, const Point& cent,
        double& edge_length) {

        double max_along_edge = -1e50;
        double min_along_edge = 1e50;

        Point mean_grad(cos(ea), sin(ea));
        Point edge_direction(sin(ea), cos(ea));
    
        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                Point d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                double dist_along_edge = d.ddot(edge_direction);
                if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                    local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                    max_along_edge = std::max(max_along_edge, dist_along_edge);
                    min_along_edge = std::min(min_along_edge, dist_along_edge);
                }
            }
        }

        edge_length = max_along_edge - min_along_edge;
    }

    inline double bin_at_angle(double ea, const map<int, scanline>& scanset, const Point& cent,
        vector<double>& sum_a, vector<double>& sum_q, vector<int>& count) {

        Point mean_grad(cos(ea), sin(ea));
        Point edge_direction(sin(ea), cos(ea));

        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                Point d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                double dist_along_edge = d.ddot(edge_direction);
                if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                    int idot = lrint(dot*4 + max_dot*4);
                    double data = img.at<uint16_t>(y,x)/256.0;
                    count[idot]++;
                    double old_sum_a = sum_a[idot];
                    sum_a[idot] += (data - sum_a[idot])/double(count[idot]);
                    sum_q[idot] += (data - old_sum_a)*(data - sum_a[idot]);
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

};

#endif
