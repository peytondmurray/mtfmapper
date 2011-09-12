#ifndef MTF_CORE_H
#define MTF_CORE_H

#include "include/component_labelling.h"
#include "include/gradient.h"
#include "include/block.h"
#include "include/rectangle.h"

#include "tbb/tbb.h"
#include "tbb/mutex.h"
using namespace tbb;

#include <fftw3.h>

// global constants for ESF-fourier MTF method
// TODO: these can be dynamic parameters, with some effort
const double max_dot = 16;
const int SAMPLES_PER_PIXEL = 32;
const size_t FFT_SIZE = int(max_dot)*2*SAMPLES_PER_PIXEL;
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
    double compute_mtf(const Point& in_cent, const map<int, scanline>& scanset);
    
    vector<Block>& get_blocks(void) {
        return detected_blocks;
    }
    
    const Component_labeller& cl;
    const Gradient&           g;
    const cv::Mat&            img;
    
    // global plan for fourier transform
    fftw_plan plan_forward;
    vector<int> valid_obj;
    
    vector<Block> detected_blocks;  
};

#endif
