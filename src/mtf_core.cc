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
#include "include/mtf_core.h"

#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include "include/peak_detector.h"

#include "include/point_helpers.h"
#include "include/mtf50_edge_quality_rating.h"
#include "include/sfr_table.h"


// global lock to prevent race conditions on detected_blocks
static tbb::mutex global_mutex;

void Mtf_core::search_borders(const Point& cent, int label) {
    
    Mrectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    if (!valid) {
        return;
    }
    
    Block block(rrect);

    if (block.get_area() > 225) {
        tbb::mutex::scoped_lock lock(global_mutex);
        shared_blocks_map[label] = block;
    } else {
        return;
    }
    
    vector<Point>& centroids = rrect.centroids;
    vector<Edge_record> edge_record(4);
    vector< map<int, scanline> > scansets(4); 
    for (size_t k=0; k < 4; k++) {
        Point mid_dir = average_dir(g, int(centroids[k].x), int(centroids[k].y));
        
        // now construct buffer around centroid, aligned with direction, of width max_dot
        Mrectangle nr(rrect, k, max_dot);
        for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
            for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                Point p(x,y);
                if (nr.is_inside(p)) {
                
                    int iy = lrint(y);
                    int ix = lrint(x);
                    if (iy > 0 && iy < img.rows && ix > 0 && ix < img.cols) {

                        edge_record[k].add_point(ix, iy, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));

                        map<int, scanline>::iterator it = scansets[k].find(iy);
                        if (it == scansets[k].end()) {
                            scanline sl(ix,ix);
                            scansets[k].insert(make_pair(iy, sl));
                        }
                        if (ix < scansets[k][iy].start) {
                            scansets[k][iy].start = ix;
                        }
                        if (ix > scansets[k][iy].end) {
                            scansets[k][iy].end = ix;
                        }
                    }
                }
            }
        }

        edge_record[k].reduce();
    }

    size_t first_edge = 0;
    vector<int> orientation_a;
    vector<int> orientation_b;
    for (size_t k=1; k < 4; k++) {
        if (edge_record[k].get_orientation() == edge_record[first_edge].get_orientation()) {
            orientation_a.push_back(first_edge);
            orientation_a.push_back(k);
        } else {
            orientation_b.push_back(k);
        }
    }
    // if we have two clear pairs of matched edges, attempt to merge them
    if (orientation_a.size() == 2 && orientation_b.size() == 2) {
        if (edge_record[orientation_a[0]].compatible(edge_record[orientation_a[1]])) {
            Edge_record merged(edge_record[orientation_a[0]], edge_record[orientation_a[1]]);
            edge_record[orientation_a[0]].set_angle_from_slope(merged.slope, true);
            edge_record[orientation_a[1]].set_angle_from_slope(merged.slope, true);
        }
        if (edge_record[orientation_b[0]].compatible(edge_record[orientation_b[1]])) {
            Edge_record merged(edge_record[orientation_b[0]], edge_record[orientation_b[1]]);
            edge_record[orientation_b[0]].set_angle_from_slope(merged.slope, true);
            edge_record[orientation_b[1]].set_angle_from_slope(merged.slope, true);
        }
    }
    
    for (size_t k=0; k < 4; k++) {
        double quality = 0;
        Point rgrad;
        vector <double> sfr(NYQUIST_FREQ*2, 0);
        double mtf50 = compute_mtf(edge_record[k].centroid, scansets[k], edge_record[k], quality, rgrad, sfr);
        
        if (mtf50 <= 1.2) { // reject mtf values above 1.2, since these are impossible, and likely to be erroneous
            tbb::mutex::scoped_lock lock(global_mutex);
            shared_blocks_map[label].set_mtf50_value(k, mtf50, quality);
            shared_blocks_map[label].set_normal(k, rgrad);
            shared_blocks_map[label].set_sfr(k, sfr);
        }
    }
    
}

bool Mtf_core::extract_rectangle(const Point& cent, int label, Mrectangle& rect) {
    
    int ix = lrint(cent.x);
    int iy = lrint(cent.y);
    
    // skip non-convex objects with centroids outside of the image 
    if (ix < 0 || ix > cl.get_width() || iy < 0 || iy > cl.get_height()) {
        return false;
    }
    
    // catch some of the non-convex object that might slip through
    if (cl(ix,iy) != label) {
        return false;
    }
    
    Pointlist points = cl.get_boundaries().find(label)->second;
    
    vector<double> thetas(points.size(), 0);
    for (size_t i=0; i < points.size(); i++) { 
        Point dir = average_dir(g, lrint(points[i].x), lrint(points[i].y));
        thetas[i] = atan2(-dir.x, dir.y); // funny ordering and signs because of average_dir conventions
    }
    vector<double> main_thetas(4,0.0);
    
    Peak_detector pd(thetas, 360/5.0);
    pd.select_best_n(main_thetas, 4);
    sort(main_thetas.begin(), main_thetas.end());
    
    Mrectangle rrect(main_thetas, thetas, points);
    rect = rrect;
    
    return rrect.valid;
}

static double angle_reduce(double x) {
    double quad1 = fabs(fmod(x, M_PI/2.0));
    if (quad1 > M_PI/4.0) {
        quad1 = M_PI/2.0 - quad1;
    }
    quad1 = quad1 / M_PI * 180;
    return quad1;
}

double Mtf_core::compute_mtf(const Point& in_cent, const map<int, scanline>& scanset,
    Edge_record& er, double& quality, Point& rgrad, vector<double>& sfr) {
    quality = 1.0; // assume this is a good edge
    
    Point cent(in_cent);
    
    Point mean_grad(0,0);
   

    double angle = er.angle;
    mean_grad.x = cos(angle);
    mean_grad.y = sin(angle);

    //printf("original angle estimate: %lf %lf\n", angle/M_PI*180, angle_reduce(angle));

    vector<Ordered_point> ordered;
    double min_sum = 1e50;
    double best_angle = angle;
    double edge_length = 0;

    // if there appears to be significant noise, refine the edge orientation estimate
    if (er.rsq >= 0.05 && angle_reduce(angle) > 0.5 && angle_reduce(angle) < 44.2 && bayer == NONE) { 

        vector<double> sum_x(32*4+1, 0);
        vector<double> sum_xx(32*4+1, 0);
        vector<int>    count(32*4+1, 0);

        double span = 1.0/180.0*M_PI;
        double step = 0.1/180.0*M_PI;

        if (er.is_pooled()) {
            span /= 3;
            step /= 3;
        }
    
        for (double ea=angle-span; ea < angle + span; ea += step) {
            for (size_t k=0; k < sum_x.size(); k++) {
                sum_x[k]  = 0;
                sum_xx[k] = 0;
                count[k]  = 0;
            }
            double varsum = bin_at_angle(ea, scanset, cent, sum_x, sum_xx, count);
            if (varsum < min_sum || (varsum == min_sum && fabs(ea-angle) < fabs(best_angle-angle)) ) {
                min_sum = varsum;
                best_angle = ea;
            }
        }
    }
    
    //printf("optimized angle estimate: %lf %lf\n", best_angle/M_PI*180, angle_reduce(best_angle));
    sample_at_angle(best_angle, ordered, scanset, cent, edge_length);
    sort(ordered.begin(), ordered.end());
    
    mean_grad.x = cos(best_angle);
    mean_grad.y = sin(best_angle);
    rgrad = mean_grad;
    
    if (ordered.size() < 10) {
        quality = 0; // this edge is not usable in any way
        return 0;
    }
    
    double* fft_in_buffer = (double*)fftw_malloc(sizeof(double)*2*(FFT_SIZE+2));
    for (size_t i=0; i < 2*(FFT_SIZE+2); i++) {
        fft_in_buffer[i] = 0.0;
    }

    double SNR = bin_fit(ordered, fft_in_buffer, FFT_SIZE, -max_dot, max_dot); // loess_fit computes the ESF derivative as part of the fitting procedure
    
    fftw_complex* fft_out_buffer = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (FFT_SIZE+1));
    fftw_execute_dft_r2c(plan_forward, fft_in_buffer, fft_out_buffer);

    double quad = angle_reduce(best_angle);
    
    double n0 = sqrt(SQR(fft_out_buffer[0][0]) + SQR(fft_out_buffer[0][1]));
    vector<double> magnitude(NYQUIST_FREQ*2);
    double sfr_area = 0;
    for (int i=0; i < NYQUIST_FREQ*2; i++) {
        magnitude[i] = sqrt(SQR(fft_out_buffer[i][0]) + SQR(fft_out_buffer[i][1])) / n0;
        if (i <= NYQUIST_FREQ) {
            sfr_area += magnitude[i];
        }
    }

    // find closest entry in sfr correction table
    const double* base_mtf = sfr_correction_table[22] + 2;
    double cdist = 1e50;
    int min_idx = 0;
    for (int i=0; i < 27; i++) {
        double dist = fabs(sfr_correction_table[i][0] - sfr_area);
        if (dist < cdist) {
            min_idx = i;
            cdist = dist;
        }
    }
    base_mtf = sfr_correction_table[min_idx] + 2;

    // critical angles: 14.036 and 26.565 -> these are 0.25 and 0.5 respectively, so they result in 
    // severe quantization of the ESF, leading to underestimates of the MTF 
    // should we detect and correct for this?

    double prev_freq = 0;
    double prev_val  = n0;
    
    double mtf50 = 0;

    bool done = false;
    for (int i=0; i < NYQUIST_FREQ*2 && !done; i++) {
        double mag = magnitude[i];
        mag /= base_mtf[i];
        if (prev_val > 0.5 && mag <= 0.5) {
            // interpolate
            double m = -(mag - prev_val)*(FFT_SIZE);
            mtf50 = -(0.5 - prev_val - m*prev_freq) / m;
            done = true;
        }
        prev_val = mag;
        prev_freq = i / double(FFT_SIZE);
    }
    if (!done) {
        mtf50 = 0.125;
    }
    mtf50 *= 8;

    if (absolute_sfr) {
        for (size_t i=0; i < size_t(NYQUIST_FREQ*2);  i++) {
            sfr[i] = (n0*magnitude[i] / base_mtf[i])/(65536*2);
        }
    } else {
        for (size_t i=0; i < size_t(NYQUIST_FREQ*2);  i++) {
            sfr[i] = magnitude[i] / base_mtf[i];
        }
    }

    // derate the quality of the known poor angles
    if (quad <= 1) {
        quality = poor_quality;
    }

    if (fabs(quad - 26) < 1.5) {
        quality = poor_quality;
    }

    if (quad >= 44) {
        quality = poor_quality;
    }
    
    if (edge_length < 25) {  // derate short edges
        quality *= poor_quality;
    }
    
    fftw_free(fft_out_buffer);
    fftw_free(fft_in_buffer);        
    
    return mtf50;
}


