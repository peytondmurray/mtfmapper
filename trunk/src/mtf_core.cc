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
    for (size_t k=0; k < 4; k++) {
        Point mid_dir = average_dir(g, int(centroids[k].x), int(centroids[k].y));
        
        // now construct buffer around centroid, aligned with direction, of width max_dot
        Mrectangle nr(rrect, k, max_dot);
        
        map<int, scanline> scanset;
        for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
            for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                Point p(x,y);
                if (nr.is_inside(p)) {
                
                    int iy = lrint(y);
                    int ix = lrint(x);
                    if (iy > 0 && iy < img.rows && ix > 0 && ix < img.cols) {

                        edge_record[k].add_point(ix, iy, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));

                        map<int, scanline>::iterator it = scanset.find(iy);
                        if (it == scanset.end()) {
                            scanline sl(ix,ix);
                            scanset.insert(make_pair(iy, sl));
                        }
                        if (ix < scanset[iy].start) {
                            scanset[iy].start = ix;
                        }
                        if (ix > scanset[iy].end) {
                            scanset[iy].end = ix;
                        }
                    }
                }
            }
        }

        edge_record[k].reduce();
        
        double quality = 0;
        Point rgrad;
        vector <double> sfr(SAMPLES_PER_PIXEL, 0);
        double mtf50 = compute_mtf(centroids[k], scanset, edge_record[k], quality, rgrad, sfr);
        
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
    if (er.rsq >= 0.05 && angle_reduce(angle) > 0.5) { 

        vector<double> sum_x(32*4+1, 0);
        vector<double> sum_xx(32*4+1, 0);
        vector<int>    count(32*4+1, 0);
    
        for (double ea=angle-1.0/180.0*M_PI; ea < angle + 1.0/180.0*M_PI; ea += 0.1/180.0*M_PI) {
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
        return 1000;
    }
    
    double* fft_in_buffer = (double*)fftw_malloc(sizeof(double)*2*(FFT_SIZE+2));
    for (size_t i=0; i < 2*(FFT_SIZE+2); i++) {
        fft_in_buffer[i] = 0.0;
    }

    //loess_fit(ordered, fft_in_buffer, FFT_SIZE, -max_dot, max_dot); // loess_fit computes the ESF derivative as part of the fitting procedure
    bin_fit(ordered, fft_in_buffer, FFT_SIZE, -max_dot, max_dot); // loess_fit computes the ESF derivative as part of the fitting procedure
    
    fftw_complex* fft_out_buffer = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (FFT_SIZE+1));
    fftw_execute_dft_r2c(plan_forward, fft_in_buffer, fft_out_buffer);

    double quad = angle_reduce(best_angle);
    size_t best_angle_idx = 0;
    
    double ba_closest_dist = 90;
    for (size_t i=0; i < 90; i++) {
        if (fabs(sfr_correction_table[i][1] - quad) < 
            fabs(sfr_correction_table[best_angle_idx][1] - quad)) {
        
            best_angle_idx = i;
            ba_closest_dist = fabs(sfr_correction_table[i][1] - quad);
        }
    }
    best_angle_idx = 0;
    //printf("best angle idx = %d, angle = %lf\n", best_angle_idx, sfr_correction_table[best_angle_idx][0]);

    double n0 = sqrt(SQR(fft_out_buffer[0][0]) + SQR(fft_out_buffer[0][1]));
    vector<double> magnitude(SAMPLES_PER_PIXEL);
    for (int i=0; i < SAMPLES_PER_PIXEL; i++) {
        magnitude[i] = sqrt(SQR(fft_out_buffer[i][0]) + SQR(fft_out_buffer[i][1])) / n0;
    }
    double mtf_area = 0.0;
    for (int i=0; i < SAMPLES_PER_PIXEL/2; i++) {
        mtf_area += magnitude[i];
    }

    // find closest area value
    int mtf_area_idx = 0;
    double min_dist = 1e50;
    for (int i=0; i < 26; i++) {
       double dist = fabs(mtf_area - sfr_correction_table[i*90][0]);
       if (dist < min_dist) {
           min_dist = dist;
           mtf_area_idx = i;
       }
    }
    //printf("min mtf idx = %d, actual area = %lf, table = %lf\n",
    //       mtf_area_idx, mtf_area, sfr_correction_table[mtf_area_idx*90][0]
    //       );

    //const double* base_mtf = sfr_correction_table[best_angle_idx + mtf_area_idx*90] + 2;

    const double base_mtf[64] = {
        1,0.999984486866557,0.999936700753812,0.999856901758387,0.999745860494588,0.999602856469437,0.999428689355937,0.999221662400338,0.998982616027967,0.998710907296687,0.998407450694172,0.998071686034368,0.997702586563836,0.997301752799086,0.99686726211316,0.996398805838579,0.99589774256414,0.995362939956322,0.99479389996896,0.994190792626164,0.993553417359351,0.992881206394793,0.99217322712681,0.991430392056121,0.990652214795646,0.989836768283707,0.988985175310973,0.988097190693832,0.987172313789717,0.986208625669173,0.985207467111431,0.984167657360914,0.983089006548719,0.981969944744486,0.980812427046191,0.979613377043908,0.978373304498,0.977092693659347,0.975769406161666,0.974403836942393,0.972993725889167,0.971542204338157,0.970044320125214,0.968503390887393,0.96691588934413,0.965282591291464,0.963601532467787,0.961873845827568,0.96009790159014,0.95827221579671,0.95639708489526,0.954471521572099,0.952494812629752,0.950466568192422,0.948385048193354,0.94625057028364,0.944060345451452,0.941815519979535,0.939514108197148,0.93715450255377,0.93473749841564,0.932260559235891,0.929723608450891,0.927125128100689
    };
    
    double prev_freq = 0;
    double prev_val  = n0;
    
    double mtf50 = 0;

    bool done = false;
    for (int i=0; i < SAMPLES_PER_PIXEL && !done; i++) {
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
        mtf50 = 2.0/SAMPLES_PER_PIXEL;
    }
    
    mtf50 *= max_dot*2; 

    for (size_t i=0; i < size_t(SAMPLES_PER_PIXEL);  i++) {
        sfr[i] = magnitude[i] / base_mtf[i];
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


