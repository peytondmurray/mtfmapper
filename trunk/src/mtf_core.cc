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
#include "include/mtf_tables.h"

// global lock to prevent race conditions on detected_blocks
static tbb::mutex global_mutex;

void Mtf_core::search_borders(const Point& cent, int label) {
    
    Mrectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    if (!valid) {
        // this may be an ellipse. check it ...
        Boundarylist::const_iterator it = cl.get_boundaries().find(label);
        Ellipse_detector e;
        int valid = e.fit(cl, g, it->second, 0, 0, 2);
        if (valid) {
            {
                tbb::mutex::scoped_lock lock(global_mutex);
                ellipses.push_back(e);
                if (e.solid) {
                    solid_ellipses.push_back(Point(e.centroid_x, e.centroid_y));
                }
            }
            for (double theta=0; theta < 2*M_PI; theta += M_PI/720.0) {
                double synth_x = e.major_axis * cos(theta);
                double synth_y = e.minor_axis * sin(theta);
                double rot_x = cos(e.angle)*synth_x - sin(e.angle)*synth_y + e.centroid_x;
                double rot_y = sin(e.angle)*synth_x + cos(e.angle)*synth_y + e.centroid_y;

                // clip to image size, just in case
                rot_x = std::max(rot_x, 0.0);
                rot_x = std::min(rot_x, (double)(od_img.cols-1));
                rot_y = std::max(rot_y, 0.0);
                rot_y = std::min(rot_y, (double)(od_img.rows-1));

                cv::Vec3b& color = od_img.at<cv::Vec3b>(lrint(rot_y), lrint(rot_x));
                color[0] = 255;
                color[1] = 255;
                color[2] = 0;
            }
        }
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
    bool reduce_success = true;
    for (size_t k=0; k < 4; k++) {
        // now construct buffer around centroid, aligned with direction, of width max_dot
        Mrectangle nr(rrect, k, max_dot);
        for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
            for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                Point p(x,y);
                if (nr.is_inside(p)) {
                
                    int iy = lrint(y);
                    int ix = lrint(x);
                    if (iy > 0 && iy < img.rows && ix > 0 && ix < img.cols) {

                        edge_record[k].add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));

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

        reduce_success &= edge_record[k].reduce();
    }
    
    if (!reduce_success) {
        printf("reduce failed, probably not a rectangle/quadrangle\n");
        return;
    }
    
    vector< pair<double, pair<int,int> > > pairs;
    for (size_t k=0; k < 3; k++) {
        for (size_t l=k+1; l < 4; l++) {
            double rel_sim = edge_record[k].relative_orientation(edge_record[l]);
            if (rel_sim > 0.80902) { // = cos(M_PI/5.0) ~ 36 degrees
                pairs.push_back( make_pair(1-rel_sim, make_pair(k, l)) );
            }
        }
    }
    sort(pairs.begin(), pairs.end());
    for (size_t i=0; i < pairs.size(); i++) {
        int e1 = pairs[i].second.first;
        int e2 = pairs[i].second.second; 
        if (edge_record[e1].compatible(edge_record[e2])) {
            if (!edge_record[e1].is_pooled() && !edge_record[e2].is_pooled()) {
                Edge_record::pool_edges(edge_record[e1], edge_record[e2]);
            } 
        } 
    }
    
    for (size_t k=0; k < 4; k++) {
        double quality = 0;
        Point rgrad;
        vector <double> sfr(NYQUIST_FREQ*2, 0);
        vector <double> esf(FFT_SIZE/2, 0);
        double mtf50 = compute_mtf(edge_record[k].centroid, scansets[k], edge_record[k], quality, rgrad, sfr, esf);
        
        if (mtf50 <= 1.2) { // reject mtf values above 1.2, since these are impossible, and likely to be erroneous
            tbb::mutex::scoped_lock lock(global_mutex);
            shared_blocks_map[label].set_mtf50_value(k, mtf50, quality);
            shared_blocks_map[label].set_normal(k, rgrad);
            shared_blocks_map[label].set_sfr(k, sfr);
            shared_blocks_map[label].set_esf(k, esf);
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
    
    Mrectangle rrect(main_thetas, thetas, points, g);
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
    Edge_record& er, double& quality, Point& rgrad, 
    vector<double>& sfr, vector<double>& esf) {
    
    quality = 1.0; // assume this is a good edge
    
    Point cent(in_cent);
    
    Point mean_grad(0,0);
   

    double angle = er.angle;
    mean_grad.x = cos(angle);
    mean_grad.y = sin(angle);

    //printf("original angle estimate: %lf %lf\n", angle/M_PI*180, angle_reduce(angle));

    vector<Ordered_point> ordered;
    double best_angle = angle;
    double edge_length = 0;

    // if there appears to be significant noise, refine the edge orientation estimate
    if (er.rsq >= 0.05 && angle_reduce(angle) > 0.5 && angle_reduce(angle) < 44.2 && bayer == NONE) { 
        double min_sum = 1e50;

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
    
    vector<double> fft_out_buffer(FFT_SIZE * 2, 0);

    double SNR = bin_fit(ordered, fft_out_buffer.data(), FFT_SIZE, -max_dot, max_dot, esf); // bin_fit computes the ESF derivative as part of the fitting procedure
    afft.realfft(fft_out_buffer.data());

    double quad = angle_reduce(best_angle);
    
    double n0 = fabs(fft_out_buffer[0]);
    vector<double> magnitude(NYQUIST_FREQ*2+9);
    double sfr_area = 0;
    double prev=0;
    double alpha=0.25;
    for (int i=0; i < NYQUIST_FREQ*2+9; i++) {
        magnitude[i] = sqrt(SQR(fft_out_buffer[i]) + SQR(fft_out_buffer[FFT_SIZE - i])) / n0;
        if (i <= NYQUIST_FREQ) {
            sfr_area += magnitude[i];
        }
        if (i == NYQUIST_FREQ*2 - 2) {
            prev = magnitude[i];
        }
        if (i >= NYQUIST_FREQ*2 - 1) {
            magnitude[i] = prev*(1-alpha) + magnitude[i]*alpha;
            prev = magnitude[i] * 0.5; // mix in some strong decay
        }
    }
    
    
    if (sfr_smoothing) {
        // perform Savitsky-Golay filtering of MTF curve
        // use filters of increasing width
        // narrow filters reduce bias in lower frequencies
        // wide filter perform the requisite strong filtering at high frequencies
        const int sgh = 7;
        vector<double> smoothed(NYQUIST_FREQ*2, 0);
        const double* sgw = 0;
        for (int idx=0; idx < NYQUIST_FREQ*2; idx++) {
            if (idx < sgh) {
                smoothed[idx] = magnitude[idx];
            } else {
                const int stride = 3;
                int filter_order = std::min(5, (idx-5)/stride);
                sgw = savitsky_golay[filter_order];
                for (int x=-sgh; x <= sgh; x++) { 
                    // since magnitude has extra samples at the end, we can safely go past the end
                    smoothed[idx] += magnitude[idx+x] * sgw[x+sgh];
                }
            }
        }
        assert(fabs(magnitude[0] - 1.0) < 1e-6);
        assert(fabs(smoothed[0] - 1.0) < 1e-6);
        for (int idx=0; idx < NYQUIST_FREQ*2; idx++) {
            magnitude[idx] = smoothed[idx]/smoothed[0];
        }
    }

    double* base_mtf = Mtf_correction::get_instance()->w.data();

    double prev_freq = 0;
    double prev_val  = n0;
    
    double mtf50 = 0;

    // first estimate mtf50 using simple linear interpolation
    bool done = false;
    int cross_idx = 0;
    for (int i=0; i < NYQUIST_FREQ*2 && !done; i++) {
        double mag = magnitude[i];
        mag /= base_mtf[i];
        if (prev_val > 0.5 && mag <= 0.5) {
            // interpolate
            double m = -(mag - prev_val)*(FFT_SIZE);
            mtf50 = -(0.5 - prev_val - m*prev_freq) / m;
            cross_idx = i;
            done = true;
        }
        prev_val = mag;
        prev_freq = i / double(FFT_SIZE);
    }
    
    if (sfr_smoothing) {
        // perform least-squares quadratic fit to compute MTF50
        const int hs = 7;
        const int npts = 2*hs + 1;
        if (done && cross_idx >= hs && cross_idx < NYQUIST_FREQ*2-hs-1) {
            const int tdim = 3;
            
            vector< vector<double> > cov(tdim, vector<double>(tdim, 0.0));
            vector<double> b(tdim, 0.0);
            
            vector<double> a(3);
            for (int r=0; r < npts; r++) {
                double y = magnitude[cross_idx-hs+r]/base_mtf[cross_idx-hs+r];
                a[0] = 1;
                a[1] = y;
                a[2] = y*y;
                for (int col=0; col < tdim; col++) { // 0,1,2
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov[col][icol] += a[col]*a[icol];
                    }
                    b[col] += a[col]*(-hs + r); // build rhs of system : A'*b
                }
            }
            
            // hardcode cholesky decomposition
            bool singular = false;
            vector<double> ldiag(tdim, 0.0);
            for (int i=0; i < tdim && !singular; i++) {
                for (int j=i; j < tdim && !singular; j++) {
                    double sum = cov[i][j];
                    for (int k=i-1; k >= 0; k--) {
                        sum -= cov[i][k]*cov[j][k];
                    }
                    if (i == j) {
                        if (sum <= 0.0) {
                            singular = true;
                        } else {
                            ldiag[i] = sqrt(sum);
                        }
                    } else {
                        cov[j][i] = sum/ldiag[i];
                    }
                }
            }
            if (!singular) {
                // hardcode backsubstitution
                vector<double> x(tdim, 0.0);
                for (int i=0; i < tdim; i++) {
                    double sum = b[i];
                    for (int k=i-1; k >= 0; k--) {
                        sum -= cov[i][k]*x[k];
                    }
                    x[i] = sum/ldiag[i];
                }
                for (int i=tdim-1; i >= 0; i--) {
                    double sum = x[i];
                    for (int k=i+1; k < tdim; k++) {
                        sum -= cov[k][i]*x[k];
                    }
                    x[i] = sum/ldiag[i];
                }
            
                double mid = x[0] + 0.5*x[1] + 0.5*0.5*x[2];
                
                mtf50 = (mid + double(cross_idx))/double(FFT_SIZE);
            }
            
        }
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

    if (fabs(quad - 26.565051) < 1) {
        quality = medium_quality;
    }

    if (quad >= 44) {
        quality = poor_quality;
    }
    
    if (edge_length < 25) {  // derate short edges
        quality *= poor_quality;
    }
    
    return mtf50;
}


