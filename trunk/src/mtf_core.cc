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
#include "include/logger.h"
#include "include/mtf_core.h"

#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include "include/peak_detector.h"

#include "include/point_helpers.h"
#include "include/mtf50_edge_quality_rating.h"
#include "include/mtf_tables.h"
#include "include/ellipse_decoder.h"

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <mutex>

// global lock to prevent race conditions on detected_blocks
static std::mutex global_mutex;

void Mtf_core::search_borders(const Point2d& cent, int label) {
    
    Mrectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    if (!valid && find_fiducials) {
        // this may be an ellipse. check it ...
        Boundarylist::const_iterator it = cl.get_boundaries().find(label);
        Ellipse_detector e;
        int valid = e.fit(cl, g, it->second, 0, 0, 2);
        if (valid) {
            Ellipse_decoder ed(e, img);
            
            // search small region near centre for mismatched labels
            // this is almost a duplicate of the foreground fraction check
            // in the Ellipse_detector class, but adds the constraint that
            // the hole must be near the centre .... 
            // TODO: consolidate fiducial validation in one spot
            bool hole_found = false;
            for (int dy=-3; dy <= 3 && !hole_found; dy++) {
                for (int dx=-3; dx <= 3 && !hole_found; dx++) {
                    if (cl(lrint(e.centroid_x+dx), lrint(e.centroid_y+dy)) != label) {
                        hole_found = true;
                    }
                }
            }
            
            if (ed.valid && ed.code >= 0 && hole_found) {
                {
                    std::lock_guard<std::mutex> lock(global_mutex);
                    e.valid = ed.valid;
                    e.set_code(ed.code);
                    ellipses.push_back(e);
                    logger.debug("Fiducial with code %d extracted at (%.2lf %.2lf)\n", e.code, e.centroid_x, e.centroid_y);
                }
                
                for (double theta=0; theta < 2*M_PI; theta += M_PI/720.0) {
                    double synth_x = e.major_axis * cos(theta);
                    double synth_y = e.minor_axis * sin(theta);
                    double rot_x = cos(e.angle)*synth_x - sin(e.angle)*synth_y + e.centroid_x;
                    double rot_y = sin(e.angle)*synth_x + cos(e.angle)*synth_y + e.centroid_y;

                    // clip to image size, just in case
                    rot_x = max(rot_x, 0.0);
                    rot_x = min(rot_x, (double)(od_img.cols-1));
                    rot_y = max(rot_y, 0.0);
                    rot_y = min(rot_y, (double)(od_img.rows-1));

                    cv::Vec3b& color = od_img.at<cv::Vec3b>(lrint(rot_y), lrint(rot_x));
                    color[0] = 255;
                    color[1] = 255;
                    color[2] = 0;
                    
                    
                }
                
                char buffer[20];
                int baseline = 0;
                sprintf(buffer, "%d", e.code);
                cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, 0.5, 1, &baseline);
                cv::Point to(-ts.width/2,  ts.height/2);
                to.x += e.centroid_x;
                to.y += e.centroid_y;
                
                cv::putText(od_img, buffer, to, 
                    cv::FONT_HERSHEY_SIMPLEX, 0.5, 
                    CV_RGB(20, 20, 20), 2.5, CV_AA
                );
                cv::putText(od_img, buffer, to, 
                    cv::FONT_HERSHEY_SIMPLEX, 0.5, 
                    CV_RGB(0, 255, 255), 1, CV_AA
                );
            }
        } 
        return;
    }
    
    Block block(rrect);

    if (block.get_area() <= 225) {
        return;
    }
    
    if (!rrect.corners_ok()) {
        logger.debug("discarding broken square (early)\n");
        return;
    }
    
    if (sliding) {
        process_with_sliding_window(rrect);
        return;
    }
    
    vector<Edge_record> edge_record(4);
    vector< map<int, scanline> > scansets(4); 
    bool reduce_success = true;
    for (size_t k=0; k < 4; k++) {
        // now construct buffer around centroid, aligned with direction, of width max_dot
        Mrectangle nr(rrect, k, max_dot+0.5);
        
        // TODO: connected-component based masking
        //cv::Mat limg(nr.br.y - nr.tl.y + 1, nr.br.x - nr.tl.x + 1, CV_8UC1, cv::Scalar::all(1));
        //cv::Mat roiimg(nr.br.y - nr.tl.y + 1, nr.br.x - nr.tl.x + 1, CV_8UC1, cv::Scalar::all(0));
        
        for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
            for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                Point2d p(x,y);
                Point2d d = p - rrect.centroids[k];
                double dot = d.x*rrect.normals[k].x + d.y*rrect.normals[k].y;
                if (nr.is_inside(p) && fabs(dot) < 14) {
                
                    int iy = lrint(y);
                    int ix = lrint(x);
                    if (iy >= 0 && iy < img.rows && ix >= 0  && ix < img.cols) {
                        edge_record[k].add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                    }
                    
                    // TODO: connected-component based masking
                    //roiimg.at<uint8_t>(iy - nr.tl.y, ix - nr.tl.x) = 255;
                }
                
                /*
                // TODO: connected-component based masking
                int iy = lrint(y);
                int ix = lrint(x);
                limg.at<uint8_t>(iy - nr.tl.y, ix - nr.tl.x) = cl(ix, iy) > 0 ? 0 : 1;
                */
            }
        }
        
        /*
        // TODO: connected-component based masking
        cv::Mat dest;
        cv::Mat labels;
        cv::distanceTransform(limg, dest, labels, CV_DIST_L1, 3);
        // The cc-labelled image works well, but a scan set is no longer convex
        // scan sets would have to be replaced by actual images ?
        // allocating a new cv::Mat for every ROI will cause a lot of memory allocations
        // but what else can we do?
        cv::imwrite("labels_" + std::to_string(k) + ".png", limg);
        cv::imwrite("dist_" + std::to_string(k) + ".png", labels);
        cv::imwrite("roi_" + std::to_string(k) + ".png", roiimg);
        */
        
        reduce_success &= edge_record[k].reduce();
    }
    
    double max_shift = 0;
    if (reduce_success) {
        // re-calculate the ROI after we have refined the edge centroid above
        Mrectangle newrect(rrect, edge_record);
        if (!newrect.corners_ok()) {
            logger.debug("discarding broken square (after updates)\n");
            return;
        }
        
        
        scansets = vector< map<int, scanline> >(4); // re-initialise
        for (size_t k=0; k < 4; k++) {
            // now construct buffer around centroid, aligned with direction, of width max_dot
            Mrectangle nr(newrect, k, max_dot+0.5);
            edge_record[k].clear();
            
            for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
                for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                    Point2d p(x,y);
                    if (nr.is_inside(p)) {
                    
                        int iy = lrint(y);
                        int ix = lrint(x);
                        if (iy >= 0 && iy < img.rows && ix >= 0  && ix < img.cols) {
                        
                            Point2d d = p - newrect.centroids[k];
                            double dot = d.x*newrect.normals[k].x + d.y*newrect.normals[k].y;
                        
                            if (fabs(dot) < 12) {
                                edge_record[k].add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                            }
                            
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
            Point2d ocx = edge_record[k].centroid;
            reduce_success &= edge_record[k].reduce();
            Point2d ncx = edge_record[k].centroid;
            double shift = sqrt(SQR(ocx.x - ncx.x) + SQR(ocx.y - ncx.y));
            max_shift = max(max_shift, shift);
        }
        rrect = newrect;
    }
    
    if (reduce_success && max_shift > 1) {
        // re-calculate the ROI after we have refined the edge centroid above
        Mrectangle newrect(rrect, edge_record);
        if (!newrect.corners_ok()) {
            logger.debug("discarding broken square (after updates)\n");
            return;
        }
        
        
        scansets = vector< map<int, scanline> >(4); // re-initialise
        for (size_t k=0; k < 4; k++) {
            // now construct buffer around centroid, aligned with direction, of width max_dot
            Mrectangle nr(newrect, k, max_dot+0.5);
            edge_record[k].clear();
            scansets[k].clear();
            
            for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
                for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                    Point2d p(x,y);
                    if (nr.is_inside(p)) {
                    
                        int iy = lrint(y);
                        int ix = lrint(x);
                        if (iy >= 0 && iy < img.rows && ix >= 0  && ix < img.cols) {
                        
                            Point2d d = p - newrect.centroids[k];
                            double dot = d.x*newrect.normals[k].x + d.y*newrect.normals[k].y;
                        
                            if (fabs(dot) < 12) {
                                edge_record[k].add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                            }

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
    }
    
    if (!reduce_success) {
        logger.debug("reduce failed, probably not a rectangle/quadrangle\n");
        return;
    }
    
    #if 1
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
    #endif
    
    bool allzero = true;
    for (size_t k=0; k < 4; k++) {
        double quality = 0;
        vector <double> sfr(mtf_width, 0);
        vector <double> esf(FFT_SIZE/2, 0);
        vector <Point2d> ridge(1);
        double mtf50 = compute_mtf(edge_record[k].centroid, scansets[k], edge_record[k], quality, sfr, esf, ridge);
        
        allzero &= fabs(mtf50) < 1e-6;
        
        if (mtf50 <= 1.2) { // reject mtf values above 1.2, since these are impossible, and likely to be erroneous
            std::lock_guard<std::mutex> lock(global_mutex);
            if (shared_blocks_map.find(label) == shared_blocks_map.end()) {
                shared_blocks_map[label] = block;
            }
            shared_blocks_map[label].set_mtf50_value(k, mtf50, quality);
            shared_blocks_map[label].set_normal(k, Point2d(cos(edge_record[k].angle), sin(edge_record[k].angle)));
            shared_blocks_map[label].set_sfr(k, sfr);
            shared_blocks_map[label].set_esf(k, esf);
            shared_blocks_map[label].set_ridge(k, ridge);
        }
    }
    if (allzero) {
        std::lock_guard<std::mutex> lock(global_mutex);
        auto it = shared_blocks_map.find(label);
        if (it != shared_blocks_map.end()) {
            shared_blocks_map.erase(it);
        }
    }
}

bool Mtf_core::extract_rectangle(const Point2d& cent, int label, Mrectangle& rect) {
    
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
        Point2d dir = average_dir(g, lrint(points[i].x), lrint(points[i].y));
        thetas[i] = atan2(-dir.x, dir.y); // funny ordering and signs because of average_dir conventions
    }
    vector<double> main_thetas(4,0.0);
    
    Peak_detector pd(thetas, 360/5.0);
    pd.select_best_n(main_thetas, 4);
    sort(main_thetas.begin(), main_thetas.end());
    
    Mrectangle rrect(main_thetas, thetas, points, g);
    rect = rrect;
    
    /*
    // TODO: rethink this test
    for (int ci=0; ci < 4; ci++) {
        if (cl(lrint(0.5*(ix+rrect.centroids[0].x)), lrint(0.5*(iy+rrect.centroids[0].y))) != label) {
            logger.debug("block with centroid (%d, %d) failed interior check\n", ix, iy);
            return false;
        }
        if (cl(lrint(0.5*(ix+rrect.corners[0].x)), lrint(0.5*(iy+rrect.corners[0].y))) != label) {
            logger.debug("block with centroid (%d, %d) failed interior (corner) check\n", ix, iy);
            return false;
        }
    }
    */
    
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

double Mtf_core::compute_mtf(const Point2d& in_cent, const map<int, scanline>& scanset,
    Edge_record& er, double& quality,  
    vector<double>& sfr, vector<double>& esf, 
    vector<Point2d>& ridge,
    bool allow_peak_shift) {
    
    quality = 1.0; // assume this is a good edge
    
    Point2d cent(in_cent);
    
    Point2d mean_grad(0,0);
   

    double angle = er.angle;
    mean_grad.x = cos(angle);
    mean_grad.y = sin(angle);

    vector<Ordered_point> ordered;
    double edge_length = 0;

    vector<double> fft_out_buffer(FFT_SIZE * 2, 0);
    
    sample_at_angle(angle, ordered, scanset, cent, edge_length, ridge);
    if (ridges_only) {
        return 0.01; // dummy MTF50 value
    }
    sort(ordered.begin(), ordered.end());
    
    if (ordered.size() < 10) {
        quality = 0; // this edge is not usable in any way
        return 0;
    }
    
    int success = bin_fit(ordered, fft_out_buffer.data(), FFT_SIZE, -max_dot, max_dot, esf, allow_peak_shift); // bin_fit computes the ESF derivative as part of the fitting procedure
    if (success < 0) {
        quality = poor_quality;
        logger.debug("failed edge\n");
        return 1.0;
    }
    afft.realfft(fft_out_buffer.data());

    double quad = angle_reduce(angle);
    
    double n0 = fabs(fft_out_buffer[0]);
    vector<double> magnitude(NYQUIST_FREQ*4);
    for (int i=0; i < NYQUIST_FREQ*4; i++) {
        magnitude[i] = sqrt(SQR(fft_out_buffer[i]) + SQR(fft_out_buffer[FFT_SIZE - i])) / n0;
    }
    
    
    if (sfr_smoothing) {
        // apply narrow SG filter to lower frequencies
        vector<double> slf(7, 0);
        const double lf_sgw[5] = {-0.086, 0.343, 0.486, 0.343, -0.086};
        for (int idx=0; idx < 7; idx++) {
            for (int x=-2; x <= 2; x++) {
                slf[idx] += magnitude[abs(idx+x)] * lf_sgw[x+2];
            }
        }
        for (int idx=0; idx < 7; idx++) {
            magnitude[idx] = slf[idx];
        }
        
        // perform Savitsky-Golay filtering of MTF curve
        // use filters of increasing width
        // narrow filters reduce bias in lower frequencies
        // wide filter perform the requisite strong filtering at high frequencies
        const int sgh = 7;
        vector<double> smoothed(NYQUIST_FREQ*4, 0);
        const double* sgw = 0;
        for (int idx=0; idx < NYQUIST_FREQ*4 - sgh; idx++) {
            if (idx < sgh) {
                smoothed[idx] = magnitude[idx];
            } else {
                const int stride = 3;
                int filter_order = min(5, (idx-5)/stride);
                sgw = savitsky_golay[filter_order];
                for (int x=-sgh; x <= sgh; x++) { 
                    // since magnitude has extra samples at the end, we can safely go past the end
                    smoothed[idx] += magnitude[idx+x] * sgw[x+sgh];
                }
            }
        }
        for (int idx = NYQUIST_FREQ * 4 - sgh; idx < NYQUIST_FREQ * 4; idx++) {
            smoothed[idx] = magnitude[idx];
        }
        assert(fabs(magnitude[0] - 1.0) < 1e-6);
        assert(fabs(smoothed[0] - 1.0) < 1e-6);
        for (int idx=0; idx < NYQUIST_FREQ*4; idx++) {
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
    
    if (!done || fabs(quad) < 0.1) {
        mtf50 = 0.125;
    }
    mtf50 *= 8;

    if (absolute_sfr) {
        for (size_t i=0; i < sfr.size();  i++) {
            sfr[i] = (n0*magnitude[i] / base_mtf[i])/(65536*2);
        }
    } else {
        for (size_t i=0; i < sfr.size();  i++) {
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
    
    if (success > 0) {  // possibly contaminated edge
        quality = very_poor_quality;
    }
    
    return mtf50;
}

void Mtf_core::process_with_sliding_window(Mrectangle& rrect) {

    double winlen = 40; // desired length of ROI along edge direction

    const vector< vector<int> >& corner_map = rrect.corner_map;
    const vector<Point2d>& corners = rrect.corners;
    
    vector<Mtf_profile_sample> local_samples;

    vector< pair<double, int> > dims;
    for (int k=0; k < 4; k++) {
        dims.push_back(make_pair(norm(corners[corner_map[k][0]] - corners[corner_map[k][1]]), k));
    }
    sort(dims.begin(), dims.end());
    
    int v1 = dims[3].second;
    int v2 = dims[2].second;
    
    vector<Point2d> b(2);
    vector<Point2d> d(2);
    
    b[0] = corners[corner_map[v1][0]];
    d[0] = corners[corner_map[v1][1]] - corners[corner_map[v1][0]];
    
    b[1] = corners[corner_map[v2][0]];
    d[1] = corners[corner_map[v2][1]] - corners[corner_map[v2][0]];
    
    if (dims[2].first < winlen) {
        logger.debug("Rectangle not really long enough for sliding mode. Skipping.\n");
        return;
    }
    
    // first, refine edge orientation using (most) of the edge
    for (int side=0; side < 1; side++) {
        Point2d dir(d[side]);
        double edge_len = norm(dir);
        Point2d start(b[side]);
        dir *= 1.0/edge_len;
        
        Point2d n(-d[side].y, d[side].x);
        n *= 1.0/norm(n);
        double cross = dir.x*n.y - dir.y*n.x;
        if (cross < 0) {
            n = -n;
        }
        
        Point2d end(b[side] + edge_len*dir);
        
        Point2d tl(start + 16*n);
        Point2d br(end - 16*n);
        if (tl.x > br.x) {
            std::swap(tl.x, br.x);
        }
        if (tl.y > br.y) {
            std::swap(tl.y, br.y);
        }
        
        Point2d tl2(start - 16*n);
        Point2d br2(end + 16*n);
        if (tl2.x > br2.x) {
            std::swap(tl2.x, br2.x);
        }
        if (tl2.y > br2.y) {
            std::swap(tl2.y, br2.y);
        }
        tl.x = min(tl.x, tl2.x);
        tl.y = min(tl.y, tl2.y);
        br.x = max(br.x, br2.x);
        br.y = max(br.y, br2.y);
        
        Edge_record edge_record;
        
        // clamp ROi to image
        tl.x = max(1.0, tl.x);
        tl.y = max(1.0, tl.y);
        br.x = min(img.cols-1.0, br.x);
        br.y = min(img.rows-1.0, br.y);
        
        for (double y=tl.y; y < br.y; y += 1.0) {
            for (double x=tl.x; x < br.x; x += 1.0) {
            
                Point2d p(x,y);
                Point2d gd = p - b[side];
                double dot = gd.x*n.x + gd.y*n.y;
                double pdot = gd.x*dir.x + gd.y*dir.y;
                
                if (fabs(dot) < 12 && pdot > 5 && pdot < (edge_len - 5)) { // TODO: making the window narrow here helps a bit ...
                    int iy = lrint(y); 
                    int ix = lrint(x);
                    edge_record.add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                }
            }
        }
        
        edge_record.reduce(); // we can now move the ROI if we have to ...
        
        Point2d nd(sin(edge_record.angle), -cos(edge_record.angle));
        double dot = nd.x * dir.x + nd.y * dir.y;
        if (dot < 0) {
            nd = -nd;
        }
        
        // use angle to set d[], keep orientation
        d[side] = nd * edge_len;
        // TODO: should we update b[] as well?
    }
    
    // now process all the windows
    for (int side=0; side < 2; side++) {
    
        const double buffer = 5; // pixels to ignore near edge of block?
        double steplen = 4;
        Point2d dir(d[side]);
        double edge_len = norm(dir);
        
        int steps = floor((edge_len - winlen)/steplen) + 1;
        
        if (samples_per_edge != 0) {
            steps = max(2, min(steps, samples_per_edge));
            steplen = (edge_len - winlen) / double(steps-1);
        }
        
        Point2d start(b[side]);
        dir *= 1.0/edge_len;
        
        
        for (int step=0; step < steps; step++) {
            Point2d n(-d[side].y, d[side].x);
            n *= 1.0/norm(n);
            double cross = dir.x*n.y - dir.y*n.x;
            if (cross < 0) {
                n = -n;
            }
            
            Point2d end(b[side] + (step*steplen + winlen)*dir);
            
            Point2d tl(start + 16*n);
            Point2d br(end - 16*n);
            if (tl.x > br.x) {
                std::swap(tl.x, br.x);
            }
            if (tl.y > br.y) {
                std::swap(tl.y, br.y);
            }
            
            Point2d tl2(start - 16*n);
            Point2d br2(end + 16*n);
            if (tl2.x > br2.x) {
                std::swap(tl2.x, br2.x);
            }
            if (tl2.y > br2.y) {
                std::swap(tl2.y, br2.y);
            }
            tl.x = min(tl.x, tl2.x);
            tl.y = min(tl.y, tl2.y);
            br.x = max(br.x, br2.x);
            br.y = max(br.y, br2.y);
            
            tl.x = max(1.0, tl.x);
            tl.y = max(1.0, tl.y);
            br.x = min(img.cols-1.0, br.x);
            br.y = min(img.rows-1.0, br.y);
            
            map<int, scanline> scanset;
            Edge_record edge_record;
            
            double min_p = 1e50;
            double max_p = -1e50;
            
            for (double y=tl.y; y < br.y; y += 1.0) {
                for (double x=tl.x; x < br.x; x += 1.0) {
                
                    Point2d p(x,y);
                    Point2d gd = p - b[side];
                    double dot = gd.x*n.x + gd.y*n.y;
                    double pdot = gd.x*dir.x + gd.y*dir.y;
                    Point2d ld = p - start;
                    double ldot = (ld.x*dir.x + ld.y*dir.y) / winlen;
                    
                    if (pdot > buffer && pdot < (edge_len - buffer) && ldot > 0 && ldot < 1) {
                    
                        int iy = lrint(y); 
                        int ix = lrint(x);
                        if (fabs(dot) < 14) {
                            edge_record.add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                        }
                        
                        min_p = min(min_p, pdot);
                        max_p = max(max_p, pdot);
                    }
                }
            }
            
            edge_record.reduce(); // we can now move the ROI if we have to ... (usually we iterate a bit here...)
            
            #if 1
            Point2d newcent(edge_record.centroid.x, edge_record.centroid.y);
            edge_record.clear();
            
            for (double y=tl.y; y < br.y; y += 1.0) {
                for (double x=tl.x; x < br.x; x += 1.0) {
                
                    Point2d p(x,y);
                    Point2d wd = p - newcent;
                    double dot = wd.x*n.x + wd.y*n.y;
                    
                    Point2d ld = p - start;
                    double ldot = (ld.x*dir.x + ld.y*dir.y) / winlen;
                    
                    Point2d gd = p - b[side];
                    double pdot = gd.x*dir.x + gd.y*dir.y;
                    
                    if (pdot > buffer && pdot < (edge_len - buffer) && ldot > 0 && ldot < 1) {
                    
                        int iy = lrint(y); 
                        int ix = lrint(x);
                        
                        if (fabs(dot) < 12) {
                            edge_record.add_point(x, y, fabs(g.grad_x(ix,iy)), fabs(g.grad_y(ix,iy)));
                        }
                        
                        if (fabs(dot) < (max_dot + 1)) {
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
            
            edge_record.reduce();
            
            #endif
            
            cv::Vec3b& color = od_img.at<cv::Vec3b>(lrint(edge_record.centroid.y), lrint(edge_record.centroid.x));
            color[0] = 255;
            color[1] = 255;
            color[2] = 0;
            
            double quality = 0;
            vector <double> sfr(mtf_width, 0);
            vector <double> esf(FFT_SIZE/2, 0);
            vector <Point2d> ridge(1);
            double mtf50 = compute_mtf(edge_record.centroid, scanset, edge_record, quality, sfr, esf, ridge);
            
            if (mtf50 < 1.0 && quality > very_poor_quality) {
                local_samples.push_back(Mtf_profile_sample(edge_record.centroid, mtf50, edge_record.angle, quality));
            }
            
            start = b[side] + (step+1)*steplen*dir;
        }
    }
    
    if (local_samples.size() > 0) {
        std::lock_guard<std::mutex> lock(global_mutex);
        samples.insert(samples.end(), local_samples.begin(), local_samples.end());
    }
}


void Mtf_core::process_image_as_roi(void) { 
    map<int, scanline> scanset;
    Edge_record er;
    for (int row=0; row < img.rows; row++) {
        for (int col=0; col < img.cols; col++) {
            er.add_point(col, row, fabs(g.grad_x(col, row)), fabs(g.grad_y(col, row)));
            
            
            if (scanset.find(row) == scanset.end()) {
                scanset[row] = scanline(col,col);
            }
            scanset[row].start = std::min(col, scanset[row].start);
            scanset[row].end   = std::max(col, scanset[row].end);
        }
    }
    
    // ROI is the whole image. 
    er.reduce();
    Point2d cent(er.centroid);
    logger.debug("ER reduce grad estimate: %lf\n", er.angle / M_PI * 180);
    
    
    // scan the ROI to identify outliers, i.e., pixels far from the
    // edge with large gradient values (indicative of contamination of the ROI)
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
                    
                    if (g.grad_magnitude(x, y) >= upper) {
                        skiplist.push_back(y*img.cols + x);
                    } 
                } 
            }
        }
    }
    sort(skiplist.begin(), skiplist.end());
    
    // iterate the edge record to both centre the edge, as well as refine the edge orientation
    cent = er.centroid;
    Point2d normal(cos(er.angle), sin(er.angle)); 
    scanset.clear();
    er.clear();

    for (int row=0; row < img.rows; row++) {
        for (int col=0; col < img.cols; col++) {
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
    er.reduce();
    logger.debug("updated: ER reduce grad estimate: %lf, centroid (%lf,%lf) -> (%lf, %lf)\n", 
        er.angle/M_PI*180, cent.x, cent.y, er.centroid.x, er.centroid.y
    );
    
    double quality;
    
    vector <double> sfr(mtf_width, 0);
    vector <double> esf(FFT_SIZE/2, 0);
    vector <Point2d> ridge(1);
    double mtf50 = compute_mtf(er.centroid, scanset, er, quality, sfr, esf, ridge, true);
    
    // add a block with the correct properties ....
    if (mtf50 <= 1.2) { 
        Mrectangle rect;
        rect.centroids[0] = er.centroid;
        Block block(rect);
        block.centroid = er.centroid; // just in case
        block.set_mtf50_value(0, mtf50, quality);
        block.set_normal(0, Point2d(cos(er.angle), sin(er.angle)));
        
        block.set_sfr(0, sfr);
        block.set_esf(0, esf);
        block.set_ridge(0, ridge);
        
        for (int k=1; k < 4; k++) {
            block.set_mtf50_value(k, 1.0, 0.0);
            block.set_sfr(k, vector<double>(NYQUIST_FREQ*2, 0));
            block.set_esf(k, vector<double>(FFT_SIZE/2, 0));
        }
        
        shared_blocks_map[1] = block;
        detected_blocks.push_back(block);
    }
}

void Mtf_core::extract_ridge(map< int, pair<double, double> >& edge_residual, 
    const Point2d& cent, const Point2d& mean_grad, const Point2d& edge_direction,
    vector<Point2d>& ridge) {
    
    // we could just drop the first two and last two samples on principle
    edge_residual.erase(edge_residual.begin());
    edge_residual.erase(edge_residual.begin());
    auto it1 = edge_residual.end();
    it1--;
    edge_residual.erase(it1);
    it1 = edge_residual.end();
    it1--;
    edge_residual.erase(it1);
    it1 = edge_residual.end();
    it1--;
    double edgelen = std::min(abs(edge_residual.begin()->first), abs(it1->first));
    
    // apply some smoothing to the initial ridge to squash any severe outliers
    vector<Ordered_point> op_ridge(edge_residual.size());
    size_t idx = 0;
    for (auto& ee: edge_residual) {
        op_ridge[idx++] = Ordered_point(ee.first, ee.second.second);
    }
    
    sort(op_ridge.begin(), op_ridge.end());
    
    // smooth the first-order differences, rather than the actual signal
    vector<double> delta_perp(op_ridge.size()-1);
    double prev = op_ridge[0].second;
    for (size_t i=0; i < op_ridge.size() - 1; i++) {
        delta_perp[i] = op_ridge[i].second - prev;
        prev = op_ridge[i].second;
    }
    
    int hw = 3;
    vector<double> smoothed_delta(delta_perp.size(), 0);
    prev = op_ridge[0].second;
    for (int i=0; i < (int)delta_perp.size(); i++) {
        double wsum = 0;
        for (int j=-hw; j <= hw; j++) {
            if (i+j >= 0 && i+j < (int)delta_perp.size()) {
                wsum += 1;
                smoothed_delta[i] += delta_perp[i+j];
            }
            
        }
        smoothed_delta[i] /= wsum;
        double recon = prev + smoothed_delta[i];
        prev = recon;
        op_ridge[i].second = recon;
    }
    
    // now estimate the ridge using the gradient image as weight
    ridge = vector<Point2d>();
    for (auto& ee: op_ridge) {
        double pd = ee.second;
        double cd = ee.first;
        
        double wsum = 0;
        double xsum = 0;
        double ysum = 0;
        Point2d recon_cent = cent + mean_grad*pd + edge_direction*cd;
        int lcx = recon_cent.x;
        int lcy = recon_cent.y;
        for (int ly = lcy-15; ly <= lcy+15; ly++) {
            for (int lx = lcx-15; lx <= lcx+15; lx++) {
            
                if (lx < 0 || ly < 0 || lx > g.width()-1 || ly > g.height()-1) continue;
                
                // TODO: we should taper the weights near the edge of the ROI
                // in the along-edge direction
                //Point2d d(lx - cent.x, ly - cent.y);
                
                Point2d rd(lx - recon_cent.x, ly - recon_cent.y);
                double pardot = rd.ddot(edge_direction);
                double perpdot = rd.ddot(mean_grad);
                
                double end_weight = 1.0;
                
                if (fabs(perpdot) >  0.15) {
                    end_weight *= 1.0/(1 + fabs(perpdot));
                }
                
                if (fabs(pardot) > 0.5) {
                    end_weight *= 2.0/(1 + fabs(pardot));
                }
                
                // downrate the edges of the ROI
                if (fabs(Point2d(lx - cent.x, ly - cent.y).ddot(edge_direction)) > edgelen+1) {
                    end_weight *= 0.5;
                }
                
                // limit the window in the along-edge direction, but keep it rather wide initially to smooth out the estimate
                if (fabs(pardot) > 5) {  
                    continue; // equivalent to setting end_weight to zero
                }
                
                                        
                double w = g.grad_magnitude(lx, ly);
                w *= w;
                xsum += lx * w * end_weight;
                ysum += ly * w * end_weight;
                wsum += w * end_weight;
            }
        }
        xsum /= wsum;
        ysum /= wsum;
        
        // keep the refined values
        ridge.push_back(Point2d(xsum, ysum));
    }
    
    // iterate a few times to converge with tighter bounds
    vector<Point2d> ridge_copy;
    for (int rep=0; rep < 4; rep++) {
        ridge_copy = ridge;
        ridge.clear();
        for (auto& p: ridge_copy) {
            double wsum = 0;
            double xsum = 0;
            double ysum = 0;
            for (int ly = int(p.y-8); ly <= int(p.y+8); ly++) {
                for (int lx = int(p.x-8); lx <= int(p.x+8); lx++) {
                
                    if (lx < 0 || ly < 0 || lx > g.width()-1 || ly > g.height()-1) continue;
                    
                    Point2d rd(lx - p.x, ly - p.y);
                    double pardot = rd.ddot(edge_direction);
                    double perpdot = rd.ddot(mean_grad);
                    
                    double end_weight = 1.0;
                    
                    if (fabs(perpdot) >  0.15) {
                        end_weight *= 1.0/(1 + fabs(perpdot));
                    }
                    
                    if (fabs(pardot) > 0.5) {
                        end_weight *= 2.0/(1 + fabs(pardot));
                    }
                    
                    // downrate the edges of the ROI
                    if (fabs(Point2d(lx - cent.x, ly - cent.y).ddot(edge_direction)) > edgelen+1) {
                        end_weight *= 0.5;
                    }
                    
                    // limit the window in the along-edge direction, but keep the window tight to reduce serial correlation
                    if (fabs(pardot) > 1.0) {  
                        continue; // equivalent to setting end_weight to zero
                    }
                    
                                            
                    double w = g.grad_magnitude(lx, ly);
                    w *= w;
                    xsum += lx * w * end_weight;
                    ysum += ly * w * end_weight;
                    wsum += w * end_weight;
                }
            }
            xsum /= wsum;
            ysum /= wsum;
            
            // keep the refined values
            ridge.push_back(Point2d(xsum, ysum));
        }
    }
    
    // construct a heavily smoothed copy of the rige 
    // detect potential outliers by their distance from this smoothed ridge
    idx = 0;
    for (auto& p: ridge) {
        double perp = (p - cent).ddot(mean_grad);
        double par = (p - cent).ddot(edge_direction);
        op_ridge[idx++] = Ordered_point(par, perp);
    }
    
    prev = op_ridge[0].second;
    for (size_t i=0; i < op_ridge.size() - 1; i++) {
        delta_perp[i] = op_ridge[i].second - prev;
        prev = op_ridge[i].second;
    }
    
    vector<double> smoothed_recon(op_ridge.size());
    vector<double> residuals(smoothed_delta.size());
    hw = 9;
    prev = op_ridge[0].second;
    smoothed_recon[0] = op_ridge[0].second;
    for (int i=0; i < (int)delta_perp.size(); i++) {
        double wsum = 0;
        smoothed_delta[i] = 0;
        for (int j=-hw; j <= hw; j++) {
            if (i+j >= 0 && i+j < (int)delta_perp.size()) {
                wsum += 1;
                smoothed_delta[i] += delta_perp[i+j];
            }
            
        }
        smoothed_delta[i] /= wsum;
        double recon = prev + smoothed_delta[i];
        prev = recon;
        smoothed_recon[i+1] = recon;
        residuals[i] = fabs(recon - op_ridge[i].second);
    }
    
    sort(residuals.begin(), residuals.end());
    double thresh = residuals[0.9*residuals.size()]; // we can use a different test, perhaps look for a step in the residuals
    
    // convert back to ridge, but skip points that look like outliers
    ridge.clear();
    for (size_t i=0; i < op_ridge.size(); i++) {
        if (fabs(smoothed_recon[i] - op_ridge[i].second) < thresh) {
            ridge.push_back(cent + op_ridge[i].first * edge_direction + op_ridge[i].second * mean_grad);
        }
    }
}


inline void Mtf_core::update_gradient_peak(map< int, pair<double, double> >& edge_residual, 
    const double& dist_along_edge, const double& perp_dist, const double& grad_mag) {
    
    // discretize dist_along_edge, then bin it
    // then keep the pair <max_grad, dot> for each bin ...
    int dae_bin = lrint(dist_along_edge); // discretize at full-pixel intervals
    auto it = edge_residual.find(dae_bin);
    if (it == edge_residual.end()) {
        edge_residual.insert(make_pair(dae_bin, make_pair(grad_mag, perp_dist)));
    } else {
        if (grad_mag > it->second.first) {
            it->second = make_pair(grad_mag, perp_dist);
        }
    }
}

void Mtf_core::sample_at_angle(double ea, vector<Ordered_point>& local_ordered, 
    const map<int, scanline>& scanset, const Point2d& cent,
    double& edge_length, vector<Point2d>& ridge) {

    double max_along_edge = -1e50;
    double min_along_edge = 1e50;
    
    if (snap_to) {
        
        double max_dot_angle = snap_to_angle;
        double max_dot = 0;
        
        double angles[4] = {snap_to_angle, -snap_to_angle, M_PI/2 - snap_to_angle, snap_to_angle - M_PI/2};
        
        for (int k=0; k < 4; k++) {
        
            double sa = angles[k];
            
            double dot = cos(ea)*cos(sa) + sin(ea)*sin(sa);
            if (dot > max_dot) {
                max_dot = dot;
                max_dot_angle = sa;
            }
        }
        
        ea = max_dot_angle;
    }

    Point2d mean_grad(cos(ea), sin(ea));
    Point2d edge_direction(-sin(ea), cos(ea));
    
    map< int, pair<double, double> > edge_residual;
    
    if (!undistort) {
        if (bayer == Bayer::NONE) {
            for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                int y = it->first;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;
                    
                    Point2d d((x) - cent.x, (y) - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                        
                        update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(x, y));
                    }
                }
            }
            
        } else {
            for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                int y = it->first;
                int rowcode = (y & 1) << 1;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    int code = rowcode | (x & 1);
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;
                    
                    Point2d d((x) - cent.x, (y) - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    
                    bool inside_roi = fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length;
                    if (!inside_roi) continue;
                    
                    // we can update the ridge on all Bayer subsets, since the gradient was computed
                    // on the interpolated/demosaiced image
                    update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(x, y));

                    // skip the appropriate sites if we are operating only on a subset
                    if (bayer == Bayer::RED && code != 0) {
                        continue;
                    } 
                    if (bayer == Bayer::BLUE && code != 3) {
                        continue;
                    } 
                    if (bayer == Bayer::GREEN && (code == 0 || code == 3)) {
                        continue;
                    }

                    local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) ));
                    max_along_edge = max(max_along_edge, dist_along_edge);
                    min_along_edge = min(min_along_edge, dist_along_edge);
                }
            }
        }
    } else {
        // we have an undistortion model (e.g., equiangular mapping)
        
        // first construct a new scanset in distorted image coordinates
        map<int, scanline> m_scanset;
        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; ++x) {
                cv::Point2i tp = undistort->transform_pixel(x, y);
                auto fit = m_scanset.find(tp.y);
                if (fit == m_scanset.end()) {
                    m_scanset.insert(make_pair(tp.y, scanline(tp.x, tp.x)));
                } else {
                    fit->second.update(tp.x);
                }
            }
        }
        
        if (bayer == Bayer::NONE) {
            // now visit each pixel in the distorted image, but use undistorted coordinates to aid ESF construction
            for (map<int, scanline>::const_iterator it=m_scanset.begin(); it != m_scanset.end(); ++it) {
                int y = it->first;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;
                    
                    // transform warped pixel to idealized rectilinear
                    cv::Point2d tp = undistort->inverse_transform_point(x, y);
                    
                    bool inside_image = lrint(tp.x) >= 0 && lrint(tp.x) < g.width() && lrint(tp.y) >= 0 && lrint(tp.y) < g.height();
                    if (!inside_image) continue;
                    
                    Point2d d(tp.x - cent.x, tp.y - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    if (!undistort->rectilinear_equivalent()) {
                        Point2d base_u = dist_along_edge*edge_direction + cent; // on the linearized edge
                        Point2d base_d = undistort->transform_point(base_u.x, base_u.y);
                        Point2d dv = Point2d(x, y) - base_d;
                        dot = norm(dv) * (dv.ddot(mean_grad) < 0 ? -1.0 : 1.0);
                    }
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) )); // TODO: hack --- we are abusing the bayer image?
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                        
                        // (x, y) are equi-angular, but we need the rectilinear
                        // coordinates to go with the gradient image (which is rectilinear)
                        // TODO: should we use rounding, or floor here ?
                        //update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(lrint(tp.x), lrint(tp.y)));
                    }
                }
            }
        } else {
            for (map<int, scanline>::const_iterator it=m_scanset.begin(); it != m_scanset.end(); ++it) {
                int y = it->first;
                int rowcode = (y & 1) << 1;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    int code = rowcode | (x & 1);
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;
                    
                    // transform warped pixel to idealized rectilinear
                    cv::Point2d tp = undistort->inverse_transform_point(x, y);
                    
                    Point2d d(tp.x - cent.x, tp.y - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    
                    bool inside_roi = fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length && 
                        lrint(tp.x) >= 0 && lrint(tp.x) < g.width() && lrint(tp.y) >= 0 && lrint(tp.y) < g.height();
                    if (!inside_roi) continue;
                    
                    //update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(lrint(tp.x), lrint(tp.y)));

                    // skip the appropriate sites if we are operating only on a subset
                    if (bayer == Bayer::RED && code != 0) {
                        continue;
                    } 
                    if (bayer == Bayer::BLUE && code != 3) {
                        continue;
                    } 
                    if (bayer == Bayer::GREEN && (code == 0 || code == 3)) {
                        continue;
                    }
                    
                    if (!undistort->rectilinear_equivalent()) {
                        Point2d base_u = dist_along_edge*edge_direction + cent; // on the linearized edge
                        Point2d base_d = undistort->transform_point(base_u.x, base_u.y);
                        Point2d dv = Point2d(x, y) - base_d;
                        dot = norm(dv) * (dv.ddot(mean_grad) < 0 ? -1.0 : 1.0);
                    }
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) ));
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                    }
                }
            }
        }
    }

    if (!undistort) extract_ridge(edge_residual, cent, mean_grad, edge_direction, ridge);
    edge_length = max_along_edge - min_along_edge;
}
