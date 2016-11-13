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
#include "include/ellipse_decoder.h"

#include <mutex>

// global lock to prevent race conditions on detected_blocks
static std::mutex global_mutex;

void Mtf_core::search_borders(const Point2d& cent, int label) {
    
    Mrectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    if (!valid && sliding) {
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
                    printf("Fiducial with code %d extracted at (%.2lf %.2lf)\n", e.code, e.centroid_x, e.centroid_y);
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
            }
        } 
        return;
    }
    
    Block block(rrect);

    if (block.get_area() <= 225) {
        return;
    }
    
    if (!rrect.corners_ok()) {
        printf("discarding broken square (early)\n");
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
                }
            }
        }
        
        reduce_success &= edge_record[k].reduce();
    }
    
    double max_shift = 0;
    if (reduce_success) {
        // re-calculate the ROI after we have refined the edge centroid above
        Mrectangle newrect(rrect, edge_record);
        if (!newrect.corners_ok()) {
            printf("discarding broken square (after updates)\n");
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
            printf("discarding broken square (after updates)\n");
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
            reduce_success &= edge_record[k].reduce(true);
        }
    }
    
    if (!reduce_success) {
        printf("reduce failed, probably not a rectangle/quadrangle\n");
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
        vector <double> sfr(NYQUIST_FREQ*2, 0);
        vector <double> esf(FFT_SIZE/2, 0);
        double mtf50 = compute_mtf(edge_record[k].centroid, scansets[k], edge_record[k], quality, sfr, esf);
        
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
            printf("block with centroid (%d, %d) failed interior check\n", ix, iy);
            return false;
        }
        if (cl(lrint(0.5*(ix+rrect.corners[0].x)), lrint(0.5*(iy+rrect.corners[0].y))) != label) {
            printf("block with centroid (%d, %d) failed interior (corner) check\n", ix, iy);
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
    vector<double>& sfr, vector<double>& esf) {
    
    quality = 1.0; // assume this is a good edge
    
    Point2d cent(in_cent);
    
    Point2d mean_grad(0,0);
   

    double angle = er.angle;
    mean_grad.x = cos(angle);
    mean_grad.y = sin(angle);

    vector<Ordered_point> ordered;
    double edge_length = 0;

    vector<double> fft_out_buffer(FFT_SIZE * 2, 0);
    
    sample_at_angle(angle, ordered, scanset, cent, edge_length);
    sort(ordered.begin(), ordered.end());
    
    if (ordered.size() < 10) {
        quality = 0; // this edge is not usable in any way
        return 0;
    }
    
    int success = bin_fit(ordered, fft_out_buffer.data(), FFT_SIZE, -max_dot, max_dot, esf); // bin_fit computes the ESF derivative as part of the fitting procedure
    if (success < 0) {
        quality = poor_quality;
        printf("failed edge\n");
        return 1.0;
    }
    afft.realfft(fft_out_buffer.data());

    double quad = angle_reduce(angle);
    
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
                int filter_order = min(5, (idx-5)/stride);
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
        printf("Rectangle not really long enough for sliding mode. Skipping.\n");
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
            vector <double> sfr(NYQUIST_FREQ*2, 0);
            vector <double> esf(FFT_SIZE/2, 0);
            double mtf50 = compute_mtf(edge_record.centroid, scanset, edge_record, quality, sfr, esf);
            
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
