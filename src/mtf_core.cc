#include "include/mtf_core.h"

#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include "include/peak_detector.h"

#include "include/point_helpers.h"

// global lock to prevent race conditions on detected_blocks
tbb::mutex global_mutex;

void Mtf_core::search_borders(const Point& cent, int label) {
    
    Rectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    Block block(rrect);
    
    if (!valid) {
        return;
    }
    
    int current_block_idx=0;
    {
        tbb::mutex::scoped_lock lock(global_mutex);
        detected_blocks.push_back(block);
        current_block_idx = detected_blocks.size() - 1;
    }
    
    vector<Point>& centroids = rrect.centroids;
    for (size_t k=0; k < 4; k++) {
        Point mid_dir = average_dir(g, int(centroids[k].x), int(centroids[k].y));
        
        // now construct buffer around centroid, aligned with direction, of width max_dot
        Rectangle nr(rrect, k, max_dot);
        
        map<int, scanline> scanset;
        for (double y=nr.tl.y; y < nr.br.y; y += 1.0) {
            for (double x=nr.tl.x; x < nr.br.x; x += 1.0) {
                Point p(x,y);
                if (nr.is_inside(p)) {
                
                    int iy = lrint(y);
                    int ix = lrint(x);
                    if (iy > 0 && iy < img.rows && ix > 0 && ix < img.cols) {
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
        double mtf50 = compute_mtf(centroids[k], scanset);
        
        if (mtf50 <= 1.0) { // reject mtf values above 1.0, since these are impossible
            detected_blocks[current_block_idx].set_mtf50_value(k, mtf50);
        }
    }
    
}

bool Mtf_core::extract_rectangle(const Point& cent, int label, Rectangle& rect) {
    
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
    
    Rectangle rrect(main_thetas, thetas, points);
    rect = rrect;
    
    return rrect.valid;
}

double Mtf_core::compute_mtf(const Point& in_cent, const map<int, scanline>& scanset) {
    
    Point cent(in_cent);
    
    Point mean_grad(0,0);
    double wsum = 0;
    int count = 0;
    // compute direction onto which points will be projected
    Point ncent(0.0, 0.0);
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
        int y = it->first;
        for (int x=it->second.start; x <= it->second.end; x++) {
            double gm = g.grad_magnitude(x,y);
            mean_grad.x += g.grad_x(x,y) * gm;
            mean_grad.y += g.grad_y(x,y) * gm;
            ncent.x += x * gm;
            ncent.y += y * gm;
            wsum += gm;
            count++;
        }
    }
    if (fabs(wsum) < 1e-6) {
        return 0;
    }
    mean_grad.x /= wsum;
    mean_grad.y /= wsum;
    ncent.x /= wsum;
    ncent.y /= wsum;
    mean_grad = normalize(mean_grad);
    
    Point old_grad = mean_grad;
    mean_grad = Point(0.0,0.0);
    wsum = 0;
    count = 0;
    
    // the exact location of the edge is probably not critical (except for the influence on apodization?)
    // TODO: investigate this
    cent = Point((cent.x + ncent.x)/2.0, (cent.y + ncent.y)/2.0);
    
    // recompute direction onto which points will be projected
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
        int y = it->first;
        for (int x=it->second.start; x <= it->second.end; x++) {
            Point d(x - cent.x, y - cent.y);
            double dot = d.ddot(old_grad);
            if (fabs(dot) < max_dot) {
                double gm = g.grad_magnitude(x,y);
                gm *= exp(-dot*dot/16.0);
                mean_grad.x += g.grad_x(x,y) * gm;
                mean_grad.y += g.grad_y(x,y) * gm;
                wsum += gm;
                count++;
            }
        }
    }
    mean_grad.x /= wsum;
    mean_grad.y /= wsum;
    mean_grad = normalize(mean_grad);
    
    
    // orientation estimates remain weak. do something about it
    double angle = atan2(mean_grad.y, mean_grad.x);
    
    vector<Ordered_point> ordered;
    double min_sum = 1e50;
    double best_angle = 0;
    for (double ea=angle-2.0/180.0*M_PI; ea < angle + 2.0/180.0; ea += 0.1/180.0*M_PI) {
        
        mean_grad.x = cos(ea);
        mean_grad.y = sin(ea);
    
        vector<Ordered_point> local_ordered;
        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                Point d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                if (fabs(dot) < max_dot) {
                    local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                }
            }
        }
        // measure variance of adjacent points
        sort(local_ordered.begin(), local_ordered.end());
        double sum = 0;
        for (size_t i=1; i < local_ordered.size(); i++) {
            sum += SQR(fabs(local_ordered[i-1].second - local_ordered[i].second));
        }
        if (sum < min_sum) {
            min_sum = sum;
            ordered = local_ordered;
            best_angle = ea;
        }
    }
    
    mean_grad.x = cos(best_angle);
    mean_grad.y = sin(best_angle);
    
    if (ordered.size() < 10) {
        return 1000;
    }
    
    const int nsteps = FFT_SIZE+1;
    double x_span = (ordered.back().first - ordered[0].first);
    double step = x_span / nsteps;
    
    double* fft_in_buffer = (double*)fftw_malloc(sizeof(double)*2*(FFT_SIZE+2));
    for (size_t i=0; i < 2*(FFT_SIZE+2); i++) {
        fft_in_buffer[i] = 0.0;
    }
    
    loess_fit(ordered, fft_in_buffer, FFT_SIZE); // loess_fit computes the ESF derivative as part of the fitting procedure
    
    int nc = (FFT_SIZE)  + 1;
    fftw_complex* fft_out_buffer = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nc);
    
    fftw_execute_dft_r2c(plan_forward, fft_in_buffer, fft_out_buffer);
    
    double n0 = sqrt(SQR(fft_out_buffer[0][0]) + SQR(fft_out_buffer[0][1]));
    double prev_freq = 0;
    double prev_val  = n0;
    
    double mtf50 = 0;

    bool done = false;
    for (size_t i=0; i < FFT_SIZE/2 && !done; i++) {
        double mag = sqrt(SQR(fft_out_buffer[i][0]) + SQR(fft_out_buffer[i][1])) / n0;
        if (prev_val > 0.5 && mag <= 0.5) {
            // interpolate
            double m = -(mag - prev_val)/(step/x_span);
            mtf50 = -(0.5 - prev_val - m*prev_freq) / m;
            done = true;
        }
        prev_val = mag;
        prev_freq = i / double(FFT_SIZE);
    }
    
    mtf50 *= 2*max_dot; // twice the span
    
    fftw_free(fft_out_buffer);
    fftw_free(fft_in_buffer);        
    
    return mtf50;
}

