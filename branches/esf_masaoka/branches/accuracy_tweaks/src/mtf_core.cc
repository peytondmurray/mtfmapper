#include "include/mtf_core.h"

#include "include/loess_fit.h"
#include "include/gaussfilter.h"
#include "include/peak_detector.h"

#include "include/point_helpers.h"

// global lock to prevent race conditions on detected_blocks
tbb::mutex global_mutex;

// coefficients of 7th order polynomials to correct mtf50 estimate bias
double mtf_correction_coeffs[5][9] = {
    {5,-0.0136553397922912,1.18082859421591,-1.11330042943354,3.65788287506631,-5.73335506231191,4.10368254171464,-0.421536084917037,-0.476675005068112},
    {10,-0.0165032895550072,1.29559924481722,-2.85364572749257,16.3563830930782,-54.4785309519634,104.302782221934,-104.578685812887,42.4647880143066},
    {20,-0.0102160182500555,1.07714702981703,-0.172221706906046,2.02150323776689,-24.1228546777887,102.451902516557,-165.112225506796,92.0016162889815},
    {30,-0.00751081552319845,0.944200769839098,2.33530303222745,-21.3027461001393,91.2709200999501,-201.337979929428,225.309979449181,-100.424743668452},
    {40,-0.0111058387325651,1.09433524520547,-0.082433986559942,-2.14050391948704,10.7239220103593,-19.0654918966613,13.5541442033553,-2.40762452892175},
};
                    

void Mtf_core::search_borders(const Point& cent, int label) {
    
    Mrectangle rrect;
    bool valid = extract_rectangle(cent, label, rrect);
    
    if (!valid) {
        return;
    }
    
    Block block(rrect);
    
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
        Mrectangle nr(rrect, k, max_dot);
        
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

void simple_binning(vector< Ordered_point  >& ordered, double* fft_in_buffer, const int fft_size, double max_dot);

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
    double best_angle = angle;
    for (double ea=angle-2.0/180.0*M_PI; ea < angle + 2.0/180.0; ea += 0.03/180.0*M_PI) {
        
        mean_grad.x = cos(ea);
        mean_grad.y = sin(ea);
        
        Point edge_direction(sin(ea), cos(ea));
    
        vector<Ordered_point> local_ordered;
        for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                Point d((x) - cent.x, (y) - cent.y);
                double dot = d.ddot(mean_grad); 
                double dist_along_edge = d.ddot(edge_direction);
                if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                    local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                }
            }
        }
        // measure variance of adjacent points
        sort(local_ordered.begin(), local_ordered.end());
        double sum = 0;
        for (size_t i=3; i < local_ordered.size(); i++) {
            sum += SQR(fabs(local_ordered[i-1].second - local_ordered[i].second));
            sum += SQR(fabs(local_ordered[i-3].second - local_ordered[i].second));
        }
        #if 1
        if (sum < min_sum) {
            min_sum = sum;
            ordered = local_ordered;
            best_angle = ea;
        }
        #else
        if (fabs(ea - angle) < 1e-5) {
            ordered = local_ordered;
            best_angle = angle;
        }
        #endif
    }
    
    printf("best angle = %lf\n", best_angle / M_PI * 180);
    
    mean_grad.x = cos(best_angle);
    mean_grad.y = sin(best_angle);
    
    if (ordered.size() < 10) {
        return 1000;
    }
    
    #if 0
    for (size_t i=0; i < ordered.size(); i++) {
        fprintf(stderr, "%lf %lf\n", ordered[i].first, ordered[i].second);
    }
    #endif
    
    double* fft_in_buffer = (double*)fftw_malloc(sizeof(double)*2*(FFT_SIZE+2));
    for (size_t i=0; i < 2*(FFT_SIZE+2); i++) {
        fft_in_buffer[i] = 0.0;
    }

    loess_fit(ordered, fft_in_buffer, FFT_SIZE, -max_dot, max_dot); // loess_fit computes the ESF derivative as part of the fitting procedure
    
    //apodize(fft_in_buffer, FFT_SIZE);
    //simple_binning(ordered, fft_in_buffer, FFT_SIZE, max_dot);
    
    fftw_complex* fft_out_buffer = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (FFT_SIZE+1));
    fftw_execute_dft_r2c(plan_forward, fft_in_buffer, fft_out_buffer);
    
    double n0 = sqrt(SQR(fft_out_buffer[0][0]) + SQR(fft_out_buffer[0][1]));
    double prev_freq = 0;
    double prev_val  = n0;
    
    double mtf50 = 0;
    
    #if 0
    for (size_t i=0; i < FFT_SIZE; i++) {
         fprintf(stderr, "%lf %lf\n", i/double(FFT_SIZE), fft_in_buffer[i]);
    }
    #endif
    
    #if 0
    for (size_t i=0; i < SAMPLES_PER_PIXEL*2; i++) {
        double mag = sqrt(SQR(fft_out_buffer[i][0]) + SQR(fft_out_buffer[i][1]));
        double delta = 0.25;
        double x = M_PI*delta*2*(i*32+1)/double(FFT_SIZE);
        double magc = mag * (sin(x)/x);
        fprintf(stderr, "%lf %lf %lf %lf\n", i*SAMPLES_PER_PIXEL/double(FFT_SIZE), mag/n0, magc/n0, sin(x)/x);
    }
    #endif

    bool done = false;
    for (size_t i=0; i < FFT_SIZE/2 && !done; i++) {
        double mag = sqrt(SQR(fft_out_buffer[i][0]) + SQR(fft_out_buffer[i][1])) / n0;
        
        //double delta = 0.25;
        //double x = M_PI*delta*2*(i*32+1)/double(FFT_SIZE);
        //mag = mag / (sin(x)/x);
        
        if (prev_val > 0.5 && mag <= 0.5) {
            // interpolate
            double m = -(mag - prev_val)*(FFT_SIZE);
            mtf50 = -(0.5 - prev_val - m*prev_freq) / m;
            done = true;
        }
        prev_val = mag;
        prev_freq = i / double(FFT_SIZE);
    }
    
    mtf50 *= SAMPLES_PER_PIXEL; 
    //mtf50 *= 4;
    
    #if 0
    // find closest angle
    double quad1 = fabs(fmod(best_angle, M_PI/2.0));
    if (quad1 > M_PI/4.0) {
        quad1 = M_PI/2.0 - quad1;
    }
    quad1 = quad1 / M_PI * 180;
    size_t angle_idx = 0;
    for (size_t i=1; i < 5; i++) {
        if (fabs(mtf_correction_coeffs[i][0] - quad1) < 
            fabs(mtf_correction_coeffs[angle_idx][0] - quad1)) {
        
            angle_idx = i;
        }    
    }
    /*
    printf("best angle = %lf, closest match is %lf, quad1 = %lf\n",
        best_angle / M_PI*180.0, 
        mtf_correction_coeffs[angle_idx][0], 
        quad1
    );
    */
    
    double s = mtf_correction_coeffs[angle_idx][1];
    double xp = mtf50;
    for (int i=2; i < 9; i++) {
        s += xp * mtf_correction_coeffs[angle_idx][i];
        xp *= mtf50;
    }
    //printf("initial = %lf, proposed = %lf, correction = %lf\n", mtf50, s, s - mtf50);
    mtf50 = s;
    #endif
    
    fftw_free(fft_out_buffer);
    fftw_free(fft_in_buffer);        
    
    return mtf50;
}

void Mtf_core::apodize(double* data, int dlen) {
    
    vector<double> ma_filtered(dlen, 0);
    
    // produce smoothed PSF
    double sum = 0;
    int count = 0;
    const int w = 72;
    for (int i=0; i < w/2; i++) {
        sum += data[i];
        count++;
    }
    int peak_idx = 0;
    for (int i=0; i < dlen; i++) {
        ma_filtered[i] = sum / double(count);
        if (i < dlen - w/2) {
            sum += data[i+w/2];
            if (count < w) {
                count++;
            } else if (i >= w/2) {
                sum -= data[i-w/2];
            }
        } else {
            sum -= data[i-w/2];
            count--;
        }
        if (ma_filtered[i] > ma_filtered[peak_idx]) {
            peak_idx = i;
        }
        
    }
    
    int lower = peak_idx;
    while (lower > 0 && ma_filtered[lower] > ma_filtered[peak_idx]*0.2) lower--;
    int upper = peak_idx;
    while (upper < dlen-1 && ma_filtered[upper] > ma_filtered[peak_idx]*0.2) upper++;
    
    int pw20 = upper - lower;
    
    int alower = lower - pw20 - 4;
    int aupper = upper + pw20 + 4;
    
    //printf("peak width = %d samples, al=%d, au=%d\n", pw20, alower, aupper);
    
    for (int i=0; i <= alower; i++) {
        data[i] = ma_filtered[i];
    }
    for (int i=aupper; i < dlen; i++) {
        data[i] = ma_filtered[i];
    }
    /*
    for (int i=0; i < dlen; i++) {
        fprintf(stderr, "%d %lf %lf\n", i, ma_filtered[i], data[i]);
    }
    */
}


void simple_binning(vector< Ordered_point  >& ordered, double* fft_in_buffer, const int fft_size, double max_dot) {
    double lower = -max_dot;
    double upper = -max_dot + 0.25;
    size_t lower_idx = 0;
    size_t idx = 0;
    for (int i=0; i < fft_size; i++) {
        fft_in_buffer[i] = 0;
    }
    for (int i=0; i < fft_size; i+=8) {
        double sample = 0;
        while (idx < ordered.size() && ordered[idx].first < upper) idx++;
        
        for (size_t j=lower_idx; j < idx; j++) {
            sample += ordered[idx].second;
        }
        if (idx > lower_idx) {
            sample /= double(idx - lower_idx);
        }
        /*
        for (size_t j=i; j < i+8; j++) {
            fft_in_buffer[j] = sample;
        }
        */
        fft_in_buffer[i/8 + fft_size/2 - 64] = sample;
        
        lower_idx = idx;
        lower = upper;
        upper += 0.25;   
    }
    for (int i=fft_size/2+62; i < fft_size; i++) {
        fft_in_buffer[i] = fft_in_buffer[i-1];
    }
    for (int i=2; i < fft_size; i++) {
        fft_in_buffer[i-2] = 0.5*(fft_in_buffer[i] - fft_in_buffer[i-2]); 
        double w = 0.54 + 0.46*cos(2*M_PI*(i - fft_size/2)/double(fft_size-1));  // Hamming window function
        fft_in_buffer[i-2] *= w;
    }
    for (int i=fft_size-8; i < fft_size; i++) {
        fft_in_buffer[i] = 0;
    }
}
