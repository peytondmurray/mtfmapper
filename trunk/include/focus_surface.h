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
#ifndef FOCUS_SURFACE_H
#define FOCUS_SURFACE_H

#include "ratpoly_fit.h"

class Focus_surface  {
  public:
  
    Focus_surface(vector< pair<Point, double> >& data, int order_n, int order_m, double cscale,
        const vector<cv::Point3d>& distance_scale) 
    : data(data), order_n(order_n), order_m(order_m), maxy(-1e50), maxx(-1e50), cscale(cscale) {
    
        double miny = 1e50;
        for (size_t i=0; i < data.size(); i++) {
            maxx = std::max(fabs(data[i].first.x), maxx);
            maxy = std::max(fabs(data[i].first.y), maxy);
            miny = std::min(fabs(data[i].first.y), miny);
        }
        
        
        map<double, double> peaks;
        
        FILE* fpeaks = fopen("fpeaks.txt", "wt");
        
        FILE* pfile = fopen("profiles.txt", "wt");
        FILE* reconf = fopen("recon.txt", "wt");
        int pindex = 0;
        
        vector<Sample> peak_pts;
        // |15 to 110| in steps of 2.5, width=5 ??
        for (int s=-1; s <= 1; s+=2) {
            for (double d=15; d <= maxy*cscale; d += 1.0) {
                double midy = s*d;
                vector<Sample> pts_row;
                vector<double> real_y;
                for (size_t i=0; i < data.size(); i++) {
                    double dy = midy - data[i].first.y*cscale;
                    if (fabs(dy) < 15 && fabs(data[i].first.y*cscale) > 10 ) { // at least 10 mm from centre of chart
                        double yw = exp(-dy*dy/(2*10*10)); // sdev of 3 mm in y direction
                        pts_row.push_back( Sample(data[i].first.x * cscale, data[i].second, yw, 0.4 + 0.6*std::max(0.01, data[i].second)) );
                        real_y.push_back(data[i].first.y*cscale);
                    } 
                }
                
                if (pts_row.size() < 50) continue;
                
                VectorXd sol;
                double lpeak;
                double peak_mtf;
                
                
                // TODO: get rid of the rep loop (and overheads)
                for (int rep=0; rep < 1; rep++) {
                    double mean_x = 0;
                    double wsum = 0;
                    fprintf(pfile, "#rep=%d\n", rep);
                    
                    vector<Sample> sample;
                    if (rep == 0) {
                        sample=pts_row;
                        for (size_t j=0; j < sample.size(); j++) {
                            mean_x += pts_row[j].weight * real_y[j];
                            wsum += pts_row[j].weight;
                        }
                    } else {
                        
                        while (sample.size() < pts_row.size()) {
                            int rindex = double(rand() * pts_row.size()) / double(RAND_MAX);
                            sample.push_back(pts_row[rindex]);
                            mean_x += pts_row[rindex].weight * real_y[rindex];
                            wsum += pts_row[rindex].weight;
                        }
                        
                    }
                    mean_x /= wsum;
                    
                    Ratpoly_fit cf(sample, order_n, order_m);
                    sol = rpfit(cf, false, true);
                    while (cf.order_n > 1 && cf.order_m > 0 && cf.has_poles(sol)) {
                        cf.order_m--;
                        sol = rpfit(cf, false, true);
                        printf("reducing order_m to %d\n", cf.order_m);
                        if (cf.has_poles(sol)) {
                            cf.order_n--;
                            printf("reducing order_n to %d\n", cf.order_n);
                            sol = rpfit(cf, false, true);
                        }
                    }
                    if (cf.has_poles(sol)) { 
                        // no solution without poles, give up, skip this sample?
                        printf("Warning: no viable RP fit. Skipping curve centred at y=%lf\n", mean_x);
                        continue;
                    }
                    printf("final: (%d, %d)\n", cf.order_n, cf.order_m);
                    double err = cf.evaluate(sol);
                    fprintf(stderr, "%lf\n", err);
                    lpeak = cf.peak(sol);
                    peak_mtf = cf.rpeval(sol, lpeak);
                    
                    #if 0
                    if (true) {
                        fprintf(pfile, "#index=%d\n", pindex);
                        for (size_t j=0; j < sample.size(); j++) {
                            fprintf(pfile, "%lf %lf %lf\n", sample[j].x, sample[j].y, sample[j].weight);
                        }
                        fprintf(pfile, "\n\n");
                        printf("index %d: lpeak=%lf (%lf mm)\n", pindex, lpeak/cscale, lpeak);
                        printf("error=%lf\n", err);
                        
                        fprintf(reconf, "#index=%d\n", pindex);
                        for (double x=-maxx*cscale; x < maxx*cscale; x += 0.1) {
                            double y = cf.rpeval(sol, x);
                            fprintf(reconf, "%lf %lf\n", x, y);
                        }
                        fprintf(reconf, "\n\n");
                        
                        
                        pindex++;
                    }
                    #endif
                    
                    
                    peaks[midy] = lpeak;
                    rpolys[midy] = sol;
                    
                    // TODO: we could weight each peak by its uncertainty ...
                    
                    
                    double pw = 0.1;
                    double xw = fabs(midy)/(0.7*maxy*cscale);
                    if (xw < 1) {
                        pw = 0.1 + 0.9*(1 - xw*xw*xw)*(1 - xw*xw*xw)*(1 - xw*xw*xw);
                    }
                    peak_pts.push_back( Sample(mean_x/cscale, lpeak/cscale, pw, 1.0/(1e-4 + err)) );
                    
                    fprintf(fpeaks, "%lf %lf %lf\n", peak_pts.back().y, peak_pts.back().x, 1.0);
                }
            }
        }
        fclose(fpeaks);
        fclose(pfile);
        fclose(reconf);
        
        Ratpoly_fit cf(peak_pts, 2,1);
        cf.base_value = 1;
        cf.pscale = 0;
        VectorXd sol = rpfit(cf, true, true);
        std::cout << "sol = " << sol.transpose() << std::endl;
        printf("evaluations=%d\n", (int)cf.evaluations());
        
        while (cf.order_n > 1 && cf.order_m > 0 && cf.has_poles(sol)) {
            cf.order_m--;
            printf("reducing order_m to %d\n", cf.order_m);
            sol = rpfit(cf, true, true);
            if (cf.has_poles(sol)) {
                cf.order_n--;
                printf("reducing order_n to %d\n", cf.order_n);
                sol = rpfit(cf, true, true);
            }
        }
        if (cf.has_poles(sol)) { 
            // no solution without poles, give up, skip this sample?
            printf("Warning: no viable RP fit to fpeaks data\n");
        }
        printf("final: (%d, %d)\n", cf.order_n, cf.order_m);
        
        
        FILE* ridgef = fopen("ridge.txt", "wt");
        
        double cmax = -1e50;
        double cmin = 1e50;
        for (double y=-maxy; y < maxy; y += 1) {
            double x = cf.rpeval(sol, y*cf.xsf)/cf.ysf;
            if (fabs(y) < maxy/4.0) { // to avoid poles near the edges
                cmax = std::max(cmax, x);
                cmin = std::min(cmin, x);
            }
            fprintf(ridgef, "%lf %lf %lf\n", x, y, 1.0);
        }
        fclose(ridgef);
        
        // TODO: if the extremum is on one of the edge, then the chart is skew
        
        printf("ridge extrema: %lf %lf\n", cmin, cmax);
        double x_inter = cf.rpeval(sol, 0)/cf.ysf;
        printf("intercept %lg %lg\n", x_inter * cscale, x_inter);
        
        if (distance_scale.size() > 0) {
            printf("distance scale is:\n");
            for (size_t i=0; i < distance_scale.size(); i++) {
                printf("pix %lf -> dist %lf\n", distance_scale[i].x, distance_scale[i].y);
            }
            
            // find the two centre-most scale markers, average their distance to estimate chart angle
            int middle = 0;
            for (int i=1; i < distance_scale.size(); i++) {
                if (fabs(distance_scale[i].x) < fabs(distance_scale[middle].x)) {
                    middle = i;
                }
            }
            double foreshortening = 0.5*(fabs(distance_scale[middle-1].x) + fabs(distance_scale[middle+1].x));
            foreshortening *= cscale/fabs(distance_scale[middle-1].y);
            
            // x_inter is in pixels, relative to centre of chart
            int scale_lower=0;
            while (scale_lower < distance_scale.size() - 1 &&
                   distance_scale[scale_lower].x < x_inter) {
                scale_lower++;
            }
            printf("scale limits: %d, %d : %lf, %lf\n", 
                scale_lower, scale_lower+1, 
                distance_scale[scale_lower].x, distance_scale[scale_lower+1].x
            );
            
            const cv::Point3d& p0 = distance_scale[scale_lower];
            const cv::Point3d& p1 = distance_scale[scale_lower+1];
            
            double slope = (p1.y - p0.y) / (p1.x - p0.x);
            double offset = p0.y - slope*p0.x;
            double focus_plane_position = offset + x_inter * slope;
            printf("foreshortening=%lf\n", foreshortening);
            printf("focus_plane %lg\n", focus_plane_position * foreshortening);
        }
    }
    
    VectorXd rpfit(Ratpoly_fit& cf, bool scale=false, bool refine=false) {
        const vector<Sample>& pts_row = cf.get_data();
        
        if (scale) {
            double xsf=0;
            double ysf=0;
            for (size_t i=0; i < pts_row.size(); i++) {
                xsf = std::max(xsf, fabs(pts_row[i].x));
                ysf = std::max(ysf, fabs(pts_row[i].y));
            }
            cf.xsf = xsf = 1.0/xsf;
            cf.ysf = ysf = 1.0/ysf;
        }
        
        int tdim = cf.dimension();
        MatrixXd cov = MatrixXd::Zero(tdim, tdim);
        VectorXd b = VectorXd::Zero(tdim);
        VectorXd a = VectorXd::Zero(tdim);
        
        VectorXd sol;
        
        for (int iter=0; iter < 1; iter++) {
            cov.setZero();
            b.setZero();
            a.setZero();
            
            for (size_t i=0; i < pts_row.size(); i++) {
                const Sample& sp = pts_row[i];
                double w = sp.weight * sp.yweight;
                a[0] = 1*w;
                double prod = sp.x * cf.xsf; // top poly
                for (int j=1; j <= cf.order_n; j++) {
                    a[j] = prod*w;
                    prod *= sp.x * cf.xsf;
                }
                prod = sp.x*cf.xsf; // bottom poly
                for (int j=1; j <= cf.order_m; j++) {
                    a[j+cf.order_n] = prod*w*sp.y*cf.ysf;
                    prod *= sp.x*cf.xsf;
                }
                
                for (int col=0; col < tdim; col++) { 
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov(col, icol) += a[col]*a[icol];
                    }
                    b[col] += cf.base_value*a[col]*sp.y*cf.ysf*w; // build rhs of system : A'*b
                }
            }
            
            for (int col=cf.order_n+1; col < cov.cols(); col++) {
                cov.col(col) = -cov.col(col);
            }
            
            sol = cov.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        }
        
        // now perform non-linear optimization
        
        if (refine) {
            sol = cf.gauss_newton_armijo(sol);
        }
        
        return sol;
    }

    double evaluate(VectorXd& v) {
        return 0;
    }
    
    int dimension(void) {
        return (order_n+1 + order_m);
    }
    
    vector< pair<Point, double> >& data;
    int order_n;
    int order_m;
    double maxy;
    double maxx;
    double cscale;
    map<double, VectorXd> rpolys;
};

#endif
