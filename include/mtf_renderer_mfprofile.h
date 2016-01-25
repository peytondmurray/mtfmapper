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
#ifndef MTF_RENDERER_MFPROFILE_H
#define MTF_RENDERER_MFPROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"

#include <stdlib.h>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#ifndef SQR 
#define SQR(x) ((x)*(x))
#endif

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

class Evaluable {
  public:
    Evaluable(void) : evaluation_count(0) {}
    virtual ~Evaluable(void) {}
    
    virtual double evaluate(VectorXd& v) = 0;  // we promise not to harm v
    virtual int dimension(void) = 0;
    
    unsigned long long evaluations(void) {
        return evaluation_count;
    }
                                
  protected:
    unsigned long long evaluation_count;
};

template < class T > inline void SWAP(T & a, T & b) {
    T t = a;
    a = b;
    b = t;
}

class Sample {
  public:
    Sample(double x, double y, double weight, double yweight=1) 
    :x(x), y(y), weight(weight), yweight(yweight) {}
    
    double x;
    double y;
    double weight;
    double yweight;
};

class Curve_fit : public Evaluable {
  public:
    typedef enum {LS, MAD, LOG} fit_type;
    
    Curve_fit(const vector<Sample>& data, int order_n, int order_m, fit_type fit=LS)
    : data(data), order_n(order_n), order_m(order_m), base_value(1.0), xsf(1), ysf(1), pscale(0.1), fit(fit) {}
    
    double evaluate(VectorXd& v) {
        double err = 0;
        for (size_t i=0; i < data.size(); i++) {
            double w = data[i].weight * data[i].yweight;
            double z = rpeval(v, data[i].x*xsf);
            double e = data[i].y*ysf - z;
            switch(fit) {
            case MAD: err += fabs(e) * w; 
                break;
            case LOG: err += log(1+e*e)*w;
                break;
            case LS: err += e*e*w;
                break;
            }
        }
        evaluation_count++;
        return err*0.5;
    }
    
    VectorXd gauss_newton_direction(VectorXd& v, VectorXd& deriv, double& fsse) {
        MatrixXd J(data.size(), v.rows());
        J.setZero();
        fsse = 0; 
        
        VectorXd r(data.size());
        for (size_t m=0; m < data.size(); m++) {
            double w = data[m].weight * data[m].yweight;
            double fx = 0;
            
            J.row(m) = rp_deriv(v, data[m].x*xsf, fx); 
            double e = fx - data[m].y*ysf;
            r[m] = e*w;
            fsse += e*e*w;
        }
        
        
        VectorXd direction = J.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-r);
        
        fsse *= 0.5;
        deriv = J.transpose() * r;
        evaluation_count++;
        return direction;
    }
    
    VectorXd evaluate_derivative(VectorXd& v) {
        // just assume "fit==LS" for now
        
        VectorXd deriv(v.rows());
        deriv.setZero();
        VectorXd d(v.rows());
        for (size_t i=0; i < data.size(); i++) {
            double w = data[i].weight * data[i].yweight;
            double z=0;
            
            d.setZero();
            if (rp_deriv_eval(v, data[i].x*xsf, d, z)) {
                double e = z - data[i].y*ysf;
                deriv += e*d*w;
            } else {
                printf("hit a pole!\n");
            }
        }
        evaluation_count++;
        return deriv;
    }
    
    VectorXd gauss_newton_armijo(VectorXd& v) {
        const double tau = 0.5;
        const double c = 1e-4;
        double fx = 0;
        
        VectorXd grad;
        VectorXd next;
        VectorXd pk;
        for (int k=0; k < 50; k++) {
            
            double alpha = 1.0;
            pk = gauss_newton_direction(v, grad, fx);
            
            double target = fx + c*alpha*pk.dot(grad);
            
            int max_steps = 30;
            next = v + alpha*pk;
            while (evaluate(next) > target && --max_steps > 0) { // iteratively step close until we have a sufficient decrease (Armijo condition)
                target = fx + c*alpha*pk.dot(grad);
                alpha *= tau;
                next = v + alpha*pk;
            }
            
            double stepsize = pk.array().abs().maxCoeff()*fabs(alpha);
            if (stepsize < 5e-8) {
                printf("stopping because improvement is only %le (at iter=%d)\n", stepsize, k);
                break;
            }
            v = next;
        }
        printf("final fx=%lf at %d evals\n", evaluate(v), (int)evaluations());
        return v;
    }
    
    virtual int dimension(void) {
        return (order_n+1 + order_m);
    }
    
    inline double rpeval(const VectorXd& v, double x) {
        double top_val = v[0];
        double p = x;
        for (int n=1; n <= order_n; n++) {
            top_val += p*v[n];
            p *= x;
        }
        double bot_val = base_value;
        p = x;
        for (int m=0; m < order_m; m++) {
            bot_val += p*v[m+order_n+1];
            p *= x;
        }
        
        return top_val / bot_val;
    }
    
    inline bool rp_deriv_eval(const VectorXd& v, double x, VectorXd& d, double& f) {
        // if x falls on a pole, we are in trouble
        // and should probably just return the zero vector?
    
        VectorXd dp = VectorXd::Zero(v.rows());
        VectorXd dq = VectorXd::Zero(v.rows());
        
        double top_val = v[0];
        double p = 1;
        dp[0] = 1;
        for (int n=1; n <= order_n; n++) {
            p *= x;
            dp[n] = p;
            top_val += p*v[n];
        }
        double bot_val = base_value;
        p = 1;
        for (int m=0; m < order_m; m++) {
            p *= x;
            dq[m+order_n+1] = p;
            bot_val += p*v[m+order_n+1];
        }
        
        double den = bot_val*bot_val;
        if (den < 1e-12) {
            d.setZero();
            return false; // return zero derivative at pole
        }
        
        f = top_val / bot_val;
        den = 1.0/den;
        
        d = (dp*bot_val - top_val*dq) * den;
        
        return true;
    }
    
    inline VectorXd rp_deriv(const VectorXd& v, double x, double& f) {
        // if x falls on a pole, we are in trouble
        // and should probably just return the zero vector?
        
        // TODO: we can probably combine this function with rp_deriv_eval ??
    
        VectorXd dp = VectorXd::Zero(v.rows());
        VectorXd dq = VectorXd::Zero(v.rows());
        
        double top_val = v[0];
        double p = 1;
        dp[0] = 1;
        for (int n=1; n <= order_n; n++) {
            p *= x;
            dp[n] = p;
            top_val += p*v[n];
        }
        double bot_val = base_value;
        p = 1;
        for (int m=0; m < order_m; m++) {
            p *= x;
            dq[m+order_n+1] = p;
            bot_val += p*v[m+order_n+1];
        }
        
        double den = bot_val*bot_val;
        if (den < 1e-12) {
            return VectorXd::Zero(v.rows());
        }
        
        f = top_val / bot_val;
        den = 1.0/den;
        
        return (dp*bot_val - top_val*dq) * den;
    }
    
    double peak(const VectorXd& v) {
        double xmin=1e50;
        double xmax=-1e50;
        for (size_t i=0; i < data.size(); i++) {
            xmin = std::min(data[i].x, xmin);
            xmax = std::max(data[i].x, xmax);
        }
        // bracket the maximum
        double peak_z = 0;
        double peak_x = (xmin + xmax)*0.5;
        double step = (xmax - xmin)/20.0;
        for (double x=xmin; x <= xmax; x += step) {
            double z = rpeval(v, x);
            if (z > peak_z) {
                peak_x = x;
                peak_z = z;
            }
        }
        
        // golden section search
        const double phi = 0.61803398874989;
        double lower = peak_x - 2*step;
        double upper = peak_x + 2*step;
        double c = upper - phi*(upper - lower);
        double d = lower + phi*(upper - lower);
        const double tol = 1e-10;
        while ((upper - lower) > tol) {
            double fc = rpeval(v, c);
            double fd = rpeval(v, d);
            if (fc > fd) {
                upper = d;
                d = c;
                c = upper - phi*(upper - lower);
            } else {
                lower = c;
                c = d;
                d = lower + phi*(upper - lower);
            }
        }
        return 0.5*(upper + lower);
    }
    
    const vector<Sample>& get_data(void) const {
        return data;
    }
    
    bool has_poles(const VectorXd& v) {
        double xmin=1e50;
        double xmax=-1e50;
        for (size_t i=0; i < data.size(); i++) {
            xmin = std::min(data[i].x, xmin);
            xmax = std::max(data[i].x, xmax);
        }
        // ensure the bounds are slightly wider than the actual data
        double span=xmax - xmin;
        xmin -= pscale*span;
        xmax += pscale*span;
        
        // now compute roots of bottom polynomial
        switch(order_m) {
        case 0:
            return false; // cannot have poles
        case 1:
            {
                double pole = -1 / v[order_n+1];
                
                if (pole >= xmin && pole <= xmax) {
                    printf("pole at %lf (v=%lf)\n", -1/v[order_n+1], v[order_n+1]);
                }
                
                return pole >= xmin && pole <= xmax;
            }
        case 2:
            {
                double a = v[order_n+2];
                double b = v[order_n+1];
                double c = base_value;
                double sb = b < 0 ? -1 : 1;
                double q = -0.5*(b + sb*sqrt(b*b - 4*a*c));
                double pole1 = q/a;
                double pole2 = c/q;
                
                if ((pole1 >= xmin && pole1 <= xmax) ||
                       (pole2 >= xmin && pole2 <= xmax)) {
                       
                       printf("poles at %lf, %lf\n", pole1, pole2);
                }
                
                return (pole1 >= xmin && pole1 <= xmax) ||
                       (pole2 >= xmin && pole2 <= xmax);
            }
        default:
            // TODO: see NR chapter 5.6 for cubic roots
            printf("Warning: no implementation to compute roots of order-%d polynomial\n",
                order_m
            );
            return false;
        };
    }
    
    const vector<Sample>& data;
    int order_n;
    int order_m;
    double base_value;
    
    double xsf;
    double ysf;
    double pscale;
    fit_type fit;
};

class Poly_fit  {
  public:
  
    Poly_fit(vector< pair<Point, double> >& data, int order_n, int order_m, double cscale,
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
                    
                    Curve_fit cf(sample, order_n, order_m);
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
        
        Curve_fit cf(peak_pts, 2,1, Curve_fit::LS);
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
    
    VectorXd rpfit(Curve_fit& cf, bool scale=false, bool refine=false) {
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
        vector<double> reweights(pts_row.size(), 1.0);
        
        for (int iter=0; iter < 1; iter++) {
            cov.setZero();
            b.setZero();
            a.setZero();
            
            for (size_t i=0; i < pts_row.size(); i++) {
                const Sample& sp = pts_row[i];
                double w = sp.weight * sp.yweight * reweights[i];
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

class Mtf_renderer_mfprofile : public Mtf_renderer {
  public:
    Mtf_renderer_mfprofile(
        const Point& zero, const Point& transverse, const Point& longitudinal,
        double chart_scale,
        const std::string& wdir, const std::string& prof_fname, 
        const std::string& gnuplot_binary,
        const cv::Mat& img, const vector<cv::Point3d>& distance_scale,
        bool lpmm_mode=false, double pixel_size=1.0) 
      :  zero(zero), transverse(transverse), longitudinal(longitudinal),
         wdir(wdir), prname(prof_fname),
         gnuplot_binary(gnuplot_binary), img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         chart_scale(chart_scale), distance_scale(distance_scale),
         gnuplot_failure(false), gnuplot_warning(true) {
      
    }
    
    
    
    void render(const vector<Block>& blocks) {
        Point centroid(0,0);
        
        if (blocks.size() < 10) { // probably not a valid image for profiles
            return;
        }
        
        
        
        vector< pair<Point, double> > points;
        double max_trans = 0;
        vector<double> mtf50;
        //FILE* rf = fopen("rawpoints.txt", "wt");
        for (size_t i=0; i < blocks.size(); i++) {
            const double angle_thresh = 20.0;
            
            for (size_t k=0; k < 4; k++) {
                Point ec = blocks[i].get_edge_centroid(k);
                Point norm = blocks[i].get_normal(k);
                double delta = transverse.x*norm.x + transverse.y*norm.y;
                
                if (acos(fabs(delta))/M_PI*180.0 < angle_thresh) { // edge perp to tangent
                    //fprintf(rf, "%lf %lf %lf\n", ec.x, ec.y, blocks[i].get_mtf50_value(k));
                    points.push_back(make_pair(ec, blocks[i].get_mtf50_value(k)));
                    if (fabs(ec.y - zero.y) > max_trans) {
                        max_trans = fabs(ec.y - zero.y);
                    }
                    mtf50.push_back(blocks[i].get_mtf50_value(k));
                } 
            }
        }
        //fclose(rf);
        
        sort(mtf50.begin(), mtf50.end());
        double p5 = mtf50[0.05*mtf50.size()];
        double p95 = mtf50[0.95*mtf50.size()];
        
        FILE* ap = fopen("axis_projected.txt", "wt");
        vector< pair<Point, double> > ap_points; // axis projected points
        for (size_t i=0; i < points.size(); i++) {
            Point d = points[i].first - zero;
            Point coord(
                d.x*longitudinal.x + d.y*longitudinal.y,
                d.x*transverse.x + d.y*transverse.y
            );
            double val = (points[i].second - p5)/(p95-p5);
            
            // TODO: we can clip the values here to [0,1]
            
            ap_points.push_back(make_pair(coord, val));
            fprintf(ap, "%lf %lf %lf\n", coord.x, coord.y, val);
        }
        fclose(ap);
        
        Poly_fit pf(ap_points, 3, 2, chart_scale, distance_scale);

        /*
        FILE* prfile = fopen((wdir+prname).c_str(), "wt");
        i=0;
        double max_med_filt = 0;
        int max_med_filt_coord = 0;
        for (map<int, double>::const_iterator it = row_max.begin(); it != row_max.end(); ++it) {
            if (med_filt_mtf[i] >= max_med_filt) {
                max_med_filt = med_filt_mtf[i];
                max_med_filt_coord = it->first;
            }
            fprintf(prfile, "%lf %lf %lf\n", 
                it->first/pixel_size, 
                it->second*pixel_size, 
                med_filt_mtf[i++]*pixel_size
            );
        }
        
        
        fprintf(prfile, "\n\n%lf %lf %lf\n", 
            transpose ? 
                 blocks[largest_block].get_edge_centroid(peak_idx).x/pixel_size :
                 blocks[largest_block].get_edge_centroid(peak_idx).y/pixel_size , 
            peak_mtf50*pixel_size,
            peak_mtf50*3*pixel_size
        );
        fclose(prfile);
        */
		       
        FILE* gpf = fopen( (wdir + string("mfprofile.gnuplot")).c_str(), "wt");
        fprintf(gpf, "set xlab \"column (%s)\"\n", lpmm_mode ? "mm" : "pixels");
        fprintf(gpf, "set ylab \"MTF50 (%s)\"\n", lpmm_mode ? "line pairs per mm" : "cycles/pixel");
        fprintf(gpf, "set term png size 1024, 768\n");
        fprintf(gpf, "set output \"%sprofile_image.png\"\n", wdir.c_str());
        fprintf(gpf, "plot [][0:%lf]\"%s\" index 0 u 1:2 t \"MTF50 (%s) raw\" w p ps 0.25, \"%s\" index 0 u 1:3 t \"MTF50 (%s) smoothed\" w l lw 3, \"%s\" index 1 u 1:2 t \"Expected focus point\" w i lc %d lw 3\n", 
            1*pixel_size, // TODO: compute effective max (95%?)
            (wdir+prname).c_str(), lpmm_mode ? "lp/mm" : "c/p",
            (wdir+prname).c_str(), lpmm_mode ? "lp/mm" : "c/p",
            (wdir+prname).c_str(), 3); // supposed to be peak quality
        fclose(gpf);
        
        char* buffer = new char[1024];
        #ifdef _WIN32
        sprintf(buffer, "\"\"%s\" \"%smfprofile.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
        #else
        sprintf(buffer, "\"%s\" \"%smfprofile.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
        #endif
        int rval = system(buffer);
        if (rval != 0) {
            printf("Failed to execute gnuplot (error code %d)\n", rval);
            printf("You can try to execute [%s] to render the plots manually\n", buffer);
            gnuplot_failure = true;
        } else {
            printf("Gnuplot plot completed successfully. Look for profile_image.png\n");
        }
        
        delete [] buffer;
        
    }

    void set_gnuplot_warning(bool gnuplot) {
        gnuplot_warning = gnuplot;
    }
    
    bool gnuplot_failed(void) {
        return gnuplot_failure;
    }

  private:

    Point zero;
    Point transverse;
    Point longitudinal;
    std::string wdir;
    std::string prname;
    std::string pfname;
    std::string gnuplot_binary;
    const cv::Mat& img;
    bool    lpmm_mode;
    double  pixel_size;
    double  chart_scale;
    const vector<cv::Point3d>& distance_scale;
    bool gnuplot_failure;
    bool gnuplot_warning;
};

#endif
