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

#ifndef DISTORTION_OPTIMIZER_H
#define DISTORTION_OPTIMIZER_H

#include "Eigen/Dense"
#include <random>
#include <limits>
#include <cmath>

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "include/block.h"

class Ridge {
  public:
    Ridge(const vector<Point2d>& ridge, const Point2d& centroid, const Point2d& normal, double weight=0) 
    : ridge(ridge), centroid(centroid), normal(normal), weight(weight) {}
    
    vector<Point2d> ridge;
    Point2d centroid;
    Point2d normal;
    double weight;
    double residual = 0;
};

class Distortion_optimizer {
  public:
    Distortion_optimizer(const vector<Block>& in_blocks, const Point2d& prin) 
    : prin(prin), max_val(1,1) {
        
        radius_norm = sqrt(prin.x*prin.x + prin.y*prin.y); // corners are at radius 1
        
        // weight each edge according to its orientation relative to the principal point
        for (const auto& block: in_blocks) {
            for (int k=0; k < 4; k++) {
            
                if (block.get_ridge(k).size() > 2) {
                
                    Point2d cent = block.get_edge_centroid(k);
                    Point2d normal = block.get_normal(k);
                    Point2d radial(cent.x - prin.x, cent.y - prin.y);
                    if (norm(radial) > 1e-6) {
                        radial *= 1.0/(norm(radial));
                        double w = fabs(radial.ddot(normal));
                        if (w > 0.17365) { // edge angle < 80 degrees w.r.t. radial direction
                            ridges.push_back(Ridge(block.get_ridge(k), block.get_edge_centroid(k), block.get_normal(k), w));
                            double wrad = norm(cent - prin) / radius_norm;
                            maxrad = std::max(wrad, maxrad);
                            max_val.x = std::max(max_val.x, fabs(cent.x - prin.x));
                            max_val.y = std::max(max_val.y, fabs(cent.y - prin.y));
                        }
                    }
                }
            }
        }
        
        logger.debug("Kept %ld of %ld edges\n", ridges.size(), 4*in_blocks.size());
        
        // pack initial parameters
        Eigen::VectorXd init(2);
        
        init[0] = 0;
        init[1] = 0;
        
        best_sol = initial = init;
    }
    
    Point2d get_max_val(void) const {
        return max_val;
    }
    
    void solve(void) {
        Eigen::VectorXd scale(best_sol.size());
        
        scale << 0.4995, 0.1665; // just inside valid region
                
        double initial_err = evaluate(best_sol, 0.0);
        logger.debug("initial distortion rmse: %lf\n", initial_err);
        
        nelder_mead_failed = false;
        seed_simplex(best_sol, scale);
        best_sol = iterate(1e-12);
        
        // now compute outliers ...
        // TODO: should take weighting into account (no need to use perpendicular edges in outlier calculation)
        vector<float> residuals;
        for (const auto& edge: ridges) {
            residuals.push_back(edge.residual);
        }
        
        if (residuals.size() > 4) {
        
            double mc = medcouple(residuals); // side effect: residuals is now sorted
            double outlier_threshold = residuals[residuals.size()*0.75] + 1.5*exp(4*mc)*(residuals[residuals.size()*0.75] - residuals[residuals.size()*0.25]);
            size_t outlier_count = 0;
            for (auto& edge: ridges) {
                if (edge.residual >= outlier_threshold) {
                    outlier_count++;
                    edge.weight = 0;
                }
            }
            logger.debug("medcouple=%lf, outlier threshold=%lf, min=%f, median=%f, max=%f\n", 
                mc, outlier_threshold, residuals.front(), residuals[residuals.size()/2], residuals.back()
            );
            logger.debug("Found %ld outliers out of %ld values\n", outlier_count, ridges.size());
            
            
            if (outlier_count > 0 && ridges.size() - outlier_count > ridges.size()*0.7) {
                logger.debug("restarting after outlier suppression:\n");
                best_sol = initial;
                seed_simplex(best_sol, scale);
                best_sol = iterate(1e-12);
            }
        }
        
        double final_err = evaluate(best_sol, 0.0);
        logger.debug("final distortion rmse: %lf\n", final_err);
    }
    
    // two-parameter division model; also works when k2(v[1]) == 0
    Point2d inv_warp(const Point2d& p, const Eigen::VectorXd& v) {
        double px = (p.x - prin.x);
        double py = (p.y - prin.y);
    
        const double rd = sqrt((px)*(px) + (py)*(py)) / radius_norm; 
        const double r2 = rd*rd;
        double ru = 1 + (v[0] + v[1]*r2)*r2;
        
        px = px/ru + prin.x;
        py = py/ru + prin.y;
        
        return cv::Point2d(px, py);
    }
    
    double model_not_invertible(const Eigen::VectorXd& v) {
        const double k1 = v[0];
        const double k2 = v[1];
        const double r1 = 1.0;
        const double r12 = r1*r1;
        const double r14 = r12*r12;
        
        if (k1*r12 <= -2) return 1.0;
        if (k1*r12 < 2) {
            if (k2*r14 <= -1 - k1*r12) return 1.0;
            if (k2*r14 >= (1 - k1*r12)/3.0) return 1.0;
            return 0.0;
        }
        // k1r12 >= 2
        if (k2*r14 <= -1 - k1*r12) return 1.0;
        if (k2*r14 >= -k1*k1*r14/12.0) return 1.0;
        return 0.0;
    }
    
    double eval_d2(const double x, const double a, const double b) {
        const double x2 = x*x;
        const double x3 = x2*x;
        const double x4 = x2*x2;
        double bottom = a*x2 + b*x4 + 1;
        bottom *= bottom*bottom;
        return (2*a*a*x3 + 6*a*x*(b*x4 - 1) + 4*b*x3*(3*b*x4 - 5)) / bottom;
    }
    
    double d2_sign_change(const Eigen::VectorXd& v, double maxrad) {
        const double r0 = 1e-4;
        double prev = eval_d2(r0, v[0], v[1]);
        double changes = 0;
        for (double rd=r0; rd < 1.05*maxrad; rd += 0.01) {
            double cur = eval_d2(rd, v[0], v[1]);
            if (prev * cur < 0) {
                changes++;
                prev = cur;
            }
        }
        return changes;
    }
    
    double medcouple(vector<float>& x) {
        
        sort(x.begin(), x.end());
        
        double xm = x[x.size()/2];
        double xscale = 2*x.back();
        
        vector<double> zplus;
        vector<double> zminus;
        
        for (size_t i=x.size()-1; i > 0 && x[i] >= xm; i--) {
            zplus.push_back((x[i] - xm)/xscale);
        }
        for (size_t i=0; i < x.size() && x[i] <= xm; i++) {
            zminus.push_back((x[i] - xm)/xscale);
        }
        
        size_t p = zplus.size();
        size_t q = zminus.size();
        
        vector<float> h(q*p);
        for (size_t i=0; i < p; i++) {
            for (size_t j=0; j < q; j++) {
                
                float& hval = h[j*p + i];
                
                if (zplus[i] == zminus[j]) {
                    int64_t dv = int64_t(p) - 1 - int64_t(i) - int64_t(j);
                    hval = dv < 0 ? -1 : (dv == 0 ? 0 : 1);
                } else {
                    hval = (zplus[i] + zminus[j]) / (zplus[i] - zminus[j]);
                }
            }
        }
        
        nth_element(h.begin(), h.begin()+h.size()/2, h.end());
        return h[h.size()/2];
    }
    
    double evaluate(const Eigen::VectorXd& v, double penalty=1.0) {
        // TODO: add support for autocropped images
        
        double count = 0;
        double merr = 0;
        for (auto& edge: ridges) {
            
            Point2d cent = inv_warp(edge.centroid, v);
            Point2d dir(-edge.normal.y, edge.normal.x);
            
            vector<Point2d> ridge(edge.ridge.size());
            
            // estimate edge length to establish scale
            Point2d first_point = inv_warp(edge.ridge.front(), v);
            Point2d last_point = inv_warp(edge.ridge.back(), v);
            double scale = std::max(fabs((first_point - cent).ddot(dir)), fabs((last_point - cent).ddot(dir)));
            
            Eigen::Matrix3d cov;
            Eigen::Vector3d yv;
            cov.setZero();
            yv.setZero();
            for (size_t ri=0; ri < edge.ridge.size(); ri++) {
                
                Point2d p = inv_warp(edge.ridge[ri], v);
                
                Point2d d = p - cent;
                double perp = d.ddot(edge.normal); 
                double par = d.ddot(dir) / scale;
                
                // keep the refined values
                ridge[ri] = Point2d(par, perp);
                
                cov(0,0) += 1;
                cov(0,1) += par;
                cov(0,2) += par*par;
                cov(1,2) += par*par*par;
                cov(2,2) += par*par*par*par;
                
                yv(0,0) += perp;
                yv(1,0) += par*perp;
                yv(2,0) += par*par*perp;
            }
            
            // complete cov matrix
            cov(1,0) = cov(0,1);
            cov(1,1) = cov(0,2);
            cov(2,0) = cov(0,2);
            cov(2,1) = cov(1,2);
            
            Eigen::LDLT<Eigen::Matrix3d> qr(cov); // LDLT is probably good enough 
            Eigen::Vector3d sol = qr.solve(yv);
            
            // for now, do another pass to calculate the standard error
            // we can probably speed this up by using the fit directly
            double ressum = 0;
            for (size_t ri=0; ri < ridge.size(); ri++) {
                double sx = ridge[ri].x;
                double res = (sol[0] + sx*sol[1] + sx*sx*sol[2]) - ridge[ri].y;
                ressum += res*res;
            }
            
            ressum /= ridge.size() - 2;
            Eigen::Matrix3d icov = cov.inverse();
            double t_const  = fabs(sol[0] / sqrt(icov(0,0)*ressum));
            double t_linear = fabs(sol[1] / sqrt(icov(1,1)*ressum));
            double t_quad = fabs(sol[2] / sqrt(icov(2,2) * ressum));
            double t0 = t_quad / std::max(t_const, t_linear);;
            
            if (std::isfinite(t0) && !std::isnan(t0)) {
                edge.residual = t0;
                merr += t0 * edge.weight;
                count += edge.weight;
            } else {
                edge.residual = 1e6;
            }
            
        }
        
        if (penalty == 0) {
            logger.debug("maxrad = %lf k1 = %lf k2 = %lf sign changes = %lf\n", maxrad, v[0], v[1], d2_sign_change(v, maxrad));
        }
        
        merr /= count;
        return merr + penalty*(model_not_invertible(v)*1e4 + merr*d2_sign_change(v, maxrad)/100.0 + merr*sqrt(std::min(fabs(v[0]),fabs(v[1])))/100.0);
        //return merr + penalty*(model_not_invertible(v)*1e4);
    }
    
    void seed_simplex(Eigen::VectorXd& v, const Eigen::VectorXd& lambda) {
        np = vector<Eigen::VectorXd>(v.size()+1);
        // seed the simplex
        for (int i = 0; i < v.size(); i++) {
            np[i] = v;
            np[i][i] += lambda[i];
        }
        np[v.size()] = v;

        ny = Eigen::VectorXd(v.size()+1);
        // now obtain their function values
        for (int i = 0; i < v.size() + 1; i++) {
            ny[i] = evaluate(np[i]);
        }
    }
    
    inline void simplex_sum(Eigen::VectorXd& psum) {
        psum.setZero();
        for (size_t m=0; m < np.size(); m++) {
            psum += np[m];
        }
    }
    
    void nelder_mead(const double ftol, int& num_evals) {
        const int max_allowed_iterations = 5000;
        const double epsilon = 1.0e-10;

        Eigen::VectorXd psum(np[0].size());
        num_evals = 0;
        simplex_sum(psum);
        
        for (;;) {
            size_t inhi;
            size_t ilo = 0;
            size_t ihi = ny[0] > ny[1] ? (inhi = 1, 0) : (inhi = 0, 1);

            for (size_t i=0; i < np.size(); i++) {
                if (ny[i] <= ny[ilo]) {
                    ilo = i;
                }
                if (ny[i] > ny[ihi]) {
                    inhi = ihi;
                    ihi = i;
                } else
                if (ny[i] > ny[inhi] && i != ihi) {
                    inhi = i;
                }
            }
            double rtol = 2.0 * fabs(ny[ihi] - ny[ilo]) / ( fabs(ny[ihi]) + fabs(ny[ilo]) + epsilon );
            if (rtol < ftol) {
                std::swap(ny[0], ny[ilo]);
                for (size_t i=0; i < (size_t)np[0].size(); i++) {
                    std::swap(np[0][i], np[ilo][i]);
                }
                break;
            }
            if (num_evals >= max_allowed_iterations) {
                nelder_mead_failed = true;
                return;
            }
            num_evals += 2;
            double ytry = try_solution(psum, ihi, -1.0);
            if (ytry <= ny[ilo]) {
                ytry = try_solution(psum, ihi, 2.0);
            } else {
                if (ytry >= ny[inhi]) {
                    double ysave = ny[ihi];
                    ytry = try_solution(psum, ihi, 0.5);
                    if (ytry >= ysave) {
                        for (size_t i=0; i < np.size(); i++) {
                            if (i != ilo) {
                                np[i] = psum = (np[i] + np[ilo]) * 0.5;
                                ny[i] = evaluate(psum);
                            }
                        }
                        num_evals += np[0].size();
                        simplex_sum(psum);
                    }
                } else {
                    num_evals--;
                }
            }
        }
    }
    
    double try_solution(Eigen::VectorXd& psum, const int ihi, const double fac) {

        double fac1 = (1.0 - fac) / double (psum.size());
        double fac2 = fac1 - fac;
        Eigen::VectorXd ptry = psum * fac1 - np[ihi] * fac2;
        double ytry = evaluate(ptry);

        if (ytry < ny[ihi]) {
            ny[ihi] = ytry;
            psum += ptry - np[ihi];
            np[ihi] = ptry;
        }
        return ytry;
    }
    
    Eigen::VectorXd iterate(double tol) {
    
        int evals = 0;
        const int tries = 2;
        for (int i = 0; i < tries; i++) {
            nelder_mead(tol, evals);
        }
        int min_idx = 0;
        for (size_t i=0; i < (size_t)ny.size(); i++) {
            if (ny[i] < ny[min_idx]) {
                min_idx = i;
            }
        }
        return np[min_idx];
    }
    
    bool optimization_failure(void) const {
        return nelder_mead_failed;
    }
    
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    vector<Ridge> ridges;
    Point2d prin;
    double radius_norm;
    Eigen::VectorXd best_sol;
    Point2d max_val;
    double maxrad;
    
    // variables used by nelder-mead
    vector<Eigen::VectorXd> np;
    Eigen::VectorXd ny;
    bool nelder_mead_failed;
    Eigen::VectorXd initial;
    double focal_lower;
    double focal_upper;
    double focal_mode_constraint;
};

#endif
