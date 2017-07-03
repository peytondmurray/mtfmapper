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
#ifndef MTF_CORE_H
#define MTF_CORE_H

#include "include/logger.h"
#include "include/component_labelling.h"
#include "include/gradient.h"
#include "include/block.h"
#include "include/rectangle.h"
#include "include/ellipse.h"
#include "include/edge_record.h"
#include "include/loess_fit.h"
#include "include/afft.h"
#include "include/mtf_profile_sample.h"
#include "include/bayer.h"
#include "include/undistort.h"

#include <map>
using std::map;

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

typedef vector<Block> block_vector;

// global constants for ESF-fourier MTF method
// TODO: these can be dynamic parameters, with some effort
const double max_dot = 28;
const int SAMPLES_PER_PIXEL = 32;
const size_t FFT_SIZE = int(16)*SAMPLES_PER_PIXEL;
const int NYQUIST_FREQ = FFT_SIZE/16;
const double max_edge_length = 200;

class Mtf_core {
  public:
    

    Mtf_core(const Component_labeller& in_cl, const Gradient& in_g, 
             const cv::Mat& in_img, const cv::Mat& in_bayer_img, std::string bayer_subset)
      : cl(in_cl), g(in_g), img(in_img), bayer_img(in_bayer_img), absolute_sfr(false),
        snap_to(false), snap_to_angle(0), sfr_smoothing(true),
        sliding(false), samples_per_edge(0), border_width(0), find_fiducials(false) {

        bayer = Bayer::from_string(bayer_subset);
        logger.debug("bayer subset is %d\n", bayer);
      
        for (Boundarylist::const_iterator it=cl.get_boundaries().begin(); it != cl.get_boundaries().end(); ++it) {
            valid_obj.push_back(it->first);
        }
        
        cv::Mat temp;
        in_img.convertTo(temp, CV_8U, 256.0/16384.0);
        cv::cvtColor(temp, od_img, CV_GRAY2RGB);
    }
    
    ~Mtf_core(void) {
        //cv::imwrite(string("detections.png"), od_img);
    }
    
    size_t num_objects(void) {
        return valid_obj.size();
    }
    
    void search_borders(const Point2d& cent, int label);
    bool extract_rectangle(const Point2d& cent, int label, Mrectangle& rect);
    double compute_mtf(const Point2d& in_cent, const map<int, scanline>& scanset, 
                       Edge_record& er,
                       double& poor, 
                       vector<double>& sfr, vector<double>& esf, 
                       vector<Point2d>& ridge,
                       bool allow_peak_shift = false);
    
    vector<Block>& get_blocks(void) {
        // make a copy into an STL container if necessary
        if (detected_blocks.size() == 0) {
            for (map<int,Block>::const_iterator it=shared_blocks_map.begin();
                 it != shared_blocks_map.end(); ++it) {
                
                bool allzero = true;
                for (int k=0; k < 4 && allzero; k++) {
                    if (fabs(it->second.get_mtf50_value(k)) > 1e-6) {
                        allzero = false;
                    }
                }
                
                if (it->second.valid && !allzero) { 
                    detected_blocks.push_back(it->second);
                }
            }
        }
        return detected_blocks;
    }
    
    vector<Mtf_profile_sample>& get_samples(void) {
        return samples;
    }
    
    void set_absolute_sfr(bool val) {
        absolute_sfr = val;
    }
    
    void set_sfr_smoothing(bool val) {
        sfr_smoothing = val;
    }
    
    void set_sliding(bool val) {
        sliding = val;
        find_fiducials = true;
    }
    
    void set_samples_per_edge(int s) {
        samples_per_edge = s;
    }
    
    void set_snap_angle(double angle) {
        snap_to = true;
        snap_to_angle = angle;
    }
    
    void set_border(int in_border_width) {
        border_width = in_border_width;
    }
    
    void set_find_fiducials(bool val) {
        find_fiducials = val;
    }
    
    void set_undistort(Undistort* u) {
        undistort = u;
    }
    
    void process_image_as_roi(void);
    
    const Component_labeller& cl;
    const Gradient&           g;
    const cv::Mat&            img;
    const cv::Mat&            bayer_img;
    Bayer::bayer_t bayer;
    
    AFFT<512> afft; // FFT_SIZE = 512 ??
    vector<int> valid_obj;
    
    vector<Block> detected_blocks;  
    map<int, Block> shared_blocks_map;
    vector<Point2d> solid_ellipses;
    vector<Ellipse_detector> ellipses;
    
    vector<Mtf_profile_sample> samples;
    
    cv::Mat od_img;
  private:
    bool absolute_sfr;
    bool snap_to;
    double snap_to_angle;
    bool sfr_smoothing;
    bool sliding;
    int samples_per_edge;
    int border_width;
    bool find_fiducials;
    Undistort* undistort = nullptr;
    
    void process_with_sliding_window(Mrectangle& rrect);
    
    void extract_ridge(map< int, pair<double, double> >& edge_residual, 
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
        double edgelen = std::min(fabs(edge_residual.begin()->first), fabs(it1->first));
        
        // now fit a quadratic polynomial to smooth out the initial centerline
        Eigen::Matrix3d cov;
        Eigen::Vector3d yv;
        cov.setZero();
        yv.setZero();
        for (const auto& ee: edge_residual) {
            double perp = ee.second.second;
            double par = ee.first;
            
            if (fabs(par) < edgelen - 1) {
            
                cov(0,0) += 1;
                cov(0,1) += par;
                cov(0,2) += par*par;
                cov(1,2) += par*par*par;
                cov(2,2) += par*par*par*par;
                
                yv(0,0) += perp;
                yv(1,0) += par*perp;
                yv(2,0) += par*par*perp;
            }
        }
        
        // complete cov matrix
        cov(1,0) = cov(0,1);
        cov(1,1) = cov(0,2);
        cov(2,0) = cov(0,2);
        cov(2,1) = cov(1,2);
        
        Eigen::LDLT<Eigen::Matrix3d> qr(cov); // LDLT is probably good enough 
        Eigen::Vector3d sol = qr.solve(yv);
        
        for (auto& ee: edge_residual) {
            double x = ee.first;
            ee.second.second = sol[0] + sol[1]*x + sol[2]*x*x;
        }
        
        // now estimate the centroid at each bin, using the gradient magnitude as weight
        ridge = vector<Point2d>();
        for (auto& ee: edge_residual) {
            double pd = ee.second.second;
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
                    
                    
                    if (fabs(pardot) > 5) {  
                        end_weight = 0;
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
    
    inline void update_gradient_peak(map< int, pair<double, double> >& edge_residual, 
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
    
    void sample_at_angle(double ea, vector<Ordered_point>& local_ordered, 
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
            
            // first construct a new scanset
            Point2d distorted_cent = undistort->transform_pixel(cent.x, cent.y);
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
                    
                    // find bounds on ROI in undistorted image 
                    Point2d d(x - cent.x, y - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                    }
                }
            }
            
            // curve is ordered by dist_along_edge, so we can binary search later
            
            vector< Ordered_point > curve;
            if (!undistort->rectilinear_equivalent()) {
                for (double dperp = lrint(min_along_edge) - 1; dperp < lrint(max_along_edge) + 1; dperp += 1.0) {
                    Point2d sample = dperp * edge_direction + cent;
                    Point2d p = undistort->transform_point(sample.x, sample.y);
                    Point2d d(p - distorted_cent);
                    
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    curve.push_back(Ordered_point(dist_along_edge, dot));
                }
            }
            
            if (bayer == Bayer::NONE) {
                // now visit each pixel in the distorted image, but use undistorted coordinates for ESF construction
                for (map<int, scanline>::const_iterator it=m_scanset.begin(); it != m_scanset.end(); ++it) {
                    int y = it->first;
                    if (y < border_width || y > img.rows-1-border_width) continue;
                    
                    for (int x=it->second.start; x <= it->second.end; ++x) {
                        
                        if (x < border_width || x > img.cols-1-border_width) continue;
                        
                        // transform warped pixel to idealized rectilinear
                        cv::Point2d tp = undistort->inverse_transform_point(x, y);
                        
                        bool inside_image = lrint(tp.x) >= 0 && lrint(tp.x) < g.width() && lrint(tp.y) >= 0 && lrint(tp.y) < g.height();
                        if (!inside_image) continue;
                        
                        double dot = 0;
                        double dist_along_edge = 0;
                        if (undistort->rectilinear_equivalent()) {
                            /*
                            cv::Vec3b& color = od_img.at<cv::Vec3b>(lrint(tp.y), lrint(tp.x));
                            color[0] = 255;
                            color[1] = 255;
                            color[2] = 0;
                            */
                            
                            Point2d d(tp.x - cent.x, tp.y - cent.y);
                            dot = d.ddot(mean_grad); 
                            dist_along_edge = d.ddot(edge_direction);
                        } else {
                            Point2d d(x - distorted_cent.x, y - distorted_cent.y);
                            dot = d.ddot(mean_grad); 
                            dist_along_edge = d.ddot(edge_direction);
                            
                            // now we have to find the closest sample
                            int fidx = lower_bound(curve.begin(), curve.end(), dist_along_edge) - curve.begin();
                            
                            // curve is the ridge line of the edge
                            // simply subtract the ridge for now
                            dot -= curve[fidx].second;
                        }
                        if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                            local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) )); // TODO: hack --- we are abusing the bayer image?
                            
                            // (x, y) are equi-angular, but we need the rectilinear
                            // coordinates to go with the gradient image (which is rectilinear)
                            // TODO: should we use rounding, or floor here ?
                            update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(lrint(tp.x), lrint(tp.y)));
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
                        
                        update_gradient_peak(edge_residual, dist_along_edge, dot, g.grad_magnitude(lrint(tp.x), lrint(tp.y)));

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
                        
                        cv::Vec3b& color = od_img.at<cv::Vec3b>(lrint(tp.y), lrint(tp.x));
                        color[0] = 255;
                        color[1] = 255;
                        color[2] = 0;
                        
                        if (!undistort->rectilinear_equivalent()) {
                            Point2d d(x - distorted_cent.x, y - distorted_cent.y);
                            dot = d.ddot(mean_grad); 
                            dist_along_edge = d.ddot(edge_direction);
                            
                            // now we have to find the closest sample
                            int fidx = lower_bound(curve.begin(), curve.end(), dist_along_edge) - curve.begin();
                            
                            // curve is the ridge line of the edge
                            // simply subtract the ridge for now
                            dot -= curve[fidx].second;
                        }
                        if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                            local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) ));
                        }
                    }
                }
            }
        }
    
        extract_ridge(edge_residual, cent, mean_grad, edge_direction, ridge);
        edge_length = max_along_edge - min_along_edge;
    }

};

#endif
