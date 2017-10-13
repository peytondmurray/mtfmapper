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

#include "include/esf_sampler_piecewise_quad.h"
#include <Eigen/Dense>

vector<double> Esf_sampler_piecewise_quad::piecewise_quadfit(const vector<Point2d>& pts) {
    /// assumes pts.size() >= 3 

    const double T = 1.0; // we could find T by looking for the peak in the second difference ??
    
    Eigen::MatrixXd A(pts.size(), 6);
    Eigen::VectorXd b(pts.size());
    Eigen::MatrixXd C(2, 6);
    
    A.setZero();
    
    // constraints: first row is f(T) = g(T), second row is f'(T) = g'(T)
    C << T*T, T, 1, -T*T, -T, -1,
         2*T, 1, 0, -2*T, -1,  0;
         
    for (size_t r=0; r < pts.size(); r++) {
        if (pts[r].x <= T) {
            A(r, 0) = pts[r].x*pts[r].x;
            A(r, 1) = pts[r].x;
            A(r, 2) = 1.0;
        } else {
            A(r, 3) = pts[r].x*pts[r].x;
            A(r, 4) = pts[r].x;
            A(r, 5) = 1.0;
        }
        b[r] = pts[r].y;
    }
    
    Eigen::MatrixXd design(6 + 2, 6 + 2);
    design.setZero();
    design.topLeftCorner(6, 6) = A.transpose() * A;
    design.bottomLeftCorner(2, 6) = C;
    design.topRightCorner(6, 2) = C.transpose();
    
    Eigen::VectorXd b_d(6 + 2);
    b_d.setZero();
    b_d.topRows(6) = A.transpose() * b;
    
    Eigen::VectorXd sol = design.colPivHouseholderQr().solve(b_d);
    
    // sol[0-5] : quadratic parameters
    // sol[6-7] : constraints values (not used)
    
    for (size_t r=0; r < pts.size(); r++) {
        double pred = 0;
        if (pts[r].x <= T) {
            pred = sol[0]*pts[r].x*pts[r].x + sol[1]*pts[r].x + sol[2];
        } else {
            pred = sol[3]*pts[r].x*pts[r].x + sol[4]*pts[r].x + sol[5];
        }
    }
    
    vector<double> parms(7);
    parms[0] = T;
    parms[1] = sol[0];
    parms[2] = sol[1];
    parms[3] = sol[2];
    parms[4] = sol[3];
    parms[5] = sol[4];
    parms[6] = sol[5];
    return parms;
}

void Esf_sampler_piecewise_quad::sample(Edge_model& edge_model, vector<Ordered_point>& local_ordered, 
    const map<int, scanline>& scanset, double& edge_length,
    const cv::Mat& geom_img, const cv::Mat& sampling_img) {
    
    double max_along_edge = -1e50;
    double min_along_edge = 1e50;
    
    vector<Point2d> local_pts;
    for (auto p: edge_model.ridge) {
        Point2d delta = p - edge_model.get_centroid();
        local_pts.push_back(Point2d(delta.ddot(edge_model.get_direction()), delta.ddot(edge_model.get_normal())));
    }
    
    vector<double> qpp = piecewise_quadfit(local_pts);
    std::array<double, 3> qp;
        
    for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
        int y = it->first;
        if (y < border_width || y > geom_img.rows-1-border_width) continue;
        int rowcode = (y & 1) << 1;
        
        for (int x=it->second.start; x <= it->second.end; ++x) {
            
            if (x < border_width || x > geom_img.cols-1-border_width) continue;
            
            int code = 1 << ( (rowcode | (x & 1)) ^ 3 );
            if ((code & cfa_mask) == 0) continue;
            
            Point2d d = Point2d(x, y) - edge_model.get_centroid();
            double perp = d.ddot(edge_model.get_normal()); 
            double par = d.ddot(edge_model.get_direction());
            
            if (par < qpp[0]) {
                qp[0] = qpp[1];
                qp[1] = qpp[2];
                qp[2] = qpp[3];
            } else {
                qp[0] = qpp[4];
                qp[1] = qpp[5];
                qp[2] = qpp[6];
            }
            
            // (par, perp) is in the local coordinate frame of the parabola qp
            // so find the closest point on qp
            vector<double> roots = quad_tangency(Point2d(par, perp), qp);
            perp = 1e20;
            for (auto r: roots) {
                Point2d on_qp(r, r*r*qp[0] + r*qp[1] + qp[2]);
                Point2d recon = on_qp.x*edge_model.get_direction() + on_qp.y*edge_model.get_normal() + edge_model.get_centroid();
                Point2d delta = Point2d(x, y) - recon;
                double dist = std::copysign(norm(delta), delta.ddot(edge_model.get_normal()));
                if (fabs(dist) < fabs(perp)) {
                    perp = dist;
                }
            }
            
            if (fabs(perp) < max_dot) {
                local_ordered.push_back(Ordered_point(perp, sampling_img.at<uint16_t>(y,x) ));
                max_along_edge = max(max_along_edge, par);
                min_along_edge = min(min_along_edge, par);
                
                /*
                if ((x&1) || (y&1)) {
                    cv::Vec3b& color = od_img.at<cv::Vec3b>(y, x);
                    color[0] = 255;
                    color[1] = 0;
                    color[2] = 0;
                }
                */
            }
        }
    }
    
    edge_length = max_along_edge - min_along_edge;
}

