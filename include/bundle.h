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

#ifndef BUNDLE_H
#define BUNDLE_H

#include "Eigen/Dense"
#include <random>
#include <limits>
#include <cmath>


class Bundle_adjuster {
  public:
    Bundle_adjuster(
        vector<Eigen::Vector2d>& img_points,
        vector<Eigen::Vector3d>& world_points,
        Eigen::Vector3d& t,
        cv::Mat in_rod_angles,
        double distortion,
        double w,
        double px, double py) 
         : img_points(img_points), world_points(world_points) {
        
        rot_mat = cv::Mat(3, 3, CV_64FC1);
        rod_angles = cv::Mat(3, 1, CV_64FC1);
        
        // pack initial parameters
        Eigen::VectorXd init(10);
        
        init[0] = t[0];
        init[1] = t[1];
        init[2] = t[2];
        init[3] = in_rod_angles.at<double>(0,0);
        init[4] = in_rod_angles.at<double>(1,0);
        init[5] = in_rod_angles.at<double>(2,0);
        init[6] = distortion;
        
        init[7] = w;
        init[8] = px;
        init[9] = py;
        
        #if 1
        FILE* bal = fopen("bal.txt", "wt");
        fprintf(bal, "1 %lu %lu\n", world_points.size(), world_points.size());
        for (size_t i=0; i < img_points.size(); i++) {
            fprintf(bal, "%d %lu %le %le\n", 0, i, img_points[i][0], -img_points[i][1]);
        }
        fprintf(bal, "%le\n%le\n%le\n", init[3], init[4], init[5]);
        fprintf(bal, "%le\n%le\n%le\n", init[0], init[1], init[2]);
        fprintf(bal, "%le\n%le\n%le\n", 1.0/init[7], init[6], 0.0);
        for (size_t i=0; i < world_points.size(); i++) {
            fprintf(bal, "%le\n%le\n%le\n", world_points[i][0], world_points[i][1], world_points[i][2]);
        }
        fclose(bal);
        #endif
        
        best_sol = solve_random_search(init);
    }
    
    void unpack(Eigen::Matrix3d& R, Eigen::Vector3d& t, double& distortion, double& w, double& px, double &py) {
        Eigen::VectorXd& v = best_sol;
        rod_angles.at<double>(0,0) = v[3];
        rod_angles.at<double>(1,0) = v[4];
        rod_angles.at<double>(2,0) = v[5];
        cv::Rodrigues(rod_angles, rot_mat);
        
        // global_scale * K * R | t
        R << rot_mat.at<double>(0,0), rot_mat.at<double>(0,1), rot_mat.at<double>(0,2),
             rot_mat.at<double>(1,0), rot_mat.at<double>(1,1), rot_mat.at<double>(1,2),
             rot_mat.at<double>(2,0), rot_mat.at<double>(2,1), rot_mat.at<double>(2,2);
        
        t = Eigen::Vector3d(v[0], v[1], v[2]);
        
        distortion = v[6];
        w = v[7];
        px = v[8];
        py = v[9];
    }
    
    double evaluate(const Eigen::VectorXd& v) {
        rod_angles.at<double>(0,0) = v[3];
        rod_angles.at<double>(1,0) = v[4];
        rod_angles.at<double>(2,0) = v[5];
        cv::Rodrigues(rod_angles, rot_mat);
        
        // global_scale * K * R | t
        Eigen::Matrix3d R;
        R << rot_mat.at<double>(0,0), rot_mat.at<double>(0,1), rot_mat.at<double>(0,2),
             rot_mat.at<double>(1,0), rot_mat.at<double>(1,1), rot_mat.at<double>(1,2),
             rot_mat.at<double>(2,0), rot_mat.at<double>(2,1), rot_mat.at<double>(2,2);
             
        // we could orthogonalize R here using [U S V] = svd(R) -> Ro = U*V'
        
        R.row(2) *= v[7]; // multiply by 1/f
        Eigen::Vector3d t(v[0], v[1], v[2]*v[7]);
        
        double bpsse = 0;
        for (size_t i=0; i < world_points.size(); i++) {
            Eigen::Vector3d bp = R*world_points[i] + t;
            bp /= bp[2];
            
            bp[0] += v[8];
            bp[1] += v[9];
            
            double rad = 1 + v[6]*(bp[0]*bp[0] + bp[1]*bp[1]);
            bp /= rad;
                        
            double bpe = (img_points[i] - Eigen::Vector2d(bp[0], bp[1])).norm();
            bpsse += bpe*bpe;
        }
        return sqrt(bpsse/double(world_points.size()));
    }
    
    double backproj(const Eigen::VectorXd& v) {
        rod_angles.at<double>(0,0) = v[3];
        rod_angles.at<double>(1,0) = v[4];
        rod_angles.at<double>(2,0) = v[5];
        cv::Rodrigues(rod_angles, rot_mat);
        
        // global_scale * K * R | t
        Eigen::Matrix3d R;
        R << rot_mat.at<double>(0,0), rot_mat.at<double>(0,1), rot_mat.at<double>(0,2),
             rot_mat.at<double>(1,0), rot_mat.at<double>(1,1), rot_mat.at<double>(1,2),
             rot_mat.at<double>(2,0), rot_mat.at<double>(2,1), rot_mat.at<double>(2,2);
        
        R.row(2) *= v[7]; // multiply by 1/f
        Eigen::Vector3d t(v[0], v[1], v[2]*v[7]);
        
        double f5 = 0;
        double bpsse = 0;
        for (size_t i=0; i < world_points.size(); i++) {
            Eigen::Vector3d bp = R*world_points[i] + t;
            
            bp /= bp[2];
            
            bp[0] += v[8];
            bp[1] += v[9];
            
            double rad = 1 + v[6]*(bp[0]*bp[0] + bp[1]*bp[1]);
            bp /= rad;
            
            printf("pt[%lu], expect (%lf %lf), got (%lf %lf)\n", i,
                img_points[i][0], img_points[i][1],
                bp[0], bp[1]
            );
            double bpe = (img_points[i] - Eigen::Vector2d(bp[0], bp[1])).norm();
            bpsse += bpe;
            if (i < 5) f5 += bpe;
        }
        printf("MAD error : %lf, first 5: %lf\n", bpsse/double(world_points.size()), f5/double(world_points.size()));
        return sqrt(bpsse);
    }
    
    Eigen::VectorXd solve_random_search(const Eigen::VectorXd& init) {
        Eigen::VectorXd scale(init.rows());
        scale << 0.01, 0.01, 0.01,    // origin
                 0.001, 0.001, 0.001, // angles
                 1e-5,                // distortion
                 1e-2,                // focal length
                 10, 10;            // principal point
        
        Eigen::VectorXd p = init;
        Eigen::VectorXd pp = init;
        Eigen::VectorXd vel = Eigen::VectorXd::Zero(init.rows());
        Eigen::VectorXd dir = Eigen::VectorXd::Zero(init.rows());
        
        std::mt19937 mt(100);
        std::uniform_real_distribution<double> dist(-1, 1);
        
        double fbest = evaluate(p);
        
        backproj(p);
        printf("initial f=%lf, distortion=%le, delta(%lf, %lf), RMSE=%lf\n", 
            1.0/init[7], init[6], init[8], init[9], fbest
        );
        
        double rho = 0.1;
        double inertia = 0.72;
        
        int count_gs = 0;
        int count_gf = 0;
        int thresh_s = 5;
        const int thresh_f = 500;
        int thresh_s5 = 5*5;
        const double contract_coeff = 0.5;
        const double expand_coeff = 1.2;
        const double search_radius = 0.1;
        const double rho_lower_bound = 1e-6;
        
        
        FILE* fout = fopen("bundle.txt", "wt");
        for (int iter=0; iter < 100000; iter++) {
            for (int r=0; r < p.rows(); r++) {
                vel[r] = dist(mt)*rho*scale[r];
            }
            vel += dir*inertia;
            pp = p + vel;
            double f = evaluate(pp);
            if (f < fbest) {
                fbest = f;
                dir = 0.9*dir + 0.1*vel;
                p = pp;
                count_gf = 0;
                count_gs++;
            } else {
                dir *= 0.9;
                count_gf++;
                count_gs = 0;
            }
            
            if (count_gs >= thresh_s) {
                rho *= expand_coeff;
                
                if (thresh_s5 < 5*100) {
                    thresh_s5++;
                }
                
            } else {
                if (count_gf >= thresh_f) {
                    rho *= contract_coeff;
                }
                
                if (thresh_s5 > 15*5) {
                    thresh_s5--;
                }
            }
            thresh_s = thresh_s5/5; // constant
  
            if (rho > search_radius) {
                rho = search_radius;  
            }
            if (rho < rho_lower_bound) {
                rho = rho_lower_bound;
            }
            
            fprintf(fout, "%lf %lf %lf\n", f, fbest, p[7]);
        }
        fclose(fout);
        
        printf("final f=%lf, distortion=%le, delta=(%lf, %lf) RMSE=%lf\n", 
            1.0/p[7], p[6], p[8], p[9], fbest
        );
        backproj(p);
        
        return p;
    }
    
    vector<Eigen::Vector2d>& img_points;
    vector<Eigen::Vector3d>& world_points;
    cv::Mat rot_mat;
    cv::Mat rod_angles;
    Eigen::VectorXd best_sol;
};

#endif
