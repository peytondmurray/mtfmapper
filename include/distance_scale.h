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
#ifndef DISTANCE_SCALE_H
#define DISTANCE_SCALE_H

#include "mtf_core.h"
#include "ellipse_decoder.h"

#include "Eigen/SVD"

#include "fiducial_positions.h"
#include "five_point_focal_length_radial_distortion.h"

#include "bundle.h"

static inline bool t_intersect(double& pix, double& piy, 
                 const double& v1x, const double& v1y,
                 const double& d1x, const double& d1y,
                 const double& v2x, const double& v2y,
                 const double& d2x, const double& d2y)  {
               
    double denom = (d2y*d1x - d2x*d1y);
    
    if (fabs(denom) < 1e-11) {
       // this happens when the lines are parallel
       // the caller handles this correctly by not
       // adding an additional intersection point
       return false;
    }
    
    double u = (d2x*(v1y - v2y) - d2y*(v1x - v2x)) / denom;
    
    pix = v1x + u*d1x;
    piy = v1y + u*d1y;
               
    return true;               
}

class Distance_scale {
  public:
    Distance_scale(void)
    : chart_scale(1.0), largest_block_index(-1), focal_length(10000), fiducials_found(false) {
    }
    
    void construct(Mtf_core& mtf_core, bool pose_based=false) {
    
        int zcount = 0;
        
        for (size_t i=0; i < mtf_core.ellipses.size() - 1; i++) {
            Ellipse_detector& e = mtf_core.ellipses[i];
            for (size_t j=i+1; j < mtf_core.ellipses.size(); j++) {
                Ellipse_detector& f = mtf_core.ellipses[j];
                double cdist = sqrt(SQR(f.centroid_x - e.centroid_x) + SQR(f.centroid_y - e.centroid_y));
                if (cdist < std::min(e.minor_axis, f.minor_axis)) {
                    // keep only the larger fiducial
                    if (e.major_axis > f.major_axis) {
                        f.valid = false;
                    } else {
                        e.valid = false;
                    }
                }
            }
        }
        
        for (const auto& e: mtf_core.ellipses) {
            if (e.valid && e.code == 0) {
                zero.x += 0.5 * e.centroid_x;
                zero.y += 0.5 * e.centroid_y;
                zcount++;
            }
        }
        
        if (zcount == 2 && pose_based) {
        
            map<int, Ellipse_detector* > by_code;
            Point2d first;
            Point2d last;
            zcount = 0;
            double max_fiducial_diameter = 0;
            for (auto& e: mtf_core.ellipses) {
                if (!e.valid) continue;
                if (e.code == 0) {
                    if (zcount == 0) {
                        first = Point2d(e.centroid_x, e.centroid_y);
                    } else {
                        last = Point2d(e.centroid_x, e.centroid_y);
                    }
                    zcount++;
                } else {
                    if (by_code.find(e.code) != by_code.end()) {
                        printf("collision: code %d found at both [%lf %lf] and [%lf %lf]\n",
                            e.code, e.centroid_x, e.centroid_y,
                            by_code[e.code]->centroid_x, by_code[e.code]->centroid_y
                        );
                    }
                    by_code[e.code] = &e;
                }
                max_fiducial_diameter = std::max(e.major_axis, max_fiducial_diameter);
            }
            
            printf("absolute centre: (%lf, %lf)\n", zero.x, zero.y);
            transverse = normalize(first - last);
            printf("transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
            // sign of transverse unknown at this stage
            longitudinal = Point2d(-transverse.y, transverse.x);
            printf("longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
            
            if (by_code.size() > 6) { // minimum of 5 unique fiducials plus one spare plus the zeros
                
                // use knowledge of the first four fiducials to orient the transverse direction
                vector< pair<int, double> > candidates {{1, 1.0}, {2, 1.0}, {3, -1.0}, {4, -1.0}};
                for (auto c: candidates) {
                    if (by_code.find(c.first) != by_code.end()) {
                        Point2d dir(by_code[c.first]->centroid_x -  zero.x, by_code[c.first]->centroid_y - zero.y);
                        dir *= c.second;
                        if (dir.x * transverse.x + dir.y * transverse.y < 0) {
                            transverse = -transverse;
                            longitudinal = Point2d(-transverse.y, transverse.x);
                        }
                        break; // one sample is enough
                    }    
                }
                
                
                prin = Point2d(mtf_core.img.cols/2.0, mtf_core.img.rows/2.0);
                img_scale = std::max(mtf_core.img.rows, mtf_core.img.cols);
                
                vector<Eigen::Vector2d> ba_img_points;
                vector<Eigen::Vector3d> ba_world_points;
                for (const auto& e: mtf_core.ellipses) {
                    if (!e.valid) continue;
                    if (e.code == 0) continue; // we do not know how to tell the two zeros appart, so just skip them
                    int main_idx = -1;
                    for (size_t i=0; i < n_fiducials && main_idx == -1; i++) {
                        if (main_fiducials[i].code == e.code) {
                            main_idx = i;
                        }
                    }
                    
                    ba_img_points.push_back(Eigen::Vector2d((e.centroid_x - prin.x)/img_scale, (e.centroid_y - prin.y)/img_scale));
                    ba_world_points.push_back(
                        Eigen::Vector3d(
                            main_fiducials[main_idx].rcoords.y, 
                            main_fiducials[main_idx].rcoords.x, 
                            1.0
                        )
                    );
                }
                
                vector<Eigen::Matrix<double, 3, 4> > projection_matrices;
                vector<vector<double> > radial_distortions;
                cv::Mat rot_matrix = cv::Mat(3, 3, CV_64FC1);
                cv::Mat rod_angles = cv::Mat(3, 1, CV_64FC1);
                Eigen::MatrixXd P;
                
                class Cal_solution {
                  public:
                    Cal_solution(double bpe=0, Eigen::MatrixXd proj=Eigen::MatrixXd(), double distort=0, vector<int> inlier_list=vector<int>())
                     : bpe(bpe), proj(proj), distort(distort), inlier_list(inlier_list) {};
                     
                    bool operator< (const Cal_solution& b) const {
                        return bpe < b.bpe;
                    }
                                 
                    double bpe;
                    Eigen::MatrixXd proj;
                    double distort;
                    vector<int> inlier_list;
                };
                
                vector<Cal_solution> solutions;
                vector<Eigen::Vector2d> feature_points(5);
                vector<Eigen::Vector3d> world_points(5);
                
                enumerate_combinations(ba_img_points.size(), 5);
                
                double inlier_threshold = max_fiducial_diameter; // this should work unless extreme distortion is present?
                double global_bpr = 1e50;
                for (size_t ri=0; ri < combinations.size(); ri++) {
                
                    for (int i=0; i < 5; i++) {
                        feature_points[i] = ba_img_points[combinations[ri][i]];
                        world_points[i] = ba_world_points[combinations[ri][i]];
                        world_points[i][2] += + 0.000001*(i+1);
                    }
                    
                    theia::FivePointFocalLengthRadialDistortion(
                        feature_points,
                        world_points,
                        1, // number of distortion parms
                        &projection_matrices,
                        &radial_distortions
                    );
                    
                    for (size_t k=0; k < projection_matrices.size(); k++) {
                    
                        // only consider solutions in front of camera
                        if (projection_matrices[k](2, 3) > 0) {
                            projection_matrices[k] *= -1;
                        }
                    
                        // check if Rodriques produces the same rotation matrix
                        double r1n = (projection_matrices[k].block(0,0,1,3)).norm();
                        double r3n = (projection_matrices[k].block(2,0,1,3)).norm();
                        double w = sqrt( (r3n*r3n)/(r1n*r1n) );
                        
                        Eigen::MatrixXd RM = projection_matrices[k].block(0,0,3,3);
                        RM.row(2) /= w;
                        RM /= RM.row(0).norm();
                        
                        for (int rr=0; rr < 3; rr++) {
                            for (int cc=0; cc < 3; cc++) {
                                rot_matrix.at<double>(rr,cc) = RM(rr,cc);
                            }
                        }
                        cv::Rodrigues(rot_matrix, rod_angles);
                        cv::Rodrigues(rod_angles, rot_matrix);
                        double rot_err = 0;
                        for (int rr=0; rr < 3; rr++) {
                            for (int cc=0; cc < 3; cc++) {
                                rot_err += fabs(rot_matrix.at<double>(rr,cc) - RM(rr,cc));
                            }
                        }
                        
                        if (rot_err < 0.01) { 
                            vector<double> residuals;
                            
                            // TODO: technically, we could re-run the five-point solver with eccentricity correction
                            // but since this should have very little effect, rather leave that for the
                            // final bundle adjustment
                            
                            Eigen::Matrix3d RMM(RM);
                            RMM.row(2) *= w;
                            Eigen::VectorXd TV = projection_matrices[k].col(3) / projection_matrices[k].block(0,0,1,3).norm();
                            
                            double r = 5.0;   // radius of circle
                            double c = 1.0/w; // focal length, but in which units?
                            Eigen::Vector3d ii = RMM.col(2); // already includes focal length scaling ... ?
                            
                            double max_err = 0;
                            vector<int> inliers;
                            for (size_t i=0; i < ba_img_points.size(); i++) {
                                
                                Eigen::Vector3d v = TV - RMM*ba_world_points[i];
                                
                                Eigen::VectorXd a(5);
                                a[0] = r*r*ii[0]*ii[0] - SQR(v[1]*ii[1] + v[2]*ii[2]) - ii[0]*ii[0]*(v[1]*v[1] + v[2]*v[2]);
                                a[1] = 2*(r*r*ii[0]*ii[1] - (v[2]*ii[0] - v[0]*ii[2])*(v[2]*ii[1] - v[1]*ii[2]) + v[0]*v[1]);
                                a[2] = r*r*ii[1]*ii[1] + SQR(v[0]*ii[2] - v[2]*ii[0]) - v[0]*v[0] - v[2]*v[2];
                                a[3] = 2*c*(v[0]*v[2]*(ii[1]*ii[1] - 1) + v[1]*(-v[0]*ii[1]*ii[2] + v[1]*ii[0]*ii[2] - v[2]*ii[0]*ii[1]) - r*r*ii[0]*ii[2]);
                                a[4] = 2*c*((v[0]*ii[2] - v[2]*ii[0])*(v[0]*ii[1] - v[1]*ii[0]) - v[1]*v[2] - r*r*ii[1]*ii[2]);
                                
                                double d = c*c*(-SQR(v[1]*ii[0] - v[0]*ii[1]) + v[0]*v[0] + v[1]*v[1] - r*r*ii[2]*ii[2]);
                                a /= d;
                                
                                double uc = (a[1]*a[4] - 2*a[2]*a[3]) / (4*a[0]*a[2] - a[1]*a[1]);
                                double vc = (a[1]*a[3] - 2*a[0]*a[4]) / (4*a[0]*a[2] - a[1]*a[1]);
                                
                                Eigen::Vector3d bp = RMM*Eigen::Vector3d(ba_world_points[i][0], ba_world_points[i][1], 1.0) + TV;
                                bp /= bp[2];
                                
                                double rad = 1 + radial_distortions[k][0]*(bp[0]*bp[0] + bp[1]*bp[1]);
                                bp /= rad;
                                
                                // correct for eccentricity error
                                Eigen::Vector2d corrected_img_point = ba_img_points[i] + Eigen::Vector2d(bp[0] - uc, bp[1] - vc)/img_scale;
                                
                                double err = (corrected_img_point - Eigen::Vector2d(bp[0], bp[1])).norm();
                                
                                if (err*img_scale < inlier_threshold) {
                                    residuals.push_back(err);
                                    inliers.push_back(i);
                                    max_err = std::max(err, max_err);
                                } 
                            }
                            
                            if (inliers.size() < 4) continue; // do not waste time on outlier-dominated solutions
                            
                            sort(residuals.begin(), residuals.end());
                            double bpr = 0;
                            // try to exclude the worst two or so specimens
                            size_t nresiduals = std::max(size_t(0.9*residuals.size()), std::min(residuals.size(), size_t(5)));
                            for (size_t i=0; i < nresiduals; i++) {
                                bpr += residuals[i]*residuals[i];
                            }
                            bpr = sqrt(bpr/double(nresiduals))*img_scale;
                            
                            if (solutions.size() < 5) { // keep 5 best solutions
                                solutions.push_back(Cal_solution(bpr, projection_matrices[k], -radial_distortions[k][0], inliers));
                            } else {
                                auto worst_sol = solutions.rbegin();
                                if (bpr < worst_sol->bpe) {
                                    *worst_sol = Cal_solution(bpr, projection_matrices[k], -radial_distortions[k][0], inliers);
                                    sort(solutions.begin(), solutions.end());
                                }
                            }
                            
                            if (bpr < global_bpr) {
                                global_bpr = bpr;
                                P = projection_matrices[k];
                                distortion = radial_distortions[k][0];
                                printf("%lu[%d]: rotation error: %lf, bpr=%lf pixels, f=%lf pixels, inliers=%lu\n", 
                                    k, ri, rot_err, bpr, img_scale/w, inliers.size()
                                );
                            }
                        }
                    }
                    projection_matrices.clear();
                    radial_distortions.clear();
                }
                combinations.clear(); // discard the combinations
                
                if (solutions.size() == 0) {
                    printf("Error: Solutions to camera calibration found. Aborting\n");
                    fiducials_found = false;
                    return;
                }
                
                printf("Retained %lu promising solutions\n", solutions.size());
                for (auto s: solutions) {
                    printf("\t solution with bpe = %lf, dist=%le\n", s.bpe, s.distort);
                }
                double w = 0;
                
                int min_idx = 0;
                
                try_next_solution:
                
                P = solutions[min_idx].proj;
                distortion = solutions[min_idx].distort;
                vector<int>& inliers = solutions[min_idx].inlier_list;
                
                double r1n = (P.block(0,0,1,3)).norm();
                double r3n = (P.block(2,0,1,3)).norm();
                w = sqrt( (r3n*r3n)/(r1n*r1n) );
                focal_length = 1.0/w;
                
                P /= P.block(0,0,1,3).norm(); // remove the arbitrary scaling factor
                P.row(2) /= w; // remove focal length from last row
                
                Eigen::MatrixXd RM = P.block(0,0,3,3);
                Eigen::Vector3d Pcop = P.col(3);
                
                std::cout << "R=\n" << RM << std::endl;
                std::cout << "T= " << Pcop.transpose() << std::endl;
                printf("focal length = %lf (times max sensor dim), or %lf pixels\n", 1.0/w, img_scale/w);
                
                for (int rr=0; rr < 3; rr++) {
                    for (int cc=0; cc < 3; cc++) {
                        rot_matrix.at<double>(rr,cc) = RM(rr,cc);
                    }
                }
                cv::Rodrigues(rot_matrix, rod_angles);
                
                
                vector<Eigen::Vector2d> inlier_feature_points(inliers.size());
                vector<Eigen::Vector3d> inlier_world_points(inliers.size());
                
                // prepate rotation matrix and 
                Eigen::Matrix3d RMM(RM);
                RMM.row(2) *= w;
                Eigen::VectorXd TV = Pcop;
                TV[2] *= w;
                
                for (size_t k=0; k < inliers.size(); k++) {
                
                    double r = 5.0;   // radius of circle
                    double c = 1.0/w; // focal length, but in which units?
                    Eigen::Vector3d i = RMM.col(2); // already includes focal length scaling ... ?
                    Eigen::Vector3d v = TV - RMM*ba_world_points[inliers[k]];
                    
                    Eigen::VectorXd a(5);
                    a[0] = r*r*i[0]*i[0] - SQR(v[1]*i[1] + v[2]*i[2]) - i[0]*i[0]*(v[1]*v[1] + v[2]*v[2]);
                    a[1] = 2*(r*r*i[0]*i[1] - (v[2]*i[0] - v[0]*i[2])*(v[2]*i[1] - v[1]*i[2]) + v[0]*v[1]);
                    a[2] = r*r*i[1]*i[1] + SQR(v[0]*i[2] - v[2]*i[0]) - v[0]*v[0] - v[2]*v[2];
                    a[3] = 2*c*(v[0]*v[2]*(i[1]*i[1] - 1) + v[1]*(-v[0]*i[1]*i[2] + v[1]*i[0]*i[2] - v[2]*i[0]*i[1]) - r*r*i[0]*i[2]);
                    a[4] = 2*c*((v[0]*i[2] - v[2]*i[0])*(v[0]*i[1] - v[1]*i[0]) - v[1]*v[2] - r*r*i[1]*i[2]);
                    
                    double d = c*c*(-SQR(v[1]*i[0] - v[0]*i[1]) + v[0]*v[0] + v[1]*v[1] - r*r*i[2]*i[2]);
                    a /= d;
                    
                    double uc = (a[1]*a[4] - 2*a[2]*a[3]) / (4*a[0]*a[2] - a[1]*a[1]);
                    double vc = (a[1]*a[3] - 2*a[0]*a[4]) / (4*a[0]*a[2] - a[1]*a[1]);
                    
                    Eigen::Vector3d bp = RMM*Eigen::Vector3d(ba_world_points[inliers[k]][0], ba_world_points[inliers[k]][1], 1.0) + TV;
                    bp /= bp[2];
                    
                    printf("point %d: [%lf %lf] -> ellipse centre [%lf %lf], point centre [%lf %lf]\n",
                        inliers[k], img_scale*bp[0] + prin.x, img_scale*bp[1] + prin.y,
                        uc*img_scale + prin.x, vc*img_scale + prin.y, 
                        (bp[0] - uc), (bp[1] - vc)
                    );
                    
                    inlier_feature_points[k] = ba_img_points[inliers[k]] + Eigen::Vector2d(bp[0] - uc, bp[1] - vc)/img_scale;
                    inlier_world_points[k] = ba_world_points[inliers[k]];
                }
                
                Bundle_adjuster ba(
                    inlier_feature_points, inlier_world_points, 
                    Pcop, 
                    rod_angles,
                    distortion,
                    w,
                    img_scale
                );
                ba.solve();
                ba.unpack(rotation, translation, distortion, w);
                bundle_rmse = ba.evaluate(ba.best_sol)*img_scale;
                printf("solution %d has rmse=%lf\n", min_idx, bundle_rmse);
                
                if (ba.optimization_failure() && min_idx < int(solutions.size() - 1)) {
                    min_idx++;
                    goto try_next_solution;
                }
                
                // take another stab, in case NM stopped at a point that
                // was not really a minimum
                // TODO: in theory, we could repeat this until RMSE stabilizes?
                ba.solve();
                ba.unpack(rotation, translation, distortion, w);
                bundle_rmse = ba.evaluate(ba.best_sol)*img_scale;
                printf("2nd solution %d has rmse=%lf\n", min_idx, bundle_rmse);
                focal_length = 1.0/w;
                
                // prepare for backprojection
                Eigen::Matrix3d K;
                K << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1.0 / focal_length;
                
                
                invP = (K*rotation).inverse();
                
                // prepare for forward and backward projection transforms    
                fwdP = K*rotation;
                fwdT = Eigen::Vector3d(translation[0], translation[1], translation[2]/focal_length);
                
                cop = Eigen::Vector3d(translation[0], translation[1], translation[2]/focal_length);
                cop = - invP * cop;
                
                centre_depth = backproject(zero.x, zero.y)[2];
                
                printf("ultimate f=%lf centre_depth=%lf distortion=%le\n", focal_length*img_scale, centre_depth, distortion);
                fiducials_found = true;
            }
            
            // construct a distance scale
            
        } else {
            printf("Warning: Manual focus profile chart mode requested, but central dots not found\n");
            const vector<Block>& blocks = mtf_core.get_blocks();
            
            double delta_1 = 1;
            double delta_2 = 1;
            vector< pair<double, int> > by_size;
            
            if (blocks.size() > 0) {
                // find largest block
                for (size_t i=0; i < blocks.size(); i++) {
                    by_size.push_back(make_pair(blocks[i].get_area(), i));
                }
                sort(by_size.begin(), by_size.end());
                
                const int sstart = 5;
                delta_1 = by_size[by_size.size()-1].first / by_size[by_size.size()-sstart-1].first;
                delta_2 = by_size[by_size.size()-sstart-1].first / by_size[by_size.size()-sstart-2].first;    
            }
            
            if (delta_1/delta_2 > 80) {
                largest_block_index = by_size.back().second;
                const Block& lblock = blocks[by_size.back().second];
                // we have a clear largest block. now determine its orientation
                vector<double> xp;
                vector<double> yp;
                for (size_t i=0; i < blocks.size(); i++) {
                    xp.push_back(blocks[i].centroid.x);
                    yp.push_back(blocks[i].centroid.y);
                }
                sort(xp.begin(), xp.end());
                sort(yp.begin(), yp.end());
                double idx_x = (find(xp.begin(), xp.end(), lblock.centroid.x) - xp.begin())/double(xp.size());
                double idx_y = (find(yp.begin(), yp.end(), lblock.centroid.y) - yp.begin())/double(yp.size());
                printf("xfrac=%lf, yfrac=%lf\n", idx_x, idx_y);
                
                Point2d median(xp[xp.size()/2], yp[yp.size()/2]); 
                
                // orientations relative to chart, not image
                int top_i=lblock.get_edge_index(Block::TOP);
                int bot_i=lblock.get_edge_index(Block::BOTTOM);
                int left_i=lblock.get_edge_index(Block::LEFT);
                int right_i=lblock.get_edge_index(Block::RIGHT);
                if (fabs(idx_x - 0.5) < fabs(idx_y - 0.5)) {
                    // outer rows arranged in columns
                    if (fabs(lblock.get_edge_centroid(top_i).y - median.y) <
                        fabs(lblock.get_edge_centroid(bot_i).y - median.y) ) {
                        
                        // chart is upside down
                        std::swap(top_i, bot_i);
                        std::swap(left_i, right_i);
                    }
                } else {
                    // outer rows arranged in rows
                    std::swap(bot_i, right_i);
                    std::swap(top_i, left_i);
                    
                    if (fabs(lblock.get_edge_centroid(top_i).x - median.x) <
                        fabs(lblock.get_edge_centroid(bot_i).x - median.x) ) {
                        
                        // chart is upside down
                        std::swap(top_i, bot_i);
                        std::swap(left_i, right_i);
                    }
                }
                transverse = normalize(lblock.get_edge_centroid(right_i) - lblock.get_edge_centroid(left_i));
                longitudinal = normalize(lblock.get_edge_centroid(bot_i) - lblock.get_edge_centroid(top_i));
                zero = lblock.get_edge_centroid(bot_i);
                double block_width = norm(lblock.get_edge_centroid(right_i) - lblock.get_edge_centroid(left_i));
                chart_scale = 62.0 / block_width; // assume central block is 62 mm wide
                printf("Warning: choosing (potentially) poor chart scale of %lf mm/pixel\n", chart_scale);
            } else { 
                printf("Warning: Could not identify largest block, choosing poor defaults\n");
                chart_scale = 0.15;
                zero = Point2d(932, 710);
                transverse = Point2d(1,0);
                longitudinal = Point2d(0,1);
            }
            printf("zero: %lf %lf\n", zero.x, zero.y);
            printf("transverse: %lf %lf\nlongitudinal: %lf %lf\n", transverse.x, transverse.y, longitudinal.x, longitudinal.y);
            
        }
    }
    
    
    inline Point2d normalize_img_coords(double pixel_x, double pixel_y) const {
        double xs = (pixel_x - prin.x)/img_scale;
        double ys = (pixel_y - prin.y)/img_scale;
        
        double rhat = sqrt(xs*xs + ys*ys);
        double r = 0;
        
        if (rhat > 1e-10) {
            double pa=distortion*rhat;
            double pb=-1;
            double pc=rhat;
            
            double q = -0.5 * (pb + sgn(pb)*sqrt(pb*pb - 4*pa*pc) );
            double r1 = q/pa;
            double r2 = pc / q;
            
            if (r1 > 0 && r2 > 0) {
                r = std::min(r1, r2);
            } else {
                if (r1 <= 0) {
                    r = r2;
                } else {
                    r = r1;
                }
            }
            
            xs = xs * rhat / r;
            ys = ys * rhat / r;
        }
        return Point2d(xs, ys);
    }
    
    Eigen::Vector3d backproject(double pixel_x, double pixel_y) const {
        Point2d ic = normalize_img_coords(pixel_x, pixel_y);
        
        Eigen::Vector3d dv(ic.x, ic.y, 1.0);
        dv = invP*dv;
        dv /= dv.norm();
        
        double s = (1 - cop[2])/dv[2];
        Eigen::Vector3d ip = cop + s*dv;
        
        // now we have ip in world coordinates, but we actually want it in camera 
        // coordinates
        
        // z-axis sign?
        ip = rotation*ip + Eigen::Vector3d(translation[0], translation[1], -translation[2]);
        
        return ip;
    }
    
    void estimate_depth_img_coords(double pixel_x, double pixel_y, double& depth) const {
        Eigen::Vector3d bp = backproject(pixel_x, pixel_y);
        depth = bp[2] - centre_depth;
    }
    
    void estimate_depth_world_coords(double world_x, double world_y, double& depth) const {
        Eigen::Vector3d bp = rotation*Eigen::Vector3d(world_x, world_y, 0) + 
            Eigen::Vector3d(translation[0], translation[1], -translation[2]);
        
        depth = bp[2] - centre_depth;
    }
    
    Point2d estimate_world_coords(double pixel_x, double pixel_y) const {
        Point2d ic = normalize_img_coords(pixel_x, pixel_y);
        
        Eigen::Vector3d dv(ic.x, ic.y, 1.0);
        dv = invP*dv;
        dv /= dv.norm();
        
        double s = (1 - cop[2])/dv[2];
        Eigen::Vector3d ip = cop + s*dv;
        
        return Point2d(ip[0], ip[1]); // ip[2] (=z) will always be in the plane
    }
    
    Point2d world_to_image(double world_x, double world_y, double world_z=1.0) const {
        Eigen::Vector3d bp = fwdP*Eigen::Vector3d(world_x, world_y, world_z) + fwdT;
            
        bp /= bp[2];
            
        double rad = 1 + distortion*(bp[0]*bp[0] + bp[1]*bp[1]);
        bp /= rad;
        bp *= img_scale;
        
        bp[0] += prin.x;
        bp[1] += prin.y;
    
        return Point2d(bp[0], bp[1]); 
    }
    
    
    Point2d zero;
    Point2d transverse;
    Point2d longitudinal;
    double chart_scale;
    int largest_block_index;
    
    Point2d prin;
    double focal_length;
    double centre_depth;
    double distortion;
    double img_scale;
    Eigen::Matrix3d rotation;
    Eigen::Vector3d translation;
    
    Eigen::Matrix3d fwdP;
    Eigen::Vector3d fwdT;
    Eigen::Matrix3d invP; // from image to world coords
    Eigen::Vector3d cop;  // centre of projection (i.e., camera)
    
    double bundle_rmse;
    bool fiducials_found;
    
  private:
    template <typename T> int sgn(T val) const {
        return (T(0) < val) - (val < T(0));
    }
    
    vector< vector<int> > combinations;
    
    int n_choose_k(int n, int k) {
        if (k == n) { 
            return 1; 
        }
        if (k == 1) { 
            return n; 
        }
        return n_choose_k(n - 1, k) + n_choose_k(n - 1, k - 1);
    }

    void sub_enumerate_n_choose_k(int n, int k, vector< vector<int> >& all, int& j, vector<int>& a, int i) {
        a[i] = n - 1;
        if (i == k - 1) {
            all[j] = a;
            j++;
            return;
        }
        for (int c=n - 1; c > 0; c--) {
            sub_enumerate_n_choose_k(c, k, all, j, a, i + 1);
        }
    }

    void enumerate_n_choose_k(int n, int k, vector< vector<int> >& arr) {
        int j = 0;
        vector<int> a(k, 0);
        for (int c=n; c >= k; c--) {
            sub_enumerate_n_choose_k(c, k, arr, j, a, 0);
        }
    }
    
    void enumerate_combinations(int n, int k) {
        int t = n_choose_k(n, k);
        
        combinations = vector< vector<int> > (t, vector<int>(k, 0));
        
        enumerate_n_choose_k(n, k, combinations);
    }
};

#endif
