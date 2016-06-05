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
    : chart_scale(1.0), largest_block_index(-1), focal_length(10000) {
    }
    
    void construct(Mtf_core& mtf_core, bool pose_based=false) {
    
        // check too see if our four main fiducials are present?
        int zcount = 0;
        for (auto e: mtf_core.ellipses) {
            if (e.code == 0) {
                zero.x += 0.5 * e.centroid_x;
                zero.y += 0.5 * e.centroid_y;
                zcount++;
            }
        }
        
        if (zcount == 2 && pose_based) { // TODO: based on above check
        
            map<int, vector<Ellipse_detector*> > by_code;
            Point2d first;
            Point2d last;
            zcount = 0;
            for (auto& e: mtf_core.ellipses) {
                by_code[e.code].push_back(&e);
                if (e.code == 0) {
                    if (zcount == 0) {
                        first = Point2d(e.centroid_x, e.centroid_y);
                    } else {
                        last = Point2d(e.centroid_x, e.centroid_y);
                    }
                    zcount++;
                }
            }
            
            printf("absolute centre: (%lf, %lf)\n", zero.x, zero.y);
            transverse = normalize(first - last);
            printf("transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
            // sign of transverse unknown at this stage
            longitudinal = Point2d(-transverse.y, transverse.x);
            printf("longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
            
            // we must find fiducials 2, 4, 6, and 8 to orient ourselves
            // code 8 is quadrant 0
            if (by_code.find(8) != by_code.end()) { // TODO: we can probably fall back on other codes here ...
                
                vector<int> to_locate{2,4,6,8};
                for (auto t: to_locate) {
                    for (int i=0; i < n_fiducials; i++) {
                        if (main_fiducials[i].code == t) {
                            fiducials[t] = vector<Fiducial>(4);
                            
                            fiducials[t][main_fiducials[i].quadrant].icoords = Point2d(by_code[t][0]->centroid_x, by_code[t][0]->centroid_y);
                            fiducials[t][main_fiducials[i].quadrant].rcoords = main_fiducials[i].rcoords;
                            
                            printf("Ellipse assigned to fiducial (code=%d, quad=%d) : (%lf %lf) -> (%lf %lf)\n",
                                t, main_fiducials[i].quadrant,
                                fiducials[t][main_fiducials[i].quadrant].icoords.x, fiducials[t][main_fiducials[i].quadrant].icoords.y,
                                fiducials[t][main_fiducials[i].quadrant].rcoords.x, fiducials[t][main_fiducials[i].quadrant].rcoords.y
                            );
                        }
                    }
                }
                
                Point2d dv(by_code[8][0]->centroid_x - zero.x, by_code[8][0]->centroid_y - zero.y);
                
                if ( (dv.x*transverse.x + dv.y*transverse.y) > 0 ) {
                    transverse = -transverse;
                    longitudinal = Point2d(-transverse.y, transverse.x);
                }
            
                printf("Oriented transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
                printf("Oriented longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
                
                
                for (auto& e: mtf_core.ellipses) {
                    Point2d pv(e.centroid_x - zero.x, e.centroid_y - zero.y);
                    double dot_l = longitudinal.x*pv.x + longitudinal.y*pv.y;
                    double dot_t = transverse.x*pv.x + transverse.y*pv.y;
                    
                    int quadrant = 0;
                    if (dot_l > 0 && dot_t <= 0) quadrant = 1;
                    if (dot_l > 0 && dot_t > 0) quadrant = 2;
                    if (dot_l <= 0 && dot_t <= 0) quadrant = 0;
                    if (dot_l <= 0 && dot_t > 0) quadrant = 3;
                    
                    // find the best match
                    for (int j=0; j < n_fiducials; j++) {
                        if (e.code == main_fiducials[j].code && quadrant == main_fiducials[j].quadrant) {
                            
                            auto it = fiducials.find(e.code);
                            if (it == fiducials.end()) {
                                fiducials[e.code] = vector<Fiducial>(4);
                            }
                            Fiducial& selected = fiducials[e.code][quadrant];
                            selected = Fiducial(
                                e.centroid_x,
                                e.centroid_y,
                                main_fiducials[j].rcoords.y, 
                                main_fiducials[j].rcoords.x,  
                                e.code,
                                quadrant
                            );
                            printf("Ellipse assigned to fiducial (code=%d, quad=%d) : (%lf %lf) -> (%lf %lf) (dot_l=%lf, dot_t=%lf)\n",
                                e.code, quadrant,
                                selected.icoords.x, selected.icoords.y,
                                selected.rcoords.x, selected.rcoords.y,
                                dot_l, dot_t
                            );
                        }
                    }
                }
                
                prin = Point2d(mtf_core.img.cols/2.0, mtf_core.img.rows/2.0);
                
                // take 2
                vector<cv::Point2i> initial{ {6,0}, {4,3}, {2,2}, {8,1}, {10,2} };
                
                vector<Eigen::Vector2d> feature_points(5);
                vector<Eigen::Vector3d> world_points(5);
                int fidcount = 0;
                Point2d pix_centroid(0, 0);
                int totalcount = 0;
                for (auto p: initial) {
                    Fiducial& fid = fiducials[p.x][p.y];
                    feature_points[fidcount] = Eigen::Vector2d(
                        (fid.icoords.x - prin.x),
                        (fid.icoords.y - prin.y)
                    );
                    world_points[fidcount] = Eigen::Vector3d(fid.rcoords.x, fid.rcoords.y, 1 + 0.000001*(fidcount+1));
                    fidcount++;
                    if (p.x <= 8) {    
                        pix_centroid.x += fid.icoords.x;
                        pix_centroid.y += fid.icoords.y;
                        totalcount++;
                    }
                }
                pix_centroid *= 1.0/double(totalcount);
                printf("pix centroid =[%lf %lf]\n", pix_centroid.x, pix_centroid.y);
                
                img_scale = std::max(mtf_core.img.rows, mtf_core.img.cols);
                
                vector<Eigen::Vector2d> ba_img_points;
                vector<Eigen::Vector3d> ba_world_points;
                vector<int> perm;
                for (auto fidv: fiducials) {
                    for (auto fid: fidv.second) {
                        if (fid.icoords.x == 0 && fid.icoords.y == 0) continue;
                        ba_img_points.push_back(Eigen::Vector2d((fid.icoords.x - prin.x)/img_scale, (fid.icoords.y - prin.y)/img_scale));
                        ba_world_points.push_back(Eigen::Vector3d(fid.rcoords.x, fid.rcoords.y, 1.0));
                        perm.push_back(perm.size());
                    }
                }
                
                
                vector<Eigen::Matrix<double, 3, 4> > projection_matrices;
                vector<vector<double> > radial_distortions;
                cv::Mat rot_matrix = cv::Mat(3, 3, CV_64FC1);
                cv::Mat rod_angles = cv::Mat(3, 1, CV_64FC1);
                Eigen::MatrixXd P;
                
                std::mt19937 mt(101);
                std::uniform_real_distribution<double> dist(0, 1);
                        
                
                double global_bpr = 1e50;
                for (int ri=0; ri < 20000; ri++) {
                
                    // shuffle perm
                    for (size_t yi=perm.size()-1; yi > 1; yi--) {
                        int j = (int)floor(dist(mt)*yi);
                        std::swap(perm[yi], perm[j]);
                    }
                
                    for (int i=0; i < 5; i++) {
                        feature_points[i] = ba_img_points[perm[i]];
                        world_points[i] = ba_world_points[perm[i]];
                        world_points[i][2] += + 0.000001*(i+1);
                    }
                    
                    theia::FivePointFocalLengthRadialDistortion(
                        feature_points,
                        world_points,
                        1, // number of distortion parms
                        &projection_matrices,
                        &radial_distortions
                    );
                    
                    double min_bpr = 1e50;
                    for (size_t k=0; k < projection_matrices.size(); k++) {
                    
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
                        
                            double bpr = 0;
                            for (size_t i=0; i < ba_img_points.size(); i++) {
                                Eigen::VectorXd rp = projection_matrices[k]*Eigen::Vector4d(ba_world_points[i][0], ba_world_points[i][1], ba_world_points[i][2], 1.0);
                                rp /= rp[2];
                                double rad = 1 + radial_distortions[k][0]*(rp[0]*rp[0] + rp[1]*rp[1]);
                                rp /= rad;
                                double err = (ba_img_points[i] - Eigen::Vector2d(rp[0], rp[1])).norm();
                                bpr += fabs(err); // more robustness against outliers is appreciated here
                            }
                            bpr /= ba_img_points.size();
                        
                            
                            
                            // only keep a rotation matrix if it repeats with the Rodrigues angle convention
                            if (bpr <= min_bpr && rot_err < 0.01) { 
                                min_bpr = bpr;
                                
                                if (bpr < global_bpr) {
                                    global_bpr = bpr;
                                    P = projection_matrices[k];
                                    distortion = radial_distortions[k][0];
                                    printf("%lu[%d]: rotation error: %lf, bpr=%lf pixels, f=%lf pixels\n", k, ri, rot_err, bpr*img_scale, img_scale/w);
                                }
                            }
                            
                        }
                    }
                    projection_matrices.clear();
                    radial_distortions.clear();
                }
                
                double r1n = (P.block(0,0,1,3)).norm();
                double r3n = (P.block(2,0,1,3)).norm();
                double w = sqrt( (r3n*r3n)/(r1n*r1n) );
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
                
                
                Point2d delta_principal(0,0);
                Bundle_adjuster ba(
                    ba_img_points, ba_world_points, 
                    Pcop, 
                    rod_angles,
                    distortion,
                    w,
                    delta_principal.x, delta_principal.y,
                    img_scale
                );
                
                ba.unpack(rotation, translation, distortion, w, delta_principal.x, delta_principal.y);
                focal_length = 1.0/w;
                
                prin.x -= delta_principal.x;
                prin.y -= delta_principal.y;
                
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
                
                centre_depth = backproject(pix_centroid.x, pix_centroid.y)[2];
            }
            
            // construct a distance scale
            
        } else {
            printf("Warning: Manual focus profile chart mode requested, but central dots not found\n");
            const vector<Block>& blocks = mtf_core.get_blocks();
            // find largest block
            vector< pair<double, int> > by_size;
            for (size_t i=0; i < blocks.size(); i++) {
                by_size.push_back(make_pair(blocks[i].get_area(), i));
            }
            sort(by_size.begin(), by_size.end());
            
            const int sstart = 5;
            double delta_1 = by_size[by_size.size()-1].first / by_size[by_size.size()-sstart-1].first;
            double delta_2 = by_size[by_size.size()-sstart-1].first / by_size[by_size.size()-sstart-2].first;
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
    
    bool bounds(double& lmin, double &lmax, double& step) {
        if (distance_scale.size() > 2) {
            lmin = distance_scale.front().x;
            lmax = distance_scale.back().x;
            step = fabs(distance_scale.back().x - distance_scale[distance_scale.size()-2].x);
            return true;
        }
        return false;
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
    
    Point2d world_to_image(double world_x, double world_y) const {
        Eigen::Matrix3d R = rotation;
    
        R.row(2) *= 1.0/focal_length; 
        Eigen::Vector3d t(translation[0], translation[1], translation[2]/focal_length);
        
        Eigen::Vector3d bp = fwdP*Eigen::Vector3d(world_x, world_y, 1.0) + fwdT;
            
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
    vector<cv::Point3d> distance_scale;
    int largest_block_index;
    
    map<int, vector<Fiducial> > fiducials;
    
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
    
  private:
    template <typename T> int sgn(T val) const {
        return (T(0) < val) - (val < T(0));
    }
};

#endif
