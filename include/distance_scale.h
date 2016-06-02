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
    
        // define reference frame from three solid circular dots in centre of chart
        printf("#solid ellipses: %d\n", (int)mtf_core.solid_ellipses.size());
        for (size_t j=0; j < mtf_core.solid_ellipses.size(); j++) {
            printf("(%lf %lf)\n", mtf_core.solid_ellipses[j].x, mtf_core.solid_ellipses[j].y);
        }
        printf("\n\n");
        if (mtf_core.solid_ellipses.size() >= 3) {
            printf("got all three solid ellipses\n");
            double maxdelta = 0;
            Point2d first;
            Point2d last;
            for (int i=0; i <= 1; i++) {
                for (int j=i+1; j < 3; j++) {
                    Point2d d = mtf_core.solid_ellipses[i] - mtf_core.solid_ellipses[j];
                    double dist = sqrt(d.x*d.x + d.y*d.y);
                    printf("considering (%lf, %lf) vs (%lf, %lf), dist=%lf ",
                        mtf_core.solid_ellipses[i].x, mtf_core.solid_ellipses[i].y,
                        mtf_core.solid_ellipses[j].x, mtf_core.solid_ellipses[j].y,
                        dist
                    );
                    if (dist > maxdelta) {
                        first = mtf_core.solid_ellipses[i];
                        last  = mtf_core.solid_ellipses[j];
                        maxdelta = dist;
                        printf("<- keeping this one\n");
                    } else printf("\n");
                }
            }
            Point2d base = first - last;
            chart_scale = 91.0 / sqrt(base.x*base.x + base.y*base.y);
            zero = Point2d(0.5*(first.x+last.x), 0.5*(first.y+last.y));
            
            printf("absolute centre: (%lf, %lf)\n", zero.x, zero.y);
            transverse = normalize(first - last);
            printf("transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
            // sign of transverse unknown at this stage
            longitudinal = Point2d(-transverse.y, transverse.x);
            printf("longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
            printf("chart scale is %lf mm per pixel\n", chart_scale);
            
            // now try to decode the ellipses
            vector< pair<double, int> > dist_by_idx;
            for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                Ellipse_decoder ed(mtf_core.ellipses[i], mtf_core.img, transverse);
                mtf_core.ellipses[i].set_code(ed.code);
                if (!mtf_core.ellipses[i].solid) {
                    Point2d d(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
                    dist_by_idx.push_back(make_pair(sqrt(d.x*d.x + d.y*d.y), i));
                }
            }
            sort(dist_by_idx.begin(), dist_by_idx.end());
            int up = dist_by_idx[0].second;
            int down = dist_by_idx[1].second;
            if (mtf_core.ellipses[up].code > 3 || mtf_core.ellipses[down].code > 3) {
                for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                    mtf_core.ellipses[i].set_code(Ellipse_decoder::reverse(mtf_core.ellipses[i].code));
                }
                transverse = -transverse;
                printf("flipping transverse direction\n");
                std::swap(first, last);
            }
            if (mtf_core.ellipses[up].code > 3 || mtf_core.ellipses[down].code > 3) {
                printf("Warning: Ellipse targets closest to zero appear to have incorrect codes %d, %d\n",
                    mtf_core.ellipses[up].code, mtf_core.ellipses[down].code
                );
            } else {
                if (mtf_core.ellipses[up].code != 1) {
                    std::swap(up, down);
                }
                if (mtf_core.ellipses[up].code == 1) {
                    printf("up ellipse identified\n");
                } else {
                    printf("Warning: could not identify 'up' ellipse (codes are up=%d, down=%d)\n", 
                        mtf_core.ellipses[up].code, mtf_core.ellipses[down].code
                    );
                }
                Point2d dv(mtf_core.ellipses[up].centroid_x - zero.x, mtf_core.ellipses[up].centroid_y - zero.y);
                double dot = longitudinal.x * dv.x + longitudinal.y * dv.y;
                if (dot < 0) {
                    longitudinal = -longitudinal;
                }
                printf("Oriented transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
                printf("Oriented longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
                
                // project ellipses onto longitudinal vector
                dist_by_idx.clear();
                for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                    Point2d delta(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
                    double dot = longitudinal.x * delta.x + longitudinal.y * delta.y;
                    double tdot = transverse.x * delta.x + transverse.y * delta.y;
                    
                    
                    if (fabs(tdot) < mtf_core.ellipses[i].major_axis * 2) {
                        dist_by_idx.push_back(make_pair(dot, i));
                    } else {
                        int ccode = mtf_core.ellipses[i].code;
                        if (ccode == 0) {
                            printf("ellipse %d (%lf, %lf) not included in distance scale\n", (int)i, mtf_core.ellipses[i].centroid_x, mtf_core.ellipses[i].centroid_y);
                        } else {
                            // these are used for focal length calibration (and probably for homography too)
                            int quadrant = ((dot < 0) ? 1 : 0) | ((tdot < 0) ? 2 : 0);
                            if (quadrant <= 1) { // ugly, but swap quadrants 0 and 1
                                quadrant = 1 - quadrant;
                            }
                            // find the best match
                            for (int j=0; j < n_fiducials; j++) {
                                if (ccode == main_fiducials[j].code && quadrant == main_fiducials[j].quadrant) {
                                    printf("Ellipse %lu assigned to fiducial (code=%d, quad=%d) : (%lf %lf) -> (%lf %lf)\n",
                                        i, ccode, quadrant,
                                        mtf_core.ellipses[i].centroid_x, mtf_core.ellipses[i].centroid_y,
                                        main_fiducials[j].rcoords.x, main_fiducials[j].rcoords.y
                                    );
                                    auto it = fiducials.find(ccode);
                                    if (it == fiducials.end()) {
                                        fiducials[ccode] = vector<Fiducial>(4);
                                        fidmap[ccode] = vector<int>(4);
                                    }
                                    Fiducial& selected = fiducials[ccode][quadrant];
                                    selected = Fiducial(
                                        mtf_core.ellipses[i].centroid_x,
                                        mtf_core.ellipses[i].centroid_y,
                                        main_fiducials[j].rcoords.x, 
                                        main_fiducials[j].rcoords.y,  
                                        ccode,
                                        quadrant
                                    );
                                    fidmap[ccode][quadrant] = i;
                                }
                            }
                        }
                    }
                }
                printf("number of ellipses used in scale = %d, total #ellipses = %d\n", (int)dist_by_idx.size(), (int)mtf_core.ellipses.size());
                sort(dist_by_idx.begin(), dist_by_idx.end());
                
                // build the preliminary distance scale
                bool centre_done = false;
                for (size_t i=0; i < dist_by_idx.size(); i++) {
                    int j=dist_by_idx[i].second;
                    printf("(%lf, %lf), dist=%lf -> code=%d, axis ratio = %lf pos3d: %lf %lf %lf\n",
                        mtf_core.ellipses[j].centroid_x, mtf_core.ellipses[j].centroid_y,
                        dist_by_idx[i].first,
                        mtf_core.ellipses[j].code,
                        mtf_core.ellipses[j].minor_axis / mtf_core.ellipses[j].major_axis,
                        mtf_core.ellipses[j].pos1[0], mtf_core.ellipses[j].pos1[1], mtf_core.ellipses[j].pos1[2]
                        //mtf_core.ellipses[j].pos2[0], mtf_core.ellipses[j].pos2[1], mtf_core.ellipses[j].pos2[2]
                    );
                    
                    if (fabs(dist_by_idx[i].first) < 15) {
                        if (!centre_done) { 
                            distance_scale.push_back(cv::Point3d(0,0,0));
                            centre_done = true;
                            if (mtf_core.ellipses[j].code != 0) {    
                                printf("Weird. Central fiducial has incorrect code %d\n", mtf_core.ellipses[j].code);
                            }
                        }
                    } else {
                        
                        distance_scale.push_back(
                            
                            cv::Point3d(
                                dist_by_idx[i].first, 
                                
                                (((dist_by_idx[i].first < 0) ? -2 : 0) + mtf_core.ellipses[j].code) * (1.5*5) * 
                                 (dist_by_idx[i].first < 0 ? -1 : 1),
                                
                                mtf_core.ellipses[j].pos2[2] // z-position of fiducial
                            )
                        );
                    }
                }
                
                
                if (!pose_based) {
                    printf("not computing pose\n");
                    return;
                }
                
                prin = Point2d(mtf_core.img.cols/2.0, mtf_core.img.rows/2.0);
                
                // take 2
                vector<Eigen::Vector2d> feature_points(5);
                vector<Eigen::Vector3d> world_points(5);
                int fidcount = 0;
                Point2d pix_centroid(0, 0);
                int totalcount = 0;
                for (auto fidv: fiducials) {
                    for (auto fid: fidv.second) {
                        if (fidcount < 5) {
                            feature_points[fidcount] = Eigen::Vector2d(
                                (fid.icoords.x - prin.x),
                                (fid.icoords.y - prin.y)
                            );
                            world_points[fidcount] = Eigen::Vector3d(fid.rcoords.x, fid.rcoords.y, 1);
                            fidcount++;
                        }
                        pix_centroid.x += fid.icoords.x;
                        pix_centroid.y += fid.icoords.y;
                        totalcount++;
                    }
                }
                pix_centroid *= 1.0/double(totalcount);
                
                vector<Eigen::Matrix<double, 3, 4> > projection_matrices;
                vector<vector<double> > radial_distortions;
                cv::Mat rot_matrix = cv::Mat(3, 3, CV_64FC1);
                cv::Mat rod_angles = cv::Mat(3, 1, CV_64FC1);
                
                theia::FivePointFocalLengthRadialDistortion(
                    feature_points,
                    world_points,
                    1, // number of distortion parms
                    &projection_matrices,
                    &radial_distortions
                );
                
                double min_bpr = 1e50;
                int min_bpr_idx = 0;
                for (size_t k=0; k < projection_matrices.size(); k++) {
                    double bpr = 0;
                    for (int j=0; j < 5; j++) {
                        Eigen::VectorXd rp = projection_matrices[k]*Eigen::Vector4d(world_points[j][0], world_points[j][1], world_points[j][2], 1.0);
                        rp /= rp[2];
                        double rad = 1 + radial_distortions[k][0]*(rp[0]*rp[0] + rp[1]*rp[1]);
                        rp /= rad;
                        bpr += (feature_points[j] - Eigen::Vector2d(rp[0], rp[1])).norm();
                    }
                    bpr /= 5.0;
                    
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
                    
                    // only keep a rotation matrix if it repeats with the Rodrigues angle convention
                    if (bpr <= min_bpr && rot_err < 0.01) { 
                        min_bpr = bpr;
                        min_bpr_idx = k;
                    }
                }
                Eigen::MatrixXd P = projection_matrices[min_bpr_idx];
                distortion = radial_distortions[min_bpr_idx][0];
                
                double r1n = (P.block(0,0,1,3)).norm();
                double r3n = (P.block(2,0,1,3)).norm();
                double w = sqrt( (r3n*r3n)/(r1n*r1n) );
                focal_length = 1.0/w;
                
                P /= P.block(0,0,1,3).norm(); // remove the arbitrary scaling factor
                P.row(2) /= w; // remove focal length from last row
                
                Eigen::MatrixXd RM = P.block(0,0,3,3);
                Eigen::Vector3d Pcop = P.col(3);
                
                for (int rr=0; rr < 3; rr++) {
                    for (int cc=0; cc < 3; cc++) {
                        rot_matrix.at<double>(rr,cc) = RM(rr,cc);
                    }
                }
                cv::Rodrigues(rot_matrix, rod_angles);
                
                vector<Eigen::Vector2d> ba_img_points;
                vector<Eigen::Vector3d> ba_world_points;
                for (auto fidv: fiducials) {
                    for (auto fid: fidv.second) {
                        ba_img_points.push_back(Eigen::Vector2d(fid.icoords.x - prin.x, fid.icoords.y - prin.y));
                        ba_world_points.push_back(Eigen::Vector3d(fid.rcoords.x, fid.rcoords.y, 1.0));
                    }
                }
                
                for (size_t i=0; i < dist_by_idx.size(); i++) {
                    int j=dist_by_idx[i].second;
                    
                    if (fabs(dist_by_idx[i].first) >= 15) {
                        ba_img_points.push_back(Eigen::Vector2d(mtf_core.ellipses[j].centroid_x - prin.x, mtf_core.ellipses[j].centroid_y - prin.y));
                        
                        double xval = (((dist_by_idx[i].first < 0) ? -2 : 0) + mtf_core.ellipses[j].code) * (1.5*5) * (dist_by_idx[i].first < 0 ? -1 : 1);
                        ba_world_points.push_back(Eigen::Vector3d(xval, 20, 1.0));
                    }
                }
                
                Point2d delta_principal(0,0);
                Bundle_adjuster ba(
                    ba_img_points, ba_world_points, 
                    Pcop, 
                    rod_angles,
                    distortion,
                    w,
                    delta_principal.x, delta_principal.y
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
    
    
    Eigen::Vector3d backproject(double pixel_x, double pixel_y) const {
        double xs = (pixel_x - prin.x);
        double ys = (pixel_y - prin.y);
        
        //double rad = 1 + distortion*(xs*xs + ys*ys); // this won't actually work --- we have to use a fixed-point iteration to solve this
        //xs *= rad;
        //ys *= rad;
        
        Eigen::Vector3d dv(xs, ys, 1.0);
        dv = invP*dv;
        dv /= dv.norm();
        
        double s = (1 - cop[2])/dv[2];
        Eigen::Vector3d ip = cop + s*dv;
        
        // now we have ip in world coordinates, but we actually want it in camera 
        // coordinates
        
        // weird hack to (apparently) fix up z-axis
        ip = rotation*ip + Eigen::Vector3d(translation[0], translation[1], -translation[2]);
        
        return ip;
    }
    
    void estimate_depth(double pixel_x, double pixel_y, double& depth) const {
        Eigen::Vector3d bp = backproject(pixel_x, pixel_y);
        depth = bp[2] - centre_depth;
    }
    
    void estimate_depth(double pixel_offset, double& depth) {
    
        if (distance_scale.size() > 0) {
            const int tdim = 3;
            const int tpts = distance_scale.size();
            MatrixXd design(tpts, tdim);
            VectorXd zv(tpts);
            VectorXd zsol;
            VectorXd yv(tpts);
            VectorXd ysol;
            
            for (int row=0; row < tpts; row++) {
                double w = fabs(distance_scale[row].y) < 50.0 ? 1.1 : 1; // prefer centre points a little
                design(row, 0) = 1*w;
                double xval = distance_scale[row].x;
                double xprod = 1.0;
                for (int col=1; col < tdim; col++) {
                    xprod *= xval;
                    design(row, col) = xprod*w;
                }
                zv[row] = distance_scale[row].z*w;
                yv[row] = distance_scale[row].y*w;
            }
            zsol = design.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(zv);
            ysol = design.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(yv);
            
            double chart_z = 0;
            double focus_plane_fs = 0; // TODO: if the values are unbalanced???
            double xval = pixel_offset;
            double xprod = 1.0;
            for (int col=1; col < tdim; col++) {
                xprod *= xval;
                focus_plane_fs += zsol(col, 0) * xprod;
                chart_z += ysol(col, 0) * xprod;
            }
            printf("interpolated focus plane position = %lf, chart scale = %lf\n", chart_z, chart_scale);
            depth = focus_plane_fs;
        } else {
            double foreshortening = sqrt(0.5); // assume 45 degree chart if no scale is provided
            depth = pixel_offset * chart_scale * foreshortening; 
        }
    }
    
    
    
    Point2d zero;
    Point2d transverse;
    Point2d longitudinal;
    double chart_scale;
    vector<cv::Point3d> distance_scale;
    int largest_block_index;
    
    map<int, vector<Fiducial> > fiducials;
    map<int, vector<int> > fidmap;
    
    Point2d prin;
    double focal_length;
    double centre_depth;
    double distortion;
    Eigen::Matrix3d rotation;
    Eigen::Vector3d translation;
    
    Eigen::Matrix3d invP; // from image to world coords
    Eigen::Vector3d cop;  // centre of projection (i.e., camera)
};

#endif
