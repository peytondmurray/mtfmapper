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
    : chart_scale(1.0), largest_block_index(-1), chart_normal(Eigen::Vector3d::Zero()), chart_constant(0), focal_length(10000) {
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
                    
                    if (norm(Point2d(mtf_core.ellipses[i].centroid_x, mtf_core.ellipses[i].centroid_y) - zero) < 10) {
                        zero_pos3d = mtf_core.ellipses[i].pos1;
                        printf("central ellipse %lf %lf %lf\n", zero_pos3d[0], zero_pos3d[1], zero_pos3d[2]);
                    }
                    
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
                
                Point2d pix_centroid(0, 0);
                Eigen::Vector3d centroid(0, 0, 0);
                int fidcount = 0;
                for (auto fidid: fidmap) {
                    for (auto fid: fidid.second) {
                        mtf_core.ellipses[fid].pose(mtf_core.img.rows, mtf_core.img.cols, 10.0);
                        centroid += mtf_core.ellipses[fid].pos1;    
                        pix_centroid += Point2d(mtf_core.ellipses[fid].centroid_x, mtf_core.ellipses[fid].centroid_y);
                        fidcount++;
                    }
                }
                centroid *= 1.0/double(fidcount);
                pix_centroid *= 1.0/double(fidcount);
                printf("1st plane centroid: [%lf %lf %lf]\n", centroid[0], centroid[1], centroid[2]);
                
                MatrixXd plane(3, fidcount);
                int plane_idx = 0;
                for (auto fidid: fidmap) {
                    for (auto fid: fidid.second) {
                         plane.col(plane_idx++) = mtf_core.ellipses[fid].pos1 - centroid;    
                    }
                }
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(plane, Eigen::ComputeThinU | Eigen::ComputeThinV);
                chart_normal = svd.matrixU().col(2);
                
                printf("Chart normal: [%lf %lf %lf]\n", chart_normal[0], chart_normal[1], chart_normal[2]);
                
                
                printf("Fiducials:\n");
                for (auto fidv: fiducials) {
                    for (auto fid: fidv.second) {
                        printf("%lf %lf ; %lf %lf 0 ; %d\n", 
                            (fid.icoords.x - prin.x)/prin.x, (fid.icoords.y - prin.y)/prin.y,
                            fid.rcoords.x, fid.rcoords.y,
                            fidv.first
                        );
                    }
                }
               
                double f =  estimate_focal_length_backproject(mtf_core);
                
                printf("best focal length = %lf\n", f);
                
                focal_length = f;
                
                centroid.setZero();
                pix_centroid = Point2d(0,0);
                fidcount = 0;
                for (auto fidid: fidmap) {
                    for (auto fid: fidid.second) {
                        mtf_core.ellipses[fid].pose(mtf_core.img.rows, mtf_core.img.cols, 10.0, f, chart_normal);
                        centroid += mtf_core.ellipses[fid].pos1;    
                        pix_centroid += Point2d(mtf_core.ellipses[fid].centroid_x, mtf_core.ellipses[fid].centroid_y);
                        fidcount++;
                    }
                }
                centroid *= 1.0/double(fidcount);
                pix_centroid *= 1.0/double(fidcount);
                printf("plane centroid: [%lf %lf %lf]\n", centroid[0], centroid[1], centroid[2]);
                
                plane = MatrixXd(3, fidcount);
                plane_idx = 0;
                for (auto fidid: fidmap) {
                    for (auto fid: fidid.second) {
                         plane.col(plane_idx++) = mtf_core.ellipses[fid].pos1 - centroid;    
                    }
                }
                svd = Eigen::JacobiSVD<Eigen::MatrixXd>(plane, Eigen::ComputeThinU | Eigen::ComputeThinV);
                
                printf("and singular values are:\n");
                std::cout << svd.singularValues() << std::endl;
                
                chart_normal = svd.matrixU().col(2);
                
                printf("Chart normal: [%lf %lf %lf]\n", chart_normal[0], chart_normal[1], chart_normal[2]);
                
                chart_constant = chart_normal.dot(centroid);
                printf("chart plane const: %lf\n", chart_constant);
                
                
                Eigen::Vector3d centre_back = backproject(pix_centroid.x, pix_centroid.y);
                centre_depth = centre_back[2];
                
                Ellipse_detector& elA = mtf_core.ellipses[fidmap[13][0]];
                Eigen::Vector3d ela_pos = backproject(elA.centroid_x, elA.centroid_y);
                printf("centre depth: %lf\n", ela_pos[2]);
                
                // take 2
                vector<Eigen::Vector2d> feature_points(5);
                vector<Eigen::Vector3d> world_points(5);
                fidcount = 0;
                
                for (auto fidv: fiducials) {
                    for (auto fid: fidv.second) {
                        if (fidcount < 5) {
                            feature_points[fidcount] = Eigen::Vector2d(
                                (fid.icoords.x - prin.x),
                                (fid.icoords.y - prin.y)
                            );
                            world_points[fidcount] = Eigen::Vector3d(fid.rcoords.x, fid.rcoords.y, 1);
                            printf("%lf %lf ; %lf %lf %lf\n", 
                                feature_points[fidcount][0], feature_points[fidcount][1],
                                world_points[fidcount][0], world_points[fidcount][1], world_points[fidcount][2]
                            );
                            fidcount++;
                        }
                    }
                }
                
                vector<Eigen::Matrix<double, 3, 4> > projection_matrices;
                vector<vector<double> > radial_distortions;
                
                bool fps = theia::FivePointFocalLengthRadialDistortion(
                    feature_points,
                    world_points,
                    1, // number of distortion parms
                    &projection_matrices,
                    &radial_distortions
                );
                
                double min_bpr = 1e50;
                int min_bpr_idx = 0;
                for (size_t k=0; k < projection_matrices.size(); k++) {
                    printf("matrix %lu:\n", k);
                    std::cout << projection_matrices[k] << std::endl;
                    printf("radial distortion coef: %lg\n", radial_distortions[k][0]);
                    printf("backprojection:\n");
                    double bpr = 0;
                    for (int j=0; j < 5; j++) {
                        Eigen::VectorXd rp = projection_matrices[k]*Eigen::Vector4d(world_points[j][0], world_points[j][1], world_points[j][2], 1.0);
                        rp /= rp[2];
                        // would the sign of radial distorion be flipped?
                        double rad = 1 - radial_distortions[k][0]*(rp[0]*rp[0] + rp[1]*rp[1]);
                        rp /= rad;
                        printf("(%lf %lf)\n", rp[0], rp[1]);
                        
                        bpr += (feature_points[j] - Eigen::Vector2d(rp[0], rp[1])).norm();
                    }
                    bpr /= 5.0;
                    printf("back-projection MAD error: %lf pixels\n", bpr);
                    if (bpr < min_bpr) {
                        min_bpr = bpr;
                        min_bpr_idx = k;
                    }
                }
                Eigen::MatrixXd P = projection_matrices[min_bpr_idx];
                distortion = radial_distortions[min_bpr_idx][0];
                
                double r1n = (P.block(0,0,1,3)).norm();
                double r3n = (P.block(2,0,1,3)).norm();
                double w = sqrt( (r3n*r3n)/(r1n*r1n) );
                printf("5pt focal len = %lf\n", 1.0/w);
                
                focal_length = 1.0/w;
                
                Eigen::MatrixXd RM = P.block(0,0,3,3);
                RM.row(2) /= w;
                
                Eigen::Vector3d cn = RM.col(2)/(RM.col(2).norm());
                printf("proposed chart normal: [%lf %lf %lf]\n", cn[0], cn[1], cn[2]);
                chart_normal = cn;
                
                chart_constant = chart_normal.dot(centroid);
                printf("updated chart plane const: %lf\n", chart_constant);
                
                centre_back = backproject(pix_centroid.x, pix_centroid.y);
                centre_depth = centre_back[2];
                printf("chart centre distance: %lf\n", centre_back[2]);
                
                printf("recomputing centre markers...\n");
                distance_scale.clear();
                centre_done = false;
                for (size_t i=0; i < dist_by_idx.size(); i++) {
                    int j=dist_by_idx[i].second;
                    mtf_core.ellipses[j].pose(mtf_core.img.rows, mtf_core.img.cols, 5.0, f, chart_normal);
                    printf("%lf  %lf dist= %lf -> code=%d, axis ratio = %lf zpos3d: %lf %lf %lf\n",
                        mtf_core.ellipses[j].centroid_x, mtf_core.ellipses[j].centroid_y,
                        dist_by_idx[i].first,
                        mtf_core.ellipses[j].code,
                        mtf_core.ellipses[j].minor_axis / mtf_core.ellipses[j].major_axis,
                        mtf_core.ellipses[j].pos1[0], mtf_core.ellipses[j].pos1[1], mtf_core.ellipses[j].pos1[2]
                    );
                    
                    if (fabs(dist_by_idx[i].first) < 15) {
                        
                        
                        if (!centre_done) { 
                            distance_scale.push_back(cv::Point3d(0,0, mtf_core.ellipses[j].pos1[2]));
                            centre_done = true;
                            if (mtf_core.ellipses[j].code != 0) {    
                                printf("Weird. Central fiducial has incorrect code %d\n", mtf_core.ellipses[j].code);
                            }
                            zero_pos3d = mtf_core.ellipses[j].pos1;
                            printf("re-setting the zero_pos variable to [%lf %lf %lf]\n", zero_pos3d[0], zero_pos3d[1], zero_pos3d[2]);       
                        }
                    } else {
                        
                        distance_scale.push_back(
                            
                            cv::Point3d(
                                dist_by_idx[i].first, 
                                
                                (((dist_by_idx[i].first < 0) ? -2 : 0) + mtf_core.ellipses[j].code) * (1.5*5) * 
                                 (dist_by_idx[i].first < 0 ? -1 : 1),
                                
                                mtf_core.ellipses[j].pos1[2] // z-position of fiducial
                            )
                            
                        );
                    }
                }
                
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
        
        /*
        rp /= rp[2];
        // would the sign of radial distorion be flipped?
        double rad = 1 - radial_distortions[k][0]*(rp[0]*rp[0] + rp[1]*rp[1]);
        rp /= rad;
        */
        
        double rad = 1 + distortion*(xs*xs + ys*ys);
        xs *= rad;
        ys *= rad;
        
        if (xs == 0 || ys == 0) {
            printf("Trying to estimate world coordinates of principal point. Undefined?\n");
            return Eigen::Vector3d(0,0,0);
        }
        
        Eigen::Vector3d pv(1.0, ys/xs, focal_length/xs);
        double xw = chart_constant/(chart_normal.dot(pv));
        double zw = focal_length*xw/xs;
        double yw = ys*zw/focal_length;    
        
        return Eigen::Vector3d(xw, yw, zw);
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
            printf("fps=%lf, zps=%lf, zsol[0]=%lf\n", focus_plane_fs, zero_pos3d[2], zsol[0]);
            depth = focus_plane_fs;
            
            FILE* res = fopen("residuals.txt", "wt");
            double err = 0;
            double wsum = 0;
            for (int row=0; row < tpts; row++) {
                double w = fabs(distance_scale[row].y) < 50.0 ? 1.1 : 1; // prefer centre points a little
                double zval = zsol[0];
                double xval = distance_scale[row].x;
                double xprod = 1.0;
                for (int col=1; col < tdim; col++) {
                    xprod *= xval;
                    zval += zsol(col, 0) * xprod;
                }
                double e = distance_scale[row].z - zval;
                fprintf(res, "%lf\n", e);
                err += e*e*w;
                wsum += w;
            }
            fclose(res);
            printf("z-fit RMSE = %lf\n", sqrt(err/double(wsum)));
            //depth = pixel_offset * chart_scale * asin(chart_normal[0]);
        } else {
            double foreshortening = sqrt(0.5); // assume 45 degree chart if no scale is provided
            depth = pixel_offset * chart_scale * foreshortening; 
        }
    }
    
    double efl_eval(double focal, Mtf_core& mtf_core) {
        double sumerr = 0;
        double maxerr = 0;
        //double perr = 0;
        int idx = 0;
        vector<double> errors;
        for (auto f: fidmap) {
            // recalculate 
            vector<Eigen::Vector3d> pos(4);
            vector<Point2d> realpos(4);
            for (size_t fi = 0; fi < f.second.size(); fi++) {
                mtf_core.ellipses[f.second[fi]].pose(
                    mtf_core.img.rows, mtf_core.img.cols, 
                    10.0, // TODO: real fiducial diameter in mm, should be configure from elsewhere
                    focal,
                    chart_normal
                ); 
                pos[fi] = mtf_core.ellipses[f.second[fi]].pos1;
                realpos[fi] = fiducials[f.first][fi].rcoords;
            }
            for (size_t fi = 0; fi < pos.size() - 1; fi++) {
                for (size_t si = fi+1; si < pos.size(); si++) {
                    double rdist = norm(realpos[fi] - realpos[si]);
                    double err = fabs((pos[fi] - pos[si]).norm() - rdist);
                    sumerr += err*err;
                    maxerr = std::max(maxerr, err); // seems to work well for this problem
                    errors.push_back(err*err*err*err);
                    idx++;
                }
            }
            /*
            // assume 4 fids per group:
            Eigen::Vector3d dab = (pos[0] - pos[1])/(pos[0] - pos[1]).norm();
            Eigen::Vector3d dcd = (pos[2] - pos[3])/(pos[2] - pos[3]).norm();
            Eigen::Vector3d dbc = (pos[2] - pos[1])/(pos[2] - pos[1]).norm();
            Eigen::Vector3d dad = (pos[3] - pos[0])/(pos[3] - pos[0]).norm();
            perr += fabs(dab.dot(dbc));
            perr += fabs(dab.dot(dad));
            perr += fabs(dcd.dot(dbc));
            perr += fabs(dcd.dot(dad));
            */
        }
        sort(errors.begin(), errors.end());
        double meanmax = (errors[errors.size()-2] + errors[errors.size()-1]) * 0.5;
        fprintf(stderr, "%lf %lf %lf\n", focal, meanmax, sumerr);
        return meanmax;
    }
    
    double estimate_focal_length_backproject(Mtf_core& mtf_core) {
        // bracket the minimum
        double xmin = 500; // really low estimate for focal length
        double xmax = 200000;
        double dip_z = 1e50;
        double dip_x = (xmin + xmax)*0.5;
        double step = (xmax - xmin)/20.0;
        for (double x=xmin; x <= xmax; x += step) {
            double z = efl_eval(x, mtf_core);
            if (z < dip_z) {
                dip_x = x;
                dip_z = z;
            }
        }
        
        // golden section search
        const double phi = 0.61803398874989;
        double lower = std::max(xmin, dip_x - 2*step);
        double upper = dip_x + 2*step;
        double c = upper - phi*(upper - lower);
        double d = lower + phi*(upper - lower);
        const double tol = 0.01; // no need to go more accurate than this
        while ((upper - lower) > tol) {
            double fc = efl_eval(c, mtf_core);
            double fd = efl_eval(d, mtf_core);
            if (fc < fd) {
                upper = d;
                d = c;
                c = upper - phi*(upper - lower);
            } else {
                lower = c;
                c = d;
                d = lower + phi*(upper - lower);
            }
        }
        printf("minimum focal length error = %lf\n", efl_eval(0.5*(upper+lower), mtf_core));
        return 0.5*(upper + lower);
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
    Eigen::Vector3d zero_pos3d;
    Eigen::Vector3d chart_normal;
    double chart_constant;
    double focal_length;
    double centre_depth;
    double distortion;
};

#endif
