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

class Distance_scale {
  public:
    Distance_scale(void)
    : chart_scale(1.0), largest_block_index(-1) {
    }
    
    void construct(Mtf_core& mtf_core) {
    
        // define reference frame from three solid circular dots in centre of chart
        printf("#solid ellipses: %d\n", (int)mtf_core.solid_ellipses.size());
        for (size_t j=0; j < mtf_core.solid_ellipses.size(); j++) {
            printf("(%lf %lf)\n", mtf_core.solid_ellipses[j].x, mtf_core.solid_ellipses[j].y);
        }
        printf("\n\n");
        if (mtf_core.solid_ellipses.size() >= 3) {
            printf("got all three solid ellipses\n");
            double maxdelta = 0;
            Point first;
            Point last;
            for (int i=0; i <= 1; i++) {
                for (int j=i+1; j < 3; j++) {
                    Point d = mtf_core.solid_ellipses[i] - mtf_core.solid_ellipses[j];
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
            Point base = first - last;
            chart_scale = 91.0 / sqrt(base.x*base.x + base.y*base.y);
            zero = Point(0.5*(first.x+last.x), 0.5*(first.y+last.y));
            printf("absolute centre: (%lf, %lf)\n", zero.x, zero.y);
            transverse = normalize(first - last);
            printf("transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
            // sign of transverse unknown at this stage
            longitudinal = Point(-transverse.y, transverse.x);
            printf("longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
            printf("chart scale is %lf mm per pixels\n", chart_scale);
            
            // now try to decode the ellipses
            vector< pair<double, int> > dist_by_idx;
            for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                Ellipse_decoder ed(mtf_core.ellipses[i], mtf_core.img, transverse);
                mtf_core.ellipses[i].set_code(ed.code);
                if (!mtf_core.ellipses[i].solid) {
                    Point d(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
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
                Point dv(mtf_core.ellipses[up].centroid_x - zero.x, mtf_core.ellipses[up].centroid_y - zero.y);
                double dot = longitudinal.x * dv.x + longitudinal.y * dv.y;
                if (dot < 0) {
                    longitudinal = -longitudinal;
                }
                printf("Oriented transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
                printf("Oriented longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
                
                // project ellipses onto longitudinal vector
                dist_by_idx.clear();
                for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                    Point delta(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
                    double dot = longitudinal.x * delta.x + longitudinal.y * delta.y;
                    double tdot = transverse.x * delta.x + transverse.y * delta.y;
                    
                    if (fabs(tdot) < mtf_core.ellipses[i].major_axis * 2) {
                        dist_by_idx.push_back(make_pair(dot, i));
                    } else {
                        printf("ellipse %d (%lf, %lf) not included in distance scale\n", (int)i, mtf_core.ellipses[i].centroid_x, mtf_core.ellipses[i].centroid_y);
                    }
                }
                printf("number of ellipses used in scale = %d, total #ellipses = %d\n", (int)dist_by_idx.size(), (int)mtf_core.ellipses.size());
                sort(dist_by_idx.begin(), dist_by_idx.end());
                bool centre_done = false;
                for (size_t i=0; i < dist_by_idx.size(); i++) {
                    int j=dist_by_idx[i].second;
                    printf("(%lf, %lf), dist=%lf -> code=%d, axis ratio = %lf\n",
                        mtf_core.ellipses[j].centroid_x, mtf_core.ellipses[j].centroid_y,
                        dist_by_idx[i].first,
                        mtf_core.ellipses[j].code,
                        mtf_core.ellipses[j].minor_axis / mtf_core.ellipses[j].major_axis
                    );
                    if (fabs(dist_by_idx[i].first) < 5) {
                        if (!centre_done) { 
                            distance_scale.push_back(cv::Point3d(0,0,0));
                            centre_done = true;
                        }
                    } else {
                        distance_scale.push_back(
                            cv::Point3d(
                                dist_by_idx[i].first, 
                                (((dist_by_idx[i].first < 0) ? -2 : 0) + mtf_core.ellipses[j].code) * (1.5*5) * 
                                 (dist_by_idx[i].first < 0 ? -1 : 1),
                                mtf_core.ellipses[j].minor_axis/mtf_core.ellipses[j].major_axis
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
                
                Point median(xp[xp.size()/2], yp[yp.size()/2]); 
                
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
                zero = Point(932, 710);
                transverse = Point(1,0);
                longitudinal = Point(0,1);
            }
            printf("zero: %lf %lf\n", zero.x, zero.y);
            printf("transverse: %lf %lf\nlongitudinal: %lf %lf\n", transverse.x, transverse.y, longitudinal.x, longitudinal.y);
            
        }
    }
    
    void estimate_depth(double pixel_offset, double& depth) {
    
        if (distance_scale.size() > 0) {
            // find the two centre-most scale markers, average their distance to estimate chart angle
            int middle = 0;
            for (int i=1; i < (int)distance_scale.size(); i++) {
                if (fabs(distance_scale[i].x) < fabs(distance_scale[middle].x)) {
                    middle = i;
                }
            }
            double foreshortening = 0.5*(fabs(distance_scale[middle-1].x) + fabs(distance_scale[middle+1].x));
            foreshortening *= chart_scale/fabs(distance_scale[middle-1].y);
            
            printf("foreshortening: %lf\n", foreshortening);
            
            // pixel_offset is in pixels, relative to centre of chart
            int scale_lower=0;
            while (scale_lower < (int)distance_scale.size() - 2 &&
                   distance_scale[scale_lower].x < pixel_offset) {
                scale_lower++;
            }
            printf("scale limits: %d, %d : %lf, %lf\n", 
                scale_lower, scale_lower+1, 
                distance_scale[scale_lower].x, distance_scale[scale_lower+1].x
            );
            
            const int tdim = 2;
            const int tpts = 4;
            MatrixXd design(tpts, tdim);
            VectorXd yv(tpts);
            
            int lindex = scale_lower;
            int uindex = scale_lower+1;
            // grab the 4 closest points
            int grabbed = 2; // scale_lower, and scale_lower + 1
            while (grabbed < tpts) {
                if (lindex >= 1) {
                    grabbed++;
                    lindex--;
                }
                if (grabbed < tpts) {
                    if (uindex < (int)distance_scale.size() - 1) {
                        uindex++;
                        grabbed++;
                    }
                }
            }
            for (int row=0; row < tpts; row++) {
                design(row, 0) = 1;
                double xval = distance_scale[lindex].x;
                double xprod = 1.0;
                for (int col=1; col < tdim; col++) {
                    xprod *= xval;
                    design(row, col) = xprod;
                }
                yv[row] = distance_scale[lindex].y;
                lindex++;
            }
            VectorXd sol = design.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(yv);
            
            double focus_plane_position = sol(0,0);
            double xval = pixel_offset;
            double xprod = 1.0;
            for (int col=1; col < tdim; col++) {
                xprod *= xval;
                focus_plane_position += sol(col, 0) * xprod;
            }
            depth = focus_plane_position * foreshortening;
        } else {
            double foreshortening = sqrt(0.5); // assume 45 degree chart if no scale is provided
            depth = pixel_offset * chart_scale * foreshortening; 
        }
    }
    
    Point zero;
    Point transverse;
    Point longitudinal;
    double chart_scale;
    vector<cv::Point3d> distance_scale;
    int largest_block_index;
    
};

#endif
