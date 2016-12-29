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
#ifndef CAMERA_DRAW_H
#define CAMERA_DRAW_H

#include "include/logger.h"
#include "include/common_types.h"
#include "include/distance_scale.h"

class Camera_draw {
  public:
    Camera_draw(const cv::Mat& img, const Distance_scale& distance_scale, double scaling_factor=1.0, int pad=0) 
    : img(img), distance_scale(distance_scale), scaling_factor(scaling_factor)
    {
        rimg = gray_to_colour(img, scaling_factor, pad);
        if (scaling_factor != 1.0) {
            this->distance_scale.img_scale *= scaling_factor;
            this->distance_scale.prin *= scaling_factor;
        }
    }
    
    void arc_with_arrow(int axis, double rad, double offset, cv::Scalar colour, bool tp=false) {
        vector<Point2d> cpts;
       
        double st = tp ? M_PI : 0;         // start angle
        double et = tp ? M_PI - M_PI*0.55: M_PI*0.55; // end angle
        double ainc = M_PI/200.0 * (tp ? -1 : 1);
        double x=0;
        double y=0;
        double z=0;
        for (double theta=st; (tp && theta > et) || (!tp && theta < et); theta += ainc) {
            switch(axis) {
            case 0: // arrow around x-axis
                x = 0;
                y = rad*cos(theta) - offset;
                z = -rad*sin(theta);
                break;
            case 1: // arrow around y-axis
                x = rad*cos(theta) - offset;
                y = 0;
                z = -rad*sin(theta);
                break;
            case 2: // arrow around z-axis
                x = rad*cos(theta) - offset;
                y = rad*sin(theta) - offset;
                z = -offset;
                break;
            }
            cpts.push_back(distance_scale.world_to_image(x, y, z));
        }
        curve(cpts, cv::Scalar(20,20,20), 4, cv::Scalar(20,20,20));
        curve(cpts, colour, 2, colour);
        {
            const double th = 0; // angle orientation, could be 180 to flip arrow?
            const double dth = 25.0/180.0*M_PI;
            const double aw=3.0;
            double rad = 10;
            
            Point2d dpts[3] = {
              { (rad+aw)*cos(th),     (rad+aw)*sin(th) },
              { (rad-aw)*cos(th+dth), (rad-aw)*sin(th+dth) },
              { (rad-aw)*cos(th-dth), (rad-aw)*sin(th-dth) }
            };
            
            cv::Point pts[3];
            for (int i=0; i < 3; i++) {
                Point2d pt;
                switch(axis) {
                case 0: // x-axis                    
                    pt = distance_scale.world_to_image(dpts[i].y, y - (tp ? -1 : 1)*(dpts[i].x - (rad - aw)), z);
                    break;
                case 1: // y-axis
                    pt = distance_scale.world_to_image(x - (tp ? -1 : 1)*(dpts[i].x - (rad - aw)), dpts[i].y, z);
                    break;
                case 2: // z-axis
                    pt = distance_scale.world_to_image(x - (tp ? -1 : 1)*(dpts[i].x - (rad - aw)), y - dpts[i].y, z);
                    break;
                }
                
                pts[i].x = lrint(pt.x);
                pts[i].y = lrint(pt.y);
            }
            for (int i=0; i < 3; i++) {
                cv::line(rimg, pts[i], pts[(i+1)%3], cv::Scalar(20,20,20), 2, CV_AA);
            }
            cv::fillConvexPoly(rimg, (const cv::Point*)&pts, 3, colour, CV_AA);
        } 
       
    }

    cv::Mat gray_to_colour(const cv::Mat& ig_img, double sfactor, int pad=0) {
        // TODO: we can downsize the image here ...
        cv::Mat channel;
        cv::Mat g_img;
        if (sfactor != 1.0) {
            cv::resize(ig_img, g_img, cv::Size(0,0), sfactor, sfactor);
        } else {
            g_img = ig_img.clone();
        }
        channel = cv::Mat(g_img.rows, g_img.cols, CV_8UC1);
       
        double imin;
        double imax;
        vector<uint64_t> histo(65537, 0);
        uint16_t* dptr = (uint16_t*)g_img.data;
        uint16_t* sentinel = dptr + g_img.rows*g_img.cols;
        while (dptr < sentinel) {
            histo[(*dptr)&0xffff]++;
            dptr += 8; // only sample every 8th pixel to speed things up??
        }
        uint64_t sum = histo[0];
        for (size_t i=1; i < histo.size(); i++) {
            histo[i-1] = sum;
            sum += histo[i];
        }
        histo.back() = sum;
        // find 2% max brightness
        uint64_t target = lrint(0.98*sum);
        imax = histo.size() - 1;
        while (imax > 0 && histo[imax] > target) imax--;
        target = lrint(0.02*sum);
        imin = 0;
        while (imin < imax && histo[imin] < target) imin++;
        double scale = 255.0/(imax - imin);
        double offset = -imin*scale;
        g_img.convertTo(channel, CV_8U, scale, offset);
       
        vector<cv::Mat> channels;
        channels.push_back(channel);
        channels.push_back(channel);
        channels.push_back(channel);
        cv::Mat merged;
        merge(channels, merged);
        
        if (pad > 0) {
            // TODO: maybe merge right column if image is too narrow
            merged.resize(merged.rows + 220); // TODO: adjustable padding size?
            rectangle(merged, Point2d(0, g_img.rows), Point2d(merged.cols, merged.rows), cv::Scalar::all(255), CV_FILLED);
        }
       
        return merged;
    }

    void curve(const vector<Point2d>& data, cv::Scalar col, double width, cv::Scalar col2=cv::Scalar(0, 0, 0)) {
        int prevx = 0;
        int prevy = 0;
       
        double total_l = 0;
        for (size_t i=1; i < data.size(); i++) {
            double dx = data[i].x - data[i-1].x;
            double dy = data[i].y - data[i-1].y;
            total_l += sqrt(dx*dx + dy*dy);
        } 
       
        bool shade = col2[0] != 0 || col2[1] != 0 || col2[1] != 0;
        cv::Scalar blended_col = col;
       
        cv::Rect bounds(0, 0, img.cols, img.rows+5);
       
        double running_l = 0;
        for (size_t i=0; i < data.size(); i++) {
            if (i > 1) {
                double dx = data[i].x - data[i-1].x;
                double dy = data[i].y - data[i-1].y;
                running_l += sqrt(dx*dx + dy*dy);
            }
           
            double progress = running_l / total_l;
           
            if (shade) {
                for (int k=0; k < 3; k++) {
                    blended_col[k] = col[k] + (col2[k] - col[k])*progress;
                }
            }
           
            cv::Point head(prevx, prevy);
            cv::Point tail(lrint(data[i].x), lrint(data[i].y));
            bool inside = cv::clipLine(bounds, head, tail);
           
            if (i > 0 && inside) {
                cv::line(rimg, head, tail, blended_col, width, CV_AA);
            }
           
            prevx = tail.x;
            prevy = tail.y;
        }
    }
   
    void chart_centre_marker(void) {
        cv::Scalar mark_col(0, 127-20, 255-20);
        Point2d czero = distance_scale.world_to_image(0, 0);
        
        for (int i=0; i < 4; i++) {
            Point2d wp(10*cos(2.0*i*M_PI/4.0), 10*sin(2.0*i*M_PI/4.0));
            Point2d pt = distance_scale.world_to_image(wp.x, wp.y);
            cv::line(rimg, czero, pt, mark_col, 2, CV_AA);
        }
        cv::circle(rimg, czero, 10, mark_col, 2, CV_AA);
        
        int j=0;
        for (double th=0; j < 4; th += 2*M_PI/4.0, j++) {
            
            const double dth = 10.0/180.0*M_PI;
            double rad = 10;
            const double aw=2.0;
            
            Point2d dpts[3] = {
              { (rad+aw)*cos(th),     (rad+aw)*sin(th) },
              { (rad-aw)*cos(th+dth), (rad-aw)*sin(th+dth) },
              { (rad-aw)*cos(th-dth), (rad-aw)*sin(th-dth) }
            };
            
            cv::Point pts[3];
            for (int i=0; i < 3; i++) {
                Point2d pt = distance_scale.world_to_image(dpts[i].x, dpts[i].y);
                pts[i].x = lrint(pt.x);
                pts[i].y = lrint(pt.y);
            }
            
            cv::fillConvexPoly(rimg, (const cv::Point*)&pts, 3, mark_col, CV_AA);
        }
    }
    
    void camera_centre_marker(double dc_x, double dc_y) {
        cv::Scalar reticle_col(0, 127-50, 255-50);
        cv::Point centre(dc_x, dc_y);
        
        
            
        const double rad = 25;
        cv::circle(rimg, centre, rad, cv::Scalar(20,20,20), 4, CV_AA);
        int i=0;
        for (double th=M_PI/2; i < 4; th += 2*M_PI/4.0, i++) {
            const double dth = 12.5/180.0*M_PI;
            
            cv::Point pts[3] = {
              { int((double)centre.x + (rad-9)*cos(th)),     int((double)centre.y + (rad-9)*sin(th)) },
              { int((double)centre.x + (rad+8)*cos(th+dth)), int((double)centre.y + (rad+8)*sin(th+dth)) },
              { int((double)centre.x + (rad+8)*cos(th-dth)), int((double)centre.y + (rad+8)*sin(th-dth)) }
            };
            
            cv::fillConvexPoly(rimg, (const cv::Point*)&pts, 3, cv::Scalar(20,20,20), CV_AA);
        }
        cv::circle(rimg, centre, rad, reticle_col, 2, CV_AA);
        i=0;
        for (double th=M_PI/2; i < 4; th += 2*M_PI/4.0, i++) {
            const double dth = 10.0/180.0*M_PI;
            
            cv::Point pts[3] = {
              { int((double)centre.x + (rad-7)*cos(th)),     int((double)centre.y + (rad-7)*sin(th)) },
              { int((double)centre.x + (rad+7)*cos(th+dth)), int((double)centre.y + (rad+7)*sin(th+dth)) },
              { int((double)centre.x + (rad+7)*cos(th-dth)), int((double)centre.y + (rad+7)*sin(th-dth)) }
            };
            
            cv::fillConvexPoly(rimg, (const cv::Point*)&pts, 3, reticle_col, CV_AA);
        }
    }
    
    void alpha_block(const Point2d& p, const cv::Size& s, const cv::Scalar& col, double alpha) {
        // trim to img bounds
        int e_width = std::min((int)p.x + (int)s.width, rimg.cols-1) - p.x;
        int s_row = std::max(0, (int)p.y - (int)s.height - 1);
        int e_row = std::min((int)p.y + 2, rimg.rows-1);
        cv::Mat roi = rimg(cv::Rect(p.x, s_row, e_width, e_row - s_row));
        cv::Mat cblock(roi.size(), CV_8UC3, col); 
        cv::addWeighted(cblock, alpha, roi, 1.0 - alpha , 0.0, roi); 
    }
    
    const cv::Mat img;
    Distance_scale distance_scale;
    double scaling_factor;
    
    cv::Mat rimg; // rendered output image
};


    

#endif
