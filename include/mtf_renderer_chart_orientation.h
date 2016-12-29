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
#ifndef MTF_RENDERER_CHART_ORIENTATION_H
#define MTF_RENDERER_CHART_ORIENTATION_H

#include "include/logger.h"
#include "mtf_renderer.h"
#include "common_types.h"
#include "distance_scale.h"
#include "include/camera_draw.h"

class Mtf_renderer_chart_orientation : public Mtf_renderer {
  public:
    Mtf_renderer_chart_orientation(
        const std::string& img_filename,
        const std::string& wdir, 
        const std::string& co_fname, 
        const cv::Mat& img,
        int gnuplot_width,
        Distance_scale& distance_scale,
        cv::Rect* dimension_correction = NULL)
      :  Mtf_renderer(img_filename),
         wdir(wdir), co_fname(co_fname), 
         img(img), gnuplot_width(gnuplot_width),
         dimension_correction(dimension_correction),
         draw(img, distance_scale, (gnuplot_width / double(img.cols))) {
    
      
    }
    
    void render(const vector<Block>& ) {
    
        Distance_scale& distance_scale = draw.distance_scale;
        
        draw.chart_centre_marker();        
        
        // draw centre-of-camera marker
        draw.camera_centre_marker(
            draw.scaling_factor*(dimension_correction->width/2 - dimension_correction->x), 
            draw.scaling_factor*(dimension_correction->height/2 - dimension_correction->y)
        );
        
        vector<Point2d> curve;
        
        cv::Scalar c_dark(20, 20, 20);
        cv::Scalar c_red(30, 30, 255);
        cv::Scalar c_green(30, 255, 30);
        cv::Scalar c_blue(255, 30, 30);
        
        cv::Scalar c_lred(80, 80, 255);
        cv::Scalar c_lgreen(130, 255, 130);
        cv::Scalar c_lblue(255, 182, 0);
        
        curve.clear();
        for (double ystep=-100; ystep <= 0; ystep += 2) {
            curve.push_back(distance_scale.world_to_image(-100, ystep));
        }
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_green, 2, c_green);
        curve.clear();
        for (double xstep=-100; xstep <= 0; xstep += 2) {
            curve.push_back(distance_scale.world_to_image(xstep, -100));
        }
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_red, 2, c_red);
        curve.clear();
        curve.push_back(distance_scale.world_to_image(-100, -100));
        curve.push_back(distance_scale.world_to_image(-100, -100, -100));
        draw.curve(curve, c_dark, 3, c_dark);
        draw.curve(curve, c_lblue, 2, c_lblue);
        
        draw.arc_with_arrow(0, 10, 100, c_lred, distance_scale.pitch_angle < 0);
        draw.arc_with_arrow(1, 10, 100, c_lgreen, distance_scale.yaw_angle < 0);
        draw.arc_with_arrow(2, 10, 100, c_lblue, distance_scale.roll_angle < 0);

        int font = cv::FONT_HERSHEY_DUPLEX; 
        char tbuffer[1024];
        int baseline;
        cv::Size ts;
        
        // TODO: degree symbol. sigh.
        sprintf(tbuffer, "Roll=%.2lf", distance_scale.roll_angle/M_PI*180);
        Point2d textpos = distance_scale.world_to_image(-115, -115, -105);
        ts = cv::getTextSize(tbuffer, font, 1, 3, &baseline);
        draw.alpha_block(textpos, ts, CV_RGB(255, 255, 255), 0.5);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(50, 50, 50), 3, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, c_lblue, 1, CV_AA);
        
        sprintf(tbuffer, "Yaw=%.2lf", distance_scale.yaw_angle/M_PI*180);
        textpos = distance_scale.world_to_image(-105, 10, 0);
        ts = cv::getTextSize(tbuffer, font, 1, 3, &baseline);
        draw.alpha_block(textpos, ts, CV_RGB(255, 255, 255), 0.5);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(50, 50, 50), 3, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, c_lgreen, 1, CV_AA);
        
        sprintf(tbuffer, "Pitch=%.2lf", distance_scale.pitch_angle/M_PI*180);
        textpos = distance_scale.world_to_image(10, -105, 0);
        ts = cv::getTextSize(tbuffer, font, 1, 3, &baseline);
        draw.alpha_block(textpos, ts, CV_RGB(255, 255, 255), 0.5);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(50, 50, 50), 3, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(draw.rimg, tbuffer, textpos, font, 1, c_lred, 1, CV_AA);
        
        imwrite(wdir + '/' + co_fname, draw.rimg);
    }
    
  private:
    string wdir;
    string co_fname; 
    const cv::Mat& img;            
    
    int gnuplot_width;
    
    cv::Rect* dimension_correction;
    
    Camera_draw draw;
};

#endif
