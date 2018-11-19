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
#ifndef MTF_RENDERER_ANNOTATE_H
#define MTF_RENDERER_ANNOTATE_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_annotate : public Mtf_renderer {
  public:
    Mtf_renderer_annotate(const cv::Mat& in_img, const std::string& fname,
      bool lpmm_mode, double pixel_size) 
      : img(in_img), ofname(fname), lpmm_mode(lpmm_mode), pixel_size(pixel_size) {
      
          unsigned int max_val = 0;
          for (int r=0; r < img.rows; r++) {
              for (int c=0; c < img.cols; c++) {
                  if (img.at<uint16_t>(r,c) > max_val) {
                      max_val = img.at<uint16_t>(r,c);
                  }
              }
          }
          cv::Mat temp_img;
          img.convertTo(temp_img, CV_8UC3, 255.0/double(max_val));
          cv::merge(vector<cv::Mat>(3, temp_img), out_img);
    }
    
    void render(const vector<Block>& blocks) {
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0) {
                    Point2d cent = blocks[i].get_edge_centroid(k);
                    write_number(out_img, lrint(cent.x), lrint(cent.y), val, blocks[i].get_quality(k));
                }
            }
        }    
        
        imwrite(ofname, out_img);
    }
    
    void write_number(cv::Mat& img, int px, int py, double val, double quality) {
        char buffer[10];
        sprintf(buffer, "%.2lf", val);
        
        int baseline = 0;
        
        cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, 0.5, 1, &baseline);
        
        ts.width -= 8;  // tweak box size slightly
        ts.height += 6; // tweak box height slighly
        
        cv::Point2d to( px - ts.width/2, py + ts.height/2 );
        
        cv::rectangle(img, 
            cv::Point2d( px - ts.width, py - ts.height), 
            cv::Point2d( px + ts.width, py + ts.height),
            CV_RGB(0,0,0), CV_FILLED
        );
        
        cv::rectangle(img, 
            cv::Point2d( px - ts.width + 1, py - ts.height + 1), 
            cv::Point2d( px + ts.width - 1, py + ts.height - 1 ),
            CV_RGB(255,255,255), 1
        );
        
        to.y -= 9;
        cv::putText(img, buffer, to, 
            cv::FONT_HERSHEY_SIMPLEX, 0.5, 
            CV_RGB(255, 255, 255), 1, 8
        );
        
        sprintf(buffer, "(%.2lf)", quality);
        ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, 0.3, 1, &baseline);
        to.x = px - ts.width/2;
        to.y += 8 + 3;
        
        cv::Scalar col = CV_RGB(255, 0, 0);
        if (quality == 1.0) {
            col = CV_RGB(0, 255, 0);
        }
        if (quality < 1.0) {
            col = CV_RGB(8, 202, 255);
        }
        if (quality < 0.8) {
            col = CV_RGB(255, 255, 0);
        }
        if (quality < 0.5) {
            col = CV_RGB(255, 187, 8);
        } 
        if (quality < 0.2) {
            col = CV_RGB(255, 0, 0);
        }
        
        cv::putText(img, buffer, to, 
            cv::FONT_HERSHEY_SIMPLEX, 0.3, 
            col, 1, 8
        );
        
    }
    
    const cv::Mat& img;
    cv::Mat out_img;
    string ofname;
    bool    lpmm_mode;
    double  pixel_size;
};

#endif
