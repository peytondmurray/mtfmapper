#ifndef MTF_RENDERER_ANNOTATE_H
#define MTF_RENDERER_ANNOTATE_H

#include "mtf_renderer.h"
#include "common_types.h"

class Mtf_renderer_annotate : public Mtf_renderer {
  public:
    Mtf_renderer_annotate(const cv::Mat& in_img, const std::string& fname) 
      : img(in_img), ofname(fname) {
      
          out_img = img.clone();
    }
    
    void render(const vector<Block>& blocks) {
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0) {
                    Point cent = blocks[i].get_edge_centroid(k);
                    write_number(out_img, lrint(cent.x), lrint(cent.y), val);
                }
            }
        }    
        
        imwrite(ofname, out_img);
    }
    
    void write_number(cv::Mat& img, int px, int py, double val) {
        char buffer[10];
        sprintf(buffer, "%.2lf", val);
        
        int baseline = 0;
        
        cv::Size ts = cv::getTextSize(buffer, cv::FONT_HERSHEY_SIMPLEX, 0.5, 1, &baseline);
        
        ts.width -= 4; // tweak box size slightly
        
        cv::Point to( px - ts.width/2, py + ts.height/2 );
        
        cv::rectangle(img, 
            cv::Point( px - ts.width, py - ts.height), 
            cv::Point( px + ts.width, py + ts.height),
            CV_RGB(0,0,0), CV_FILLED
        );
        
        cv::rectangle(img, 
            cv::Point( px - ts.width + 1, py - ts.height + 1), 
            cv::Point( px + ts.width - 1, py + ts.height - 1 ),
            CV_RGB(65535,65535,65535), 1
        );
        
        cv::putText(img, buffer, to, 
            cv::FONT_HERSHEY_SIMPLEX, 0.5, 
            CV_RGB(65535, 65535, 65535), 1, 8
        );
        
    }
    
    const cv::Mat& img;
    cv::Mat out_img;
    string ofname;
};

#endif
