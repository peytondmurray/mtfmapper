#ifndef SVG_PAGE_H
#define SVG_PAGE_H

#include <stdio.h>

#include <cmath>
#include <string>
using std::string;

#include <cv.h>
using namespace cv;

typedef cv::Point_<double> dPoint;
typedef cv::Point_<int> iPoint;

class Svg_page {
  public:
    Svg_page(const string& page_spec, const string& fname) 
      : fname(fname)  {
        
        if (page_spec == "A4" || page_spec == "a4") {
            width_mm = 210;
            height_mm = 297;
        } 
        if (page_spec == "A3" || page_spec == "a3") {
            width_mm = 297;
            height_mm = 420;
        }
        if (page_spec == "A2" || page_spec == "a2") {
            width_mm = 420;
            height_mm = 594;
        }
        
        width  = width_mm * 100;
        height = height_mm * 100;
        
        fout = fopen(fname.c_str(), "wt");
        emit_header();
        
        // set default style
        style="fill:black;stroke:black;stroke-width:1";
    }
    
    ~Svg_page(void) {
        fprintf(fout, "\n\n</svg>\n");
        fclose(fout);
    }
    
    
    
    virtual void render(void) = 0;
    
    
  protected:  
    void emit_header(void) {
        fprintf(fout, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fout, "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.0\" ");
        fprintf(fout, "width=\"%Zdmm\" height=\"%Zdmm\" viewBox=\"0 0 %Zd %Zd\">\n", width_mm, height_mm, width, height);
        
        // draw a white rectangle to cover the background
        fprintf(fout, "  <rect x=\"%Zd\" y=\"%Zd\" width=\"%Zd\" height=\"%Zd\" style=\"fill:white;stroke:none\"/>\n", 
            size_t(0), size_t(0), width, height
        );

        // draw a black border rectangle
        fprintf(fout, "  <rect x=\"%Zd\" y=\"%Zd\" width=\"%Zd\" height=\"%Zd\" style=\"fill:none;stroke:black;stroke-width:1\"/>\n", 
            size_t(4), size_t(4), width-4, height-4
        );
    }
    
    void square(size_t tlx, size_t tly, size_t width) {
        fprintf(fout, "  <rect x=\"%Zd\" y=\"%Zd\" width=\"%Zd\" height=\"%Zd\" style=\"%s\"/>", tlx, tly, width, width, style.c_str());
    }

    void rect(size_t tlx, size_t tly, size_t width, size_t height) {
        fprintf(fout, "  <rect x=\"%Zd\" y=\"%Zd\" width=\"%Zd\" height=\"%Zd\" style=\"%s\"/>", tlx, tly, width, height, style.c_str());
    }
    
    virtual iPoint project(double x, double y) {
        
        x = floor(x*width);
        y = floor(y*height);
        
        return iPoint(int(x), int(y));
    }
    
    
    void rotated_square(double tlx, double tly, double bwidth, double angle) {
        iPoint p = project(tlx + bwidth*cos(angle), tly + bwidth*sin(angle));
        fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+M_PI/2.0), tly + bwidth*sin(angle+M_PI/2.0));
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+M_PI), tly + bwidth*sin(angle+M_PI));
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(tlx + bwidth*cos(angle+1.5*M_PI), tly + bwidth*sin(angle+1.5*M_PI));
        fprintf(fout, "%d,%d\" style=\"%s\"/>", p.x, p.y, style.c_str());
    }
    
    
    string fname;
    string style;
    FILE*  fout;
    
    size_t width_mm;   // in mm
    size_t height_mm;  // in mm
    size_t width;   // in svg units
    size_t height;  // in svg units
    
};

#endif