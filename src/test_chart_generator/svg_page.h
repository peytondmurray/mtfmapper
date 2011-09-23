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
#ifndef SVG_PAGE_H
#define SVG_PAGE_H

#include <stdio.h>

#include <cmath>
#include <string>
using std::string;

#ifdef _MSC_VER
	typedef unsigned short int uint16_t;

	#define M_PI 3.14159265358979
	#define lrint(x) ( (x < 0) ? int(floor(x-0.5)) : int(floor(x+0.5)) )
#endif

#include "opencv/cv.h"
using namespace cv;

typedef cv::Point_<double> dPoint;
typedef cv::Point_<int> iPoint;

class Svg_page {
  public:
    Svg_page(const string& page_spec, const string& fname) 
      : fname(fname), page_size(page_spec)  {
        
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
        if (page_spec == "A1" || page_spec == "a1") {
            width_mm = 594;
            height_mm = 841;
        }
        if (page_spec == "A0" || page_spec == "a0") {
            width_mm = 841;
            height_mm = 1189;
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
        fprintf(fout, "width=\"%dmm\" height=\"%dmm\" viewBox=\"0 0 %d %d\">\n", int(width_mm), int(height_mm), int(width), int(height));
        
        // draw a white rectangle to cover the background
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:white;stroke:none\"/>\n", 
            int(0), int(0), int(width), int(height)
        );

        // draw a black border rectangle
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:none;stroke:black;stroke-width:1\"/>\n", 
            int(4), int(4), int(width-4), int(height-4)
        );
    }
    
    void square(size_t tlx, size_t tly, size_t width) {
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"%s\"/>", int(tlx), int(tly), int(width), int(width), style.c_str());
    }

    void rect(size_t tlx, size_t tly, size_t width, size_t height) {
        fprintf(fout, "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"%s\"/>", int(tlx), int(tly), int(width), int(height), style.c_str());
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
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", p.x, p.y, style.c_str());
    }
    
    void text(const std::string& s, double x, double y, int font_size=20) { // x, y are fractional positions
        fprintf(fout, "  <text x=\"%ld\" y=\"%ld\" font-family=\"Verdana\" font-size=\"%ldmm\" fill=\"black\" > %s </text>\n",
            lrint(x*width), lrint(y*height), lrint(font_size), s.c_str()
        );   
    }
    
    
    string fname;
    string style;
    FILE*  fout;
    
    size_t width_mm;   // in mm
    size_t height_mm;  // in mm
    size_t width;   // in svg units
    size_t height;  // in svg units
    std::string page_size;
    
};

#endif
