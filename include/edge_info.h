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
#ifndef EDGE_INFO_H
#define EDGE_INFO_H

#include <string>
using std::string;

#include <opencv2/core/core.hpp>

class Edge_info {
  public:
    Edge_info(void) {
    }
    
    // serialize multiple edges into a string
    static string serialize(
        const vector<cv::Point2d>& in_centroid,
        const vector<double>& in_angle,
        const vector<double>& in_mtf50,
        const vector<double>& /*in_quality*/,
        const vector<std::shared_ptr<vector<double>>>& in_sfr,
        const vector<std::shared_ptr<vector<double>>>& in_esf,
        const vector<cv::Point2d> in_snr,
        const vector<cv::Point2d> in_chromatic_aberration,
        const vector<bool> valid_edge
    ) {
        // there is no need to keep data grouped by block at this point, so each record is independent
        string s = "";
        vector<char> buffer(20*1024);
        for (size_t k=0; k < in_centroid.size(); k++) {
            if (!valid_edge[k]) continue;
            
            char* buf = buffer.data();
            
            buf += sprintf(buf, "%.3lf %.3lf %lf %lf %.2lf %.3lf %.3lf %.3lf ", 
                in_centroid[k].x, in_centroid[k].y, 
                in_angle[k], // TODO: could be relative angle ...
                in_mtf50[k], // MTF-XX, in c/p
                in_snr[k].x, // mean CNR
                in_snr[k].y,  // oversampling factor
                in_chromatic_aberration[k].x, in_chromatic_aberration[k].y // CA red-green, blue-green
            );
            
            buf += sprintf(buf, "%d ", (int)in_sfr[k]->size());
            for (size_t j=0; j < in_sfr[k]->size(); j++) {
                buf += sprintf(buf, "%lf ", (*in_sfr[k])[j]);
            }
            
            buf += sprintf(buf, "%d ", (int)in_esf[k]->size());
            for (size_t j=0; j < in_esf[k]->size(); j++) {
                buf += sprintf(buf, "%.3lf ", (*in_esf[k])[j]);
            }
            
            buf += sprintf(buf, "\n");
            
            s += string(buffer.data());
        }
        return s;
    }
    
    // deserialize a single edge from string s
    static Edge_info deserialize(string s) {
        Edge_info b;
        
        // TODO: we expect s to contain only one row (one entry). Some checking would be good
        
        char* s_buf = const_cast<char*>(s.data());
        FILE* fin = fmemopen(s_buf, s.length(), "r");
        
        int nread = 0;
        
        nread += fscanf(fin, "%lf %lf", &b.centroid.x, &b.centroid.y);
        nread += fscanf(fin, "%lf", &b.angle);
        nread += fscanf(fin, "%lf", &b.mtf50);
        nread += fscanf(fin, "%lf %lf", &b.snr.x, &b.snr.y);
        nread += fscanf(fin, "%lf %lf", &b.chromatic_aberration.x, &b.chromatic_aberration.y);
        
        int nsfr = 0;
        nread += fscanf(fin, "%d", &nsfr);
        
        b.sfr = std::shared_ptr<vector<double>>(new vector<double>(nsfr, 0.0));
        for (int j=0; j < nsfr; j++) {
            nread += fscanf(fin, "%lf", &(*b.sfr)[j]);
        }
        
        int nesf = 0;
        nread += fscanf(fin, "%d", &nesf);
        
        b.esf = std::shared_ptr<vector<double>>(new vector<double>(nesf, 0.0));
        for (int j=0; j < nesf; j++) {
            nread += fscanf(fin, "%lf", &(*b.esf)[j]);
        }
        
        fclose(fin);
        
        return b;
    }
    
    void set_mtf_contrast(double c) {
        mtf_contrast = c;
    }
    
    void set_pixel_pitch(double p) {
        pixel_pitch = p;
    }
    
    cv::Point2d centroid;
    double angle;
    double mtf50;
    double quality;
    std::shared_ptr<vector<double>> sfr;
    std::shared_ptr<vector<double>> esf;
    cv::Point2d snr;
    cv::Point2d chromatic_aberration;
    double pixel_pitch = 1.0;
    double mtf_contrast = 0.5;
    
  private:
};

#endif
