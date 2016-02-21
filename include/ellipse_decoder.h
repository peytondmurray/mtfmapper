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
#ifndef ELLIPSE_DECODER_H
#define ELLIPSE_DECODER_H

#include "include/common_types.h"
#include "include/ellipse.h"
#include "include/component_labelling.h"

const int bitreverse4[16] = {
    0, 8, 4, 12, 2, 10, 6, 14,
    1, 9, 5, 13, 3, 11, 7, 15
};

class Ellipse_decoder {
  public:
    
    Ellipse_decoder(const Ellipse_detector& e, const cv::Mat& img, 
        const Point& trans) : origin(e.centroid_x, e.centroid_y), 
        trans(trans), code(-1), valid(false), ratio(e.minor_axis/e.major_axis) {
        
        if (e.solid) {
            valid = true;
            code = 0;
        } else {
            string bitvec;
            extract_bitvector(e, img, bitvec);
            decode_bitvector(bitvec);
        }
    }
    
    void extract_bitvector(const Ellipse_detector&e, 
        const cv::Mat& img, string& bitvec) {
        
        // collect histogram stats inside ellipse
        map<int, int> histo;
        double total = 0;
        for (map<int, scanline>::const_iterator it=e.scanset.begin(); it != e.scanset.end(); it++) {
            int y = it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                int val = img.at<uint16_t>(y, x);
            
                map<int, int>::iterator hi=histo.find(val);
                if (hi == histo.end()) {
                    histo[val] = 0;
                }
                histo[val]++;
                total++;
            }
        }
        double sum = 0;
        // compute Otsu threshold
        for (map<int, int>::iterator it=histo.begin(); it != histo.end(); it++) {
            sum += it->first * it->second;
        }
        double sum_b = 0;
        double w_b = 0;
        double max = 0;
        double thresh1 = 0;
        double thresh2 = 0;
        for (map<int, int>::iterator it=histo.begin(); it != histo.end(); it++) {
            w_b += it->second;
            double w_f = total - w_b;
            if (w_f == 0) break;
            
            sum_b += it->first * it->second;
            double m_b = sum_b / w_b;
            double m_f = (sum - sum_b) / w_f;
            double between = w_b * w_f * (m_b - m_f) * (m_b - m_f);
            if (between >= max) {
                thresh1 = it->first;
                if (between > max) {
                    thresh2 = it->first;
                }
                max = between;
            }
        }
        int otsu = lrint((thresh1 + thresh2) * 0.5);
        
        int first_one = -1;
        int last_one = -1;
        for (double du=-e.major_axis-1; du <= e.major_axis+1; du += 1) {
            int sx = lrint(origin.x - du*trans.x);
            int sy = lrint(origin.y - du*trans.y);
            int bit = img.at<uint16_t>(sy, sx) < otsu ? 1 : 0;
            if (first_one < 0 && bit) {
                first_one = 0;
            }
            if (first_one >= 0) {
                bitvec.push_back(bit + '0');
                first_one++;
            }
            if (bit) {
                last_one = first_one;
            }
        }
        bitvec = bitvec.substr(0, last_one);    
    }
    
    void extract_bitvector(const Ellipse_detector&e, 
        const Component_labeller& cl, string& bitvec) {
        
        int first_one = -1;
        int last_one = -1;
        for (double du=-e.major_axis-1; du <= e.major_axis+1; du += 1) {
            int bit = sample(origin - du*trans, cl);
            if (first_one < 0 && bit) {
                first_one = 0;
            }
            if (first_one >= 0) {
                bitvec.push_back(bit + '0');
                first_one++;
            }
            if (bit) {
                last_one = first_one;
            }
        }
        bitvec = bitvec.substr(0, last_one);    
    }
    
    void decode_bitvector(const string& bitvec) {
        vector<int> bits(4,1);
        // synthesize all the codes, pick one with minimum hamming distance
        int rad = bitvec.size()/2;
        int inner_lim = lrint(0.5*bitvec.size()/1.4);
        int outer_lim = lrint(0.5*bitvec.size()*0.6);
        int segments[7] = {
            0,
            rad - outer_lim,
            inner_lim,
            (int)bitvec.size()/2,
            (int)bitvec.size() - inner_lim,
            (int)bitvec.size() - (rad - outer_lim),
            (int)bitvec.size()
        };
        
        int mindist = bitvec.size() + 1;
        int bestcode = 15; // all ones
        for (int i=0; i < 16; i++) {
            int dist = 0;
            int total = 0;
            for (int j=0; j < 6; j++) {
                char match = '1';
                if (j >= 1 && j < 5) {
                    match = (i & (1 << (3 - (j-1)))) ? '0' : '1';
                }
                for (int k=segments[j]; k < segments[j+1]; k++) {
                    total++;
                    dist += (bitvec[k] == match) ? 0 : 1;
                }
            }
            if (dist < mindist) {
                bestcode = i;
                mindist = dist;
            }
        }
        for (int i=0; i < 4; i++) {
            bits[i] = ((1 << i) & bestcode) ? 1 : 0;
        }
        
        code = bestcode;
        valid = mindist < 0.15*bitvec.size(); 
                    
        printf("%lf (%lf %lf) -> bits=[%d%d%d%d] / %d, mindist=%d/%d, valid=%d, ratio=%lf\n", origin.y, origin.x, origin.y, 
            bits[0], bits[1], bits[2], bits[3], bestcode, mindist, (int)bitvec.size(), valid, ratio);
        
    }
    
    inline int sample(const Point& p, const Component_labeller& cl) {
        return cl(lrint(p.x), lrint(p.y)) <= 0 ? 0 : 1;
    }
    
    void reverse(void) {
        if (code >= 0 && code <= 15) {
            code = bitreverse4[code];
        } else {
            printf("Error: tried to reverse a value outside of the valid range (ellipse decoder): %d\n", code);
        }
    }
    
    static int reverse(int incode) {
        if (incode >= 0 && incode <= 15) {
            return bitreverse4[incode];
        } 
        printf("Error: tried to reverse a value outside of the valid range (ellipse decoder): %d\n", incode);
        return -1;
    }
    
    Point origin;
    Point trans;
    
    int code;
    bool valid;
    
    double ratio;
};

#endif
