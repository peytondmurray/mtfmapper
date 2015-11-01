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
#ifndef AFFT_H
#define AFFT_H

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

template< int n >
class AFFT {
  public:
    AFFT(void) 
    : cstab(2*n), cstable(2*n) {
        
        for (int i=0; i < n; i++) {
            cstab[2*i] = sin(-2*M_PI*i/double(n));
            cstab[2*i+1] = cos(-2*M_PI*i/double(n));
        }
        
        unsigned int nn = n >> 1; // number of complex values
        int j = 1;
        for (int i = 1; i < n; i += 2) {
            if (j > i) {
                bitrev.push_back(make_pair(j - 1, i - 1));
                bitrev.push_back(make_pair(j, i));
            }
            int m = nn;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        };
        
        double twiddlebase = 2*M_PI/double(n/2);
        double y = 0.0;
        for (int i=0; i < n; i++, y += 1.0) {
            double x = y*twiddlebase;
            cstable[2*i] = cos(x);
            cstable[2*i+1] = -sin(x); 
        }
        
        power = -1;
        unsigned int t = n/2;
        while (t > 0) {
            power++;
            t >>= 1;
        }
    }
    
    void forward(double* data) {
        // apply bit-reversing 
        for (size_t j=0; j < bitrev.size(); j++) {
            std::swap(data[bitrev[j].first], data[bitrev[j].second]);
        }
        
        
        // perform first stage where bflys==1 and cos=1, sin=0
        unsigned int node_space = 2; // start at 2 because we interleave real and complex parts
        double* x0p = data;
        for (unsigned int group = 0; group < n/4; group++) {
            double twidf[2] = {*(x0p + node_space), *(x0p + node_space+1)};
            
            *(x0p + node_space)   = *(x0p) - twidf[0];
            *(x0p + node_space+1) = *(x0p+1) - twidf[1];
            *(x0p) += twidf[0];
            *(x0p+1) += twidf[1];

            x0p += 2 + node_space;
        }                       //end group loop
        node_space <<= 1;
        
        unsigned int bflys_per_group = 2; // start at two because first stage is hardcoded
        unsigned int num_groups = n/8; // start at n/8 because first stage is hardcoded
            
        //stage loop
        for (unsigned int stage = 1; stage < power; stage++) {
            double* x0p = data;
            //group loop
            for (unsigned int group = 0; group < num_groups; group++) {

                double* c_ptr = cstable.data();
                //butterfly loop
                for (unsigned int bflys = 0; bflys < bflys_per_group; bflys++) {
                    
                    double twidf[2] = {
                        (*c_ptr) * *(x0p + node_space) - (*(c_ptr+1)) * *(x0p + node_space+1),
                        (*c_ptr) * *(x0p + node_space+1) + (*(c_ptr+1)) * *(x0p + node_space)
                    };
                    
                    *(x0p + node_space)   = *(x0p) - twidf[0];
                    *(x0p + node_space+1) = *(x0p+1) - twidf[1];
                    *(x0p) += twidf[0];
                    *(x0p+1) += twidf[1];

                    x0p += 2;
                    c_ptr += 2*num_groups;
                }                   //end butterfly loop
                x0p += node_space;
            }                       //end group loop
            num_groups >>= 1;
            bflys_per_group = node_space;
            node_space <<= 1;
            
        }                           //end stage loop
    }

    
    void realfft(double *data) {
        forward(data);
        
        double* sp = cstab.data() + 2;
        
        for (int i=1; i < n/2; i++) {
            int k = 2*i;
            int ik = 2*(n/2 - i);
            
            double cr = (data[k] + data[ik]);
            double ci = (data[k+1] - data[ik + 1]);
            
            double br = (data[k] - data[ik]);
            double bi = (data[k+1] + data[ik + 1]);
            
            cr += br*(*sp) + bi*(*(sp+1));
            ci += (*sp)*bi - (*(sp+1))*br;
            data[k] = 0.5*cr;
            data[k+1] = 0.5*ci;
            sp += 2;
        }
        // TODO: deal with data[N/2]
        // TODO: deal with data[0]
    }
    
    
  
    vector<double> cstab;
    vector< pair<int,int> > bitrev;
    vector<double> cstable;
    unsigned int power;
};

#endif
