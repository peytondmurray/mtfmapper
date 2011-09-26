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
#ifndef MTF_RENDERER_STATS_H
#define MTF_RENDERER_STATS_H

#include "mtf_renderer.h"
#include "common_types.h"

#include <vector>
using std::vector;

#include <cmath>

class Mtf_renderer_stats : public Mtf_renderer {
  public:
    Mtf_renderer_stats(void)  {
      
    }
    
    void render(const vector<Block>& blocks) {
        
        if (blocks.size() == 0) {
            return;
        }
    
        vector<double> unfiltered;
        vector<double> filtered;
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
            
                double val = blocks[i].get_mtf50_value(k);
                unfiltered.push_back(val);
                
                if (blocks[i].get_quality(k) >= 0.5) {
                    filtered.push_back(val);
                }
            }
        }
        
        if (unfiltered.size() < 2) {
            return;
        }
        
        sort(filtered.begin(), filtered.end());
        sort(unfiltered.begin(), unfiltered.end());
        
        printf("    Quantiles ->                   %4d%% %4d%% %4d%% %4d%% %4d%%\n", 5, 25, 50, 75, 95);
        printf("Statistics on all edges:           %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf  (total=%u)\n", 
            quantile(unfiltered, 0.05),
            quantile(unfiltered, 0.25),
            quantile(unfiltered, 0.5),
            quantile(unfiltered, 0.75),
            quantile(unfiltered, 0.95),
            (unsigned int)unfiltered.size()
        );
        
        if (filtered.size() < 2) {
            return;
        }
        
        printf("Statistics on all filtered edges : %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf  (total=%u)\n", 
            quantile(filtered, 0.05),
            quantile(filtered, 0.25),
            quantile(filtered, 0.5),
            quantile(filtered, 0.75),
            quantile(filtered, 0.95),
            (unsigned int)filtered.size()
        );
    }
    
  private:
    double quantile(const vector<double>& d, double q) {
        size_t idx = (int)floor(d.size()*q);
        return d[idx];
    }
};

#endif
