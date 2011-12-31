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
#ifndef EDGE_RECORD_H
#define EDGE_RECORD_H

#include "common_types.h"


class Edge_record {
  public:
    Edge_record(void) {
    }

    typedef enum {HORIZONTAL, VERTICAL} orientation_t;

    void add_point(int x, int y, double gx, double gy) {
        if (col_mean.find(x) == col_mean.end()) {
            col_mean[x] = y * gy;
            col_weight[x] = gy;
        } else {
            col_mean[x] += y * gy;
            col_weight[x] += gy;
        }

        if (row_mean.find(y) == row_mean.end()) {
            row_mean[y] = x * gx;
            row_weight[y] = gx;
        } else {
            row_mean[y] += x * gx;
            row_weight[y] += gx;
        }
    }

    void reduce(void) { // compute orientation, and remove weak points

        map<int, double>* mean = 0;
        map<int, double>* weight = 0;
        map<int, double>* other_mean = 0;
        map<int, double>* other_weight = 0;
        if (col_mean.size() > row_mean.size()) {
            mean = &col_mean;
            weight = &col_weight;
            other_mean = &row_mean;
            other_weight = &row_weight;
            orientation = VERTICAL;
        } else {
            mean = &row_mean;
            weight = &row_weight;
            other_mean = &col_mean;
            other_weight = &col_weight;
            orientation = HORIZONTAL;
        }

        vector<double> sweight;
        for (map<int, double>::const_iterator it=weight->begin(); it != weight->end(); it++) {
            if (it->second > 0) {
                sweight.push_back(it->second);
            }
        }
        sort(sweight.begin(), sweight.end());

        if (sweight.size() == 0) {
            printf("Major error: no edge points found!\n");
            return;
        }

        double weight_thresh = 0;
        if (sweight.size() > 5) {
            weight_thresh = sweight[int(0.2*sweight.size())];
        } else {
            weight_thresh = sweight[0];
        }

        for (map<int, double>::const_iterator it=weight->begin(); it != weight->end(); it++) {
            if (it->second >= weight_thresh) {
                points.push_back(make_pair<double, double>(it->first, (*mean)[it->first] / it->second));
                weights.push_back(it->second);
            }
        }

        double ratio = double(mean->size())/double(other_mean->size());
        if (ratio < 1.4) {
            for (map<int, double>::const_iterator it=other_weight->begin(); it != other_weight->end(); it++) {
                if (it->second >= weight_thresh) {
                    points.push_back(make_pair<double, double>((*other_mean)[it->first] / it->second, it->first));
                    weights.push_back(it->second/ratio);
                }
            }
        }

        //printf("got %d points after filtering\n", (int)weights.size());
        
        rsq = lsq_fit(points, weights, slope, offset);
        //printf("slope=%lf, offset=%lf\n", slope, offset);

        const double il_thresh = 2;
        for (size_t i=0; i < points.size(); i++) {
            double ey = fabs(offset + slope*points[i].first - points[i].second);
            if (ey > il_thresh) {
                weights[i] /= 100000;
            }
        }
        rsq = lsq_fit(points, weights, slope, offset);
        printf("slope=%lf, offset=%lf, rsq=%lf\n", slope, offset, rsq);


        double dx = points[0].first - points[points.size()-1].first;
        double dy = points[0].second - points[points.size()-1].second;
        double angle_offset = 0;

        if (orientation == VERTICAL) {
            dx = points[0].second - points[points.size()-1].second;
            dy = points[0].first - points[points.size()-1].first;
            slope = 1.0/slope;
        }

        if (dx > 0) {
            angle = -atan(slope) + angle_offset;
        } else {
            if (dy >= 0) {
                angle = -atan(slope) + M_PI + angle_offset;
            } else {
                angle = -atan(slope) - M_PI + angle_offset;
            }
        }
        if (angle < -M_PI) {
            angle += M_PI;
        }
        if (angle > M_PI) {
            angle -= M_PI;
        }

        //printf("slope estimate is %lf, %lf degrees, rsq = %lf, offset=%lf\n", slope, angle/M_PI*180, rsq, offset);
    }

    double slope;
    double offset;
    double angle;
    double rsq;

  private:

    double lsq_fit(const vector< pair<double, double> >& points,
        const vector< double >& weights, double& slope, double& offset) {

        double rsq = 0;
    
        int n = weights.size();
        
        double sx  = 0;
        double sy  = 0;
        double ss = n;
        
        ss = 0;
        for (int i=0; i < n; i++) {
            double weight = SQR(weights[i]);
            ss += weight;
            sx += points[i].first * weight;
            sy += points[i].second * weight;
        }
        double sxoss = sx / ss;
        
        double st2 = 0;
        double b = 0;
        for (int i=0; i < n; i++) {
            double t = (points[i].first - sxoss) * weights[i];
            st2 += t*t;
            b += t*points[i].second * weights[i];
        }
        b /= st2;
        double a = (sy - sx*b)/ss;
        offset = a;
        slope = b;
        for (int i=0; i < n; i++) {
            double r = (points[i].first*slope + offset) - points[i].second;
            rsq += fabs(r); // m-estimate of goodness-of-fit 
        }
        return rsq/double(n);
    }

    orientation_t orientation;

    vector< pair<double, double> > points;
    vector< double > weights;

    map<int, double> row_mean;
    map<int, double> row_weight;
    map<int, double> col_mean;
    map<int, double> col_weight;

};

#endif 


