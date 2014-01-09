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
    Edge_record(void) : pooled(false) {
    }

    // constructor for merging two data sets
    Edge_record(const Edge_record& a, const Edge_record& b) {
        double cx = a.centroid.x;
        double cy = a.centroid.y;
        if (a.orientation == HORIZONTAL) {
            double temp = cx;
            cx = cy;
            cy = temp;
        }
        for (size_t i=0; i < a.points.size(); i++) {
            points.push_back(
                make_pair(
                    a.points[i].first - cx, 
                    a.points[i].second - cy
                )
            );
            weights.push_back(a.weights[i]);
        }

        cx = b.centroid.x;
        cy = b.centroid.y;
        if (b.orientation == HORIZONTAL) {
            double temp = cx;
            cx = cy;
            cy = temp;
        }
        for (size_t i=0; i < b.points.size(); i++) {
            points.push_back(
                make_pair(
                    b.points[i].first - cx, 
                    b.points[i].second - cy
                )
            );
            weights.push_back(b.weights[i]);
        }

        rsq = lsq_fit(points, weights, slope, offset);
        //printf("got %d points in merged sample, rsq=%lf\n", (int)points.size(), rsq);

        const double il_thresh = rsq*1.1;
        int omitted = 0;
        for (size_t i=0; i < points.size(); i++) {
            double ey = fabs(offset + slope*points[i].first - points[i].second);
            if (ey > il_thresh) {
                weights[i] = 0;
                omitted++;
            }
        }
        if (omitted < 0.6*points.size()) {
            rsq = lsq_fit(points, weights, slope, offset);
        }

        //printf("slope=%lf, offset=%lf\n", slope, offset);
        //printf("*got %d points in merged sample, rsq=%lf\n", (int)points.size() - omitted, rsq);
    }

    typedef enum {HORIZONTAL, VERTICAL} orientation_t;

    bool compatible(const Edge_record& b) {
        double Z = (slope - b.slope)/sqrt(sB*sB + b.sB*b.sB);

        return fabs(Z) < 1.66; // ~90% confidence interval on t-distribution with ~80 dof
    }

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
        for (map<int, double>::const_iterator it=weight->begin(); it != weight->end(); ++it) {
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

        for (map<int, double>::const_iterator it=weight->begin(); it != weight->end(); ++it) {
            if (it->second >= weight_thresh) {
                points.push_back(make_pair(it->first, (*mean)[it->first] / it->second));
                weights.push_back(it->second);
            }
        }

        double ratio = double(mean->size())/double(other_mean->size());
        if (ratio < 1.6) {
            for (map<int, double>::const_iterator it=other_weight->begin(); it != other_weight->end(); ++it) {
                if (it->second >= weight_thresh) {
                    points.push_back(make_pair((*other_mean)[it->first] / it->second, it->first));
                    weights.push_back(it->second/ratio);
                }
            }
        }

        
        rsq = lsq_fit(points, weights, slope, offset);

        const double il_thresh = 2;
        centroid = Point(0,0);
        int c_count = 0;
        for (size_t i=0; i < points.size(); i++) {
            double ey = fabs(offset + slope*points[i].first - points[i].second);
            if (ey > il_thresh) {
                weights[i] /= 100000;
            } else {
                centroid.x += points[i].first;
                centroid.y += points[i].second;
                c_count++;
            }
        }
        rsq = lsq_fit(points, weights, slope, offset);
        //printf("n=%d, slope=%lf, offset=%lf, rsq=%lf\n", (int)points.size(), slope, offset, rsq);
        centroid.x /= double(c_count);
        centroid.y /= double(c_count);
        mse = 0;
        double ss_x = 0;
        int e_count = 0;
        for (size_t i=0; i < points.size(); i++) {
            double ey = fabs(offset + slope*points[i].first - points[i].second);
            if (ey <= il_thresh) {
                mse += ey*ey;
                ss_x += (points[i].first - centroid.x)*(points[i].first - centroid.x);
                e_count++;
            }
        }
        mse /= e_count - 2;
        sB = sqrt(mse/ss_x);

        dx = points[0].first - points[points.size()-1].first;
        dy = points[0].second - points[points.size()-1].second;

        if (orientation == VERTICAL) {
            double temp = dx;
            dx = dy;
            dy = temp;
        } else {
            double temp = centroid.x;
            centroid.x = centroid.y;
            centroid.y = temp;
        }

        set_angle_from_slope(slope);
        
        //printf("slope estimate is %lf, %lf degrees, rsq = %lf, offset=%lf\n", slope, angle/M_PI*180, rsq, offset);
    }

    inline orientation_t get_orientation(void) const {
        return orientation;
    }

    inline bool is_pooled(void) {
        return pooled;
    }

    inline void set_angle_from_slope(double slope, bool is_pooled=false) {
        if (orientation == VERTICAL) {
            slope = 1.0/slope;
        }

        if (dx > 0) {
            angle = -atan(slope);
        } else {
            if (dy >= 0) {
                angle = -atan(slope) + M_PI;
            } else {
                angle = -atan(slope) - M_PI;
            }
        }
        if (angle < -M_PI) {
            angle += M_PI;
        }
        if (angle > M_PI) {
            angle -= M_PI;
        }

        pooled = is_pooled;
    }

    double slope;
    double offset;
    double angle;
    double rsq;
    Point  centroid;

  private:

    double lsq_fit(const vector< pair<double, double> >& points,
        const vector< double >& weights, double& slope, double& offset) {

        double rsq = 0;
    
        int n = weights.size();
        
        double sx  = 0;
        double sy  = 0;
        double ss  = 0;
        
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

    double mse;
    double sB; // standard error in slope estimate

    double dx;
    double dy;

    bool pooled;
};

#endif 


