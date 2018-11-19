#ifndef PEAK_DETECTOR_H
#define PEAK_DETECTOR_H

#include <cmath>
#include <vector>
using std::vector;

#include <queue>
using std::priority_queue;

#include <stdlib.h>

#define DBG(x) 

class Peak_detector {
  public:
    Peak_detector(const vector<double>& in_data, size_t in_bins, double in_min_data=-M_PI, double in_max_data=M_PI)
      : nbins(in_bins), data(in_data), counts(in_bins,0), min_data(in_min_data), max_data(in_max_data) {
          
        range = max_data - min_data;
        
        // build histogram
        DBG(FILE* rawpts = fopen("raw_points.txt", "wt");)
        for (size_t i=0; i < data.size(); i++) {
            DBG(fprintf(rawpts, "%lf\n", data[i]);)
            size_t idx = lrint((nbins*(data[i] - min_data)/range));
            if (idx > counts.size() - 1) {
                idx = counts.size() - 1;
            }
            counts[idx]++;
        }
        DBG(fclose(rawpts);)
  
        // perform non-minimum suppression
        vector<size_t> nm(nbins, 0);
        const int delta = 3;
        for (size_t i=0; i < nbins; i++) {
            int start = int(i) - delta;
            int end   = int(i) + delta;
            
            size_t max_count = 0;
            for (int j=start; j < end; j++) {
                int eidx = j;
                if (eidx < 0) eidx = (eidx + nbins) % nbins;
                eidx = eidx % nbins;
                
                if (counts[eidx] > max_count) {
                    max_count = counts[eidx];
                }
            }
            if (counts[i] < max_count) {
                nm[i] = 0;
            } else {
                nm[i] = counts[i];
            }
        }
        
        // now extract the peaks, noting both their value, and their counts
        DBG(
          FILE* fout = fopen("peaks.txt", "wt");
        )
        for (size_t i=0; i < nbins; i++) {
            DBG(
              fprintf(fout, "%Zd %lf %Zd %Zd\n", i, i*range/double(nbins), nm[i], counts[i]);
            )
            if (nm[i] > 0) {
                pq.push( pair<size_t, double>(counts[i], min_data + double(range*i)/double(nbins)) );
            }
        } 
        DBG(
          fclose(fout);
        )
    }
      
    void select_best_n(vector<double>& best, size_t in_n, double il_thresh=0.1) { // 1% of the range
      
        best.clear();
        for (size_t k=0; k < in_n && !pq.empty(); k++) {
            best.push_back(pq.top().second);
            pq.pop();
        }
        
        // now go back to the data, and refine the values of the best ones
        vector< vector<double> > members(4);
        for (size_t i=0; i < data.size(); i++) {
            for (size_t k=0; k < best.size(); k++) {
                if (angular_diff(data[i], best[k]) < il_thresh * range) {
                    if ( angular_diff(best[k], min_data) < il_thresh*range && data[i] < 0) { // close to minimum
                        // shift the origin if the peak is anywhere near -pi 
                        members[k].push_back(data[i] + 2*M_PI);
                    } else {
                        members[k].push_back(data[i]);
                    }
                    
                }
            }
        }
        
        // compute trimmed mean, discarding <5 and >95 percentiles
        for (size_t k=0; k < best.size(); k++) {
            sort(members[k].begin(), members[k].end());
            double best_sum = 0.0;
            double best_count = 0.0;
            size_t lower = members[k].size() / 20;
            size_t upper = (members[k].size() * 19) / 20;
            for (size_t j=lower; j < upper; j++) {
                double weight = 1.0;
                best_sum += members[k][j] * weight;
                best_count += weight;
            }
            best[k] = best_sum / best_count;
        }
    }
      
    static double angular_diff(double a, double b) {
        return acos(cos(a)*cos(b) + sin(a)*sin(b));
    }
      
  private:
    size_t nbins;
    const vector<double>& data;
    vector <size_t> counts;
    double min_data;
    double max_data;
    double range;
    
    priority_queue < pair<size_t, double> > pq;
};

#endif
