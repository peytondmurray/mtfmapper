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
#ifndef ENTRY_VIEW_H
#define ENTRY_VIEW_H

#include "include/logger.h"
#include "include/edge_info.h"

#include <string>
using std::string;

class Plot_details {
  public:
    Plot_details(string title, string xlabel, string ylabel)
    : title(title), xlabel(xlabel), ylabel(ylabel) {}
    
    string title;
    string xlabel;
    string ylabel;
  private:
};

class Entry_view {
  public:
    typedef enum {
        VIEW_SFR = 0,
        VIEW_ESF,
        VIEW_LSF
    } view_t;
    
    Entry_view(bool lp_mm_mode=false) 
    : lp_mm_mode(lp_mm_mode), default_pixel_pitch(1000.0), view(VIEW_SFR), previous_view(VIEW_SFR) {
        plot_details = { 
          {"SFR / MTF curve", "Frequency (c/p)", "Contrast"},
          {"SFR / MTF curve", "Frequency (lp/mm)", "Contrast"},
          {"ESF profile", "Distance (pixels)", "Intensity"},
          {"LSF profile ", "Distance (pixels)", "Intensity"}
        };
    }
    
    void set_view(view_t nview) {
        view = nview;
    }
    
    view_t get_view(void) const {
        return view;
    }
    
    void push_view(void) {
        previous_view = view;
    }
    
    void pop_view(void) {
        view = previous_view;
    }
    
    const string& title(void) const {
        return plot_details[view].title;
    }
    
    const string& xlabel(void) const {
        return plot_details[view].xlabel;
    }
    
    const string& ylabel(void) const {
        return plot_details[view].ylabel;
    }
    
    const string mtf_xx(const Sfr_entry& entry) const {
        char buffer[200];
        if (lp_mm_mode) {
            sprintf(buffer, "MTF%2d=%.2lf lp/mm", (int)rint(100*entry.info.mtf_contrast), entry.info.mtf50*frequency_scale(entry));
        } else {
            sprintf(buffer, "MTF%2d=%.3lf c/p", (int)rint(100*entry.info.mtf_contrast), entry.info.mtf50);
        }
        return string(buffer);
    }
    
    
    string x_format(double val) const {
        char buffer[200];
        switch(view) {
        case VIEW_SFR: 
            if (lp_mm_mode) {
                sprintf(buffer, "frequency: %5.1lf ", val);
            } else {
                sprintf(buffer, "frequency: %5.3lf ", val); 
            }
            break;
        case VIEW_ESF: 
        case VIEW_LSF: sprintf(buffer, "distance: %5.1lf ", val); break;
        }
        
        return string(buffer);
    }
    
    string y_format(double val) const {
        char buffer[200];
        switch(view) {
        case VIEW_SFR: sprintf(buffer, "contrast: %5.3lf ", val); break;
        case VIEW_ESF: 
        case VIEW_LSF: sprintf(buffer, "intensity: %5.0lf ", val); break;
        }
        
        return string(buffer);
    }
    
    void update(const vector<Sfr_entry>& entries, vector<QLineSeries*>& series,
        QChart& chart, QValueAxis& x_axis, QValueAxis& y_axis) {
        
        chart.removeAllSeries();
        series.clear();
        for (const auto& e : entries) {
            series.push_back(new QLineSeries());
            populate(e, *series.back());
        }
        
        axis_setup(series, x_axis, y_axis);
        
        // note the order: first gather min/max, then add series to chart, then attach axis
        for (auto& s : series) {
            chart.addSeries(s);
            s->attachAxis(&x_axis);
            s->attachAxis(&y_axis);
        }
    }
    
    bool set_lp_mm_mode(bool new_lp_mm_mode) {
        bool changed = lp_mm_mode != new_lp_mm_mode;
        lp_mm_mode = new_lp_mm_mode;
        return changed;
    }
    
    bool set_default_pixel_pitch(double new_pixel_pitch) {
        bool changed = default_pixel_pitch != new_pixel_pitch;
        default_pixel_pitch = new_pixel_pitch;
        return changed;
    }
    
  private:
    void axis_setup(vector<QLineSeries*>& series, QValueAxis& x_axis, QValueAxis& y_axis) {
        double min_x = 1e50;
        double max_x = -1e50;
        double min_y = 1e50;
        double max_y = -1e50;
        for (auto s : series) {
            QVector<QPointF> points = s->pointsVector();
            for (const QPointF& p : points) {
                max_x = std::max(max_x, p.x());
                min_x = std::min(min_x, p.x());
                max_y = std::max(max_y, p.y());
                min_y = std::min(min_y, p.y());
            }
        }
        x_axis.setRange(min_x, max_x);
        x_axis.setTitleText(xlabel().c_str());
        y_axis.setTitleText(ylabel().c_str());
        
        if (view == VIEW_SFR) {
            x_axis.setTickCount(11);
            if (lp_mm_mode) {
                x_axis.setLabelFormat("%3.0f");
            } else {
                x_axis.setLabelFormat("%3.1f");
            }
            
            double roundup_max = multiple(max_y);
            y_axis.setRange(0, roundup_max);
            y_axis.setTickCount(roundup_max*5 + 1);
            y_axis.setLabelFormat("%3.1f");
        } else {
            x_axis.setTickCount(17);
            x_axis.setLabelFormat("%2.0f");
            
            if (view == VIEW_ESF) {
                y_axis.setRange(0, max_y*1.05);
            } else {
                y_axis.setRange(min_y < 0 ? 1.05*min_y : min_y, max_y*1.05);
            }
            
            y_axis.setTickCount(10);
            y_axis.setLabelFormat("%5.0f");
        }
    }
  
    void populate(const Sfr_entry& entry, QLineSeries& s) {
        QVector<QPointF> points;
        
        if (view == VIEW_SFR) {
            // least-squares fit (inverse of design matrix) of a cubic polynomial through 4 points [-1..2]/64.0
            constexpr  double sfr_cubic_weights[4][4] = {
              { 0.00000000000000e+00,   1.00000000000000e+00,   0.00000000000000e+00,   0.00000000000000e+00},
              {-2.13333333333333e+01,  -3.20000000000000e+01,   6.40000000000000e+01,  -1.06666666666667e+01},
              { 2.04800000000000e+03,  -4.09600000000000e+03,   2.04800000000000e+03,   0.00000000000000e+00},
              {-4.36906666666667e+04,   1.31072000000000e+05,  -1.31072000000000e+05,   4.36906666666667e+04}
            };
        
            vector<double>& sfr = *entry.info.sfr;
            for (size_t i=0; i < sfr.size(); i++) {
                double coef[4] = {0,0,0,0};
                
                for (int si = int(i) - 1, ri = 0; si <= int(i) + 2; si++, ri++) {
                    int ei = si < 0 ? 0 : si;
                    double eiv = 0;
                    if (ei > (int)sfr.size() - 1) {
                        eiv = (sfr[sfr.size() - 1] - sfr[sfr.size() - 2]) * (ei - (sfr.size() - 1)) + sfr[sfr.size() - 1]; // extend last point linearly
                    } else {
                        eiv = sfr[ei];
                    }

                    for (int c = 0; c < 4; c++) {
                        coef[c] += eiv * sfr_cubic_weights[c][ri];
                    }
                }
                // now we can evaluate points at position x in [0,1) using the fitted polynomial
                double freq_scale = frequency_scale(entry);
                for (int xi=0; xi < 20; xi++) {
                    double dx = xi*(1.0/(20.0*64.0));
                    double iy = coef[0] + coef[1]*dx + coef[2]*dx*dx + coef[3]*dx*dx*dx;
                    
                    points.push_back(QPointF(freq_scale*(i*(1.0/64.0) + dx), iy));
                }
            }
        }
        
        
        if (view == VIEW_ESF || view == VIEW_LSF) {
            bool reverse = false;
            vector<double>& esf = *entry.info.esf;
            size_t left_lower = 0;
            for (size_t i=0; i < 2*8; i++) {
                left_lower += esf[i] < esf[esf.size()-1 - i];
            }
            reverse = left_lower < 8;
        
            if (view == VIEW_ESF) {
                vector<double>& esf = *entry.info.esf;
                if (reverse) {
                    for (size_t i=0; i < esf.size(); i++) {
                        points.push_back(QPointF(float(i)*0.125 - 16.0 - 0.25, esf[esf.size()-1 - i]));
                    }
                } else {
                    for (size_t i=0; i < esf.size(); i++) {
                        points.push_back(QPointF(float(i)*0.125 - 16.0 + 3*0.125, esf[i]));
                    }
                }
            }
            
            if (view == VIEW_LSF) {
                vector<double>& esf = *entry.info.esf;
                if (reverse) {
                    for (size_t i=1; i < esf.size() - 1; i++) {
                        points.push_back(QPointF(float(i)*0.125 - 16.0 - 0.25, (esf[esf.size()-1 - (i+1)] - esf[esf.size()-1 - (i-1)])/0.25));
                    }
                } else {
                    for (size_t i=1; i < esf.size() - 1; i++) {
                        points.push_back(QPointF(float(i)*0.125 - 16.0 + 3*0.125, (esf[i+1] - esf[i-1])/0.25));
                    }
                }
            }
        }
        
        s.replace(points);
    }
    
    double multiple(double x) const {
        double rval = 0;
        while (rval < (x+5e-4)) {
            rval += 0.2;
        }
        return rval;
    }
    
    double frequency_scale(const Sfr_entry& entry) const {
        double freq_scale = 1.0;
        if (lp_mm_mode) {
            if (fabs(entry.info.pixel_pitch - 1000.0) < 1e-6) {
                freq_scale = 1000.0/default_pixel_pitch;
            } else {
                freq_scale = 1000.0/entry.info.pixel_pitch;
            }
        }
        return freq_scale;
    }
    
    bool lp_mm_mode;
    double default_pixel_pitch;
    vector<Plot_details> plot_details;
    view_t view;
    view_t previous_view;
};


#endif



