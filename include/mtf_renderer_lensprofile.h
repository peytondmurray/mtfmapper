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
#ifndef MTF_RENDERER_LENSPROFILE_H
#define MTF_RENDERER_LENSPROFILE_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "include/loess_fit.h"
#include "include/mtf_tables.h"

class Mtf_renderer_lensprofile : public Mtf_renderer {
  public:
    Mtf_renderer_lensprofile(const std::string& wdir, 
        const std::string& prof_fname, 
        const std::string& gnuplot_binary,
        const cv::Mat& img, 
        bool lpmm_mode=false, double pixel_size=1.0) 
      :  wdir(wdir), prname(prof_fname), 
         gnuplot_binary(gnuplot_binary), img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         gnuplot_failure(false), gnuplot_warning(true) {
    
      
    }
    
    void render(const vector<Block>& blocks) {
        Point centr(img.cols/2, img.rows/2);
        
        vector<Ordered_point> sagittal;
        vector<Ordered_point> meridional;
        for (size_t i=0; i < blocks.size(); i++) {
        
            
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                Point ec = blocks[i].get_edge_centroid(k);
                
                Point udir = ec - centr;
                double radial_len = norm(udir);
                Point dir = udir * (1.0/radial_len);

                Point norm = blocks[i].get_normal(k);
                double delta = dir.x*norm.x + dir.y*norm.y;
                
                const vector<double>& sfr = blocks[i].get_sfr(k);
                
                double contrast = sfr[12]; // arbitrary, TODO: select some (several?) resolutions
                
                double angle_to_radial = acos(fabs(delta))/M_PI*180.0;
                
                if (fabs(delta) < cos(15.0/180.0*M_PI)) { // edge perp to tangent
                    sagittal.push_back(Ordered_point(radial_len, contrast));
                } 
                if (fabs(fabs(delta) - 1) < cos(15.0/180.0*M_PI)) { // edge perp to radial : TODO: check math
                    meridional.push_back(Ordered_point(radial_len, contrast));
                }
                
            }
        }    
        sort(sagittal.begin(), sagittal.end());
        sort(meridional.begin(), meridional.end());
        
        FILE* fraw = fopen((wdir + "lp_raw.txt").c_str(), "wt");
        for (size_t i=0; i < sagittal.size(); i++) {
            fprintf(fraw, "%lf %lf\n", sagittal[i].first, sagittal[i].second);
        }
        fprintf(fraw, "\n\n");
        for (size_t i=0; i < meridional.size(); i++) {
            fprintf(fraw, "%lf %lf\n", meridional[i].first, meridional[i].second);
        }
        fclose(fraw);
        
        FILE* fout = fopen((wdir + prname).c_str(), "wt");
        
        printf("got %d sagittal / %d meridional samples\n", (int)sagittal.size(), (int)meridional.size());
        
        vector<Ordered_point> s_fitted;
        vector<Ordered_point> s_spread;
        
        lsfit(sagittal, s_fitted, s_spread);
        
        vector<Ordered_point> m_fitted;
        vector<Ordered_point> m_spread;
        
        lsfit(meridional, m_fitted, m_spread);
        
        fprintf(fout, "#sagittal curve\n");
        for (size_t i=0; i < s_fitted.size(); i++) {
            fprintf(fout, "%lf %lf\n", s_fitted[i].first, s_fitted[i].second);
        }
        
        fprintf(fout, "\n\n#meridional curve\n");
        for (size_t i=0; i < m_fitted.size(); i++) {
            fprintf(fout, "%lf %lf\n", m_fitted[i].first, m_fitted[i].second);
        }
        
        
        fprintf(fout, "\n\n#sagittal bounds\n");
        for (size_t i=0; i < s_fitted.size(); i++) {
            fprintf(fout, "%lf %lf\n", s_fitted[i].first, s_spread[i].first);
        }
        for (int i=(int)s_fitted.size()-1; i >= 0; i--) {
            fprintf(fout, "%lf %lf\n", s_fitted[i].first, s_spread[i].second);
        }
        
        
        fprintf(fout, "\n\n#meridional bounds\n");
        for (size_t i=0; i < m_fitted.size(); i++) {
            fprintf(fout, "%lf %lf\n", m_fitted[i].first, m_spread[i].first);
        }
        for (int i=(int)m_fitted.size()-1; i >= 0; i--) {
            fprintf(fout, "%lf %lf\n", m_fitted[i].first, m_spread[i].second);
        }
        fclose(fout);
    }
    
    void lsfit(const vector<Ordered_point>& in_data, vector<Ordered_point>& recon, 
        vector<Ordered_point>& spread, int recon_samples=64) {
        
        recon  = vector<Ordered_point> (recon_samples);
        spread  = vector<Ordered_point> (recon_samples, Ordered_point(0,-1.0));
        
        vector<Ordered_point> data(in_data);
        
        double x_span = data.back().first - data.front().first;
        
        for (int i=0; i < 5; i++) {
            data.push_back(Ordered_point(in_data[i].first - 0.01*x_span, in_data[i].second));
            data.push_back(Ordered_point(in_data[i].first - 0.015*x_span, in_data[i+1].second));
            data.push_back(Ordered_point(in_data[in_data.size()-1-i].first + 0.01*x_span, in_data[in_data.size()-1-i].second));
            data.push_back(Ordered_point(in_data[in_data.size()-1-i].first + 0.015*x_span, in_data[in_data.size()-1-i-1].second));
        }
        sort(data.begin(), data.end());
        
        
        int tdim = 2;
        double scale = 5.0;
        double h = x_span/scale; 
        double bin_width = x_span/double(recon.size()-1);
        
        size_t upper_idx = 0;
        size_t lower_idx = 0;
        double prev_scale = scale;
        for (size_t q=0; q < recon.size(); q++) {
            prev_scale = scale;
            if (q < recon.size()*0.2 || q >= recon.size()*0.8) {
                tdim = 3;
                scale = 5.0;
            } else {
                tdim = 4;
                scale = 5.0;
            }
            h = x_span/scale;
            if (scale != prev_scale) {
                lower_idx = 0;
                upper_idx = 0;
            }
        
            double xmid = (double(q)-0.5)*x_span/double(recon.size()-2) + in_data.front().first;
            
            for (; lower_idx < data.size() && data[lower_idx].first < (xmid-h); lower_idx++);
            for (; upper_idx < data.size() && data[upper_idx].first < (xmid+h); upper_idx++);
            
            vector< vector<double> > cov(tdim, vector<double>(tdim, 0.0));
            vector<double> b(tdim, 0.0);
            vector<double> a(tdim);
            
            for (size_t r=lower_idx; r < upper_idx; r++) { 
                double sh = 1.0/scale;
                double x = (data[r].first - xmid)/x_span; 
                if (fabs(x) > 1) x /= fabs(x);
                double w = fabs(x)*fabs(x)*fabs(x)/(sh*sh*sh);
                w = fabs(x) < sh ? (1 - w)*(1 - w)*(1 - w) : 0;
                
                a[0] = w;
                double px = x;
                for (int jj=1; jj < tdim; jj++) {
                    a[jj] = w*px;
                    px *= x;
                }
                for (int col=0; col < tdim; col++) {
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov[col][icol] += a[col]*a[icol];
                    }
                    b[col] += a[col]*data[r].second*w; // build rhs of system : A'*b
                }
            }
            
            // hardcode cholesky decomposition
            bool singular = false;
            vector<double> ldiag(tdim, 0.0);
            for (int i=0; i < tdim && !singular; i++) {
                for (int j=i; j < tdim && !singular; j++) {
                    double sum = cov[i][j];
                    for (int k=i-1; k >= 0; k--) {
                        sum -= cov[i][k]*cov[j][k];
                    }
                    if (i == j) {
                        if (sum <= 0.0) {
                            singular = true;
                        } else {
                            ldiag[i] = sqrt(sum);
                        }
                    } else {
                        cov[j][i] = sum/ldiag[i];
                    }
                }
            }
            
            if (!singular) {
                // hardcode backsubstitution
                vector<double> x(tdim, 0.0);
                for (int i=0; i < tdim; i++) {
                    double sum = b[i];
                    for (int k=i-1; k >= 0; k--) {
                        sum -= cov[i][k]*x[k];
                    }
                    x[i] = sum/ldiag[i];
                }
                for (int i=tdim-1; i >= 0; i--) {
                    double sum = x[i];
                    for (int k=i+1; k < tdim; k++) {
                        sum -= cov[k][i]*x[k];
                    }
                    x[i] = sum/ldiag[i];
                }
            
                
                recon[q].first = xmid;
                recon[q].second = x[0]; // x is centred on xmid, so parabola constant is midpoint extimate
                
                vector<double> residuals;
                for (size_t r=lower_idx; r < upper_idx; r++) {
                    // only compute residuals on really close points?
                    double cent_x = data[r].first - xmid;
                    
                    if (fabs(cent_x) < 2*bin_width) {
                        double delta = data[r].second - polyeval(cent_x/x_span, x);
                        residuals.push_back(fabs(delta));
                    }
                }
                if (residuals.size() < 10) {
                    for (size_t r=lower_idx; r < upper_idx; r++) {
                        double cent_x = data[r].first - xmid;
                        
                        if (fabs(cent_x) < 4*bin_width) {
                            double delta = data[r].second - polyeval(cent_x/x_span, x);
                            residuals.push_back(fabs(delta));
                        }
                    }
                }
                if (residuals.size() > 2) {
                    sort(residuals.begin(), residuals.end());
                    double p90 = residuals[lrint((residuals.size() - 1)*0.9)];
                    double p25 = residuals[lrint((residuals.size() - 1)*0.25)];
                    spread[q].first = p90 - p25;
                    spread[q].second = p90 - p25;
                }
            }
        }
        
        size_t nz_idx = 0;
        while (nz_idx < spread.size() && spread[nz_idx].second == -1) nz_idx++;
        if (nz_idx == spread.size()) { // no nonzero points at all!
            nz_idx = 0;
            spread.front() = spread.back() = Ordered_point(1,1);
        } else {
            for (size_t q=0; q < nz_idx; q++) { // pad out leading zeros with first nonzero
                spread[q] = spread[nz_idx];
            }
        }
        size_t lnz_idx = spread.size() - 1;
        while (lnz_idx > 0 && spread[lnz_idx].second == -1) lnz_idx--;
        for (size_t q=lnz_idx; q < spread.size(); q++) { // pad out trailing zeros with last nonzero
            spread[q] = spread[lnz_idx];
        }
        // now we have valid endpoints, so we can interpolate
        size_t prev_nz = nz_idx;
        for (size_t q=nz_idx+1; q < lnz_idx; q++) {
            if (spread[q].second == -1) {
                // find next nonzero index
                size_t next_nz = q+1;
                while (next_nz < lnz_idx && spread[next_nz].second == -1) next_nz++;
                double slope = (spread[next_nz].second - spread[prev_nz].second) / (recon[next_nz].first - recon[prev_nz].first);
                double offset = spread[prev_nz].second - slope*recon[prev_nz].first;
                while (q < next_nz) {
                    spread[q].first = spread[q].second = recon[q].first * slope + offset;
                    q++;
                }
            } 
            prev_nz = q;
        }
        
        // now hit recon with SG smoothing
        const int sgh = 7;
        vector<double> recon_smoothed(recon.size(), 0);
        const double* sgw = 0;
        for (int q=0; q < (int)recon.size(); q++) {
            if (q < sgh || q > (int(spread.size()) - 1 - sgh)) {
                recon_smoothed[q] = recon[q].second;
            } else {
                //const int stride = 3;
                int filter_order = 5; //std::min(5, (q-5)/stride);
                sgw = savitsky_golay[filter_order];
                for (int x=-sgh; x <= sgh; x++) { 
                    // since magnitude has extra samples at the end, we can safely go past the end
                    recon_smoothed[q] += recon[q+x].second * sgw[x+sgh];
                }
            }
        }
        for (int q=0; q < (int)recon.size(); q++) {
            recon[q].second = recon_smoothed[q];
        }
        
        // apply some exponential smoothing to spread
        vector<double> smoothed(spread.size(), 0);
        double fsmooth = spread.front().first;
        double bsmooth = spread.back().first;
        double alpha = 0.35;
        vector<double> fs(spread.size(), 0);
        vector<double> bs(spread.size(), 0);
        for (size_t q=0; q < spread.size(); q++) {
            fsmooth = (1 - alpha)*fsmooth + alpha * spread[q].first;
            bsmooth = (1 - alpha)*bsmooth + alpha * spread[spread.size()-1 - q].first;
            fs[q] = fsmooth;
            bs[spread.size()-1 - q] = bsmooth;
        }
        for (size_t q=0; q < spread.size(); q++) {
            smoothed[q] = 0.5*fs[q] + 0.5*bs[q];
            if (q < 10) {
                smoothed[q] = bs[q];
            }
            if (q > (spread.size()-1-10)) {
                smoothed[q] = fs[q];
            }
        }
        
        // every element of spread is now interpolated, so add spread to recon
        for (size_t q=0; q < spread.size(); q++) {
            spread[q].first = recon[q].second + smoothed[q];
            spread[q].second = recon[q].second - smoothed[q];
        }
    }
    
  private:
    double angle_reduce(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    inline double polyeval(double x, const vector<double>& a) const {
        double px = x;
        double s = a[0];
        for (size_t i=1; i < a.size(); i++) {
            s += a[i]*px;
            px *= x;
        }
        return s;
    }
    
    string wdir;
    string prname; 
    string pfname;
    string gnuplot_binary;
    const cv::Mat& img;            
    bool    lpmm_mode;
    double  pixel_size;
    bool gnuplot_failure;
    bool gnuplot_warning;
};

#endif
