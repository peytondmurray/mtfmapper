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
#ifndef MTF_RENDERER_GRID_H
#define MTF_RENDERER_GRID_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "gaussfilter.h"
#include <cvwimage.h>

class Mtf_renderer_grid : public Mtf_renderer {
  public:
    Mtf_renderer_grid(const std::string& wdir, const std::string& in_fname, 
        const std::string& gnuplot_binary, const cv::Mat& img,
        bool lpmm_mode, double pixel_size)
      :  wdir(wdir), fname(in_fname), gnuplot_binary(gnuplot_binary), 
         img_y(img.rows), img_x(img.cols), img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         gnuplot_failure(false), gnuplot_warning(true),
         m_lower(0), m_upper(0) {

        const int base_grid_size = 200;
        if (img.rows > img.cols) {
            grid_y = base_grid_size;
            grid_x = base_grid_size * img.cols / img.rows;
        } else {
            grid_x = base_grid_size;
            grid_y = base_grid_size * img.rows / img.cols;
        }
    }
    
    void set_gnuplot_warning(bool gnuplot) {
        gnuplot_warning = gnuplot;
    }
    
    void render(const vector<Block>& blocks) {

        // first check if we have enough good data to generate a plot
        vector<double> allvals;
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0) {
                    allvals.push_back(val);
                }
            }
        }

        if (allvals.size() < 20) {
            printf("Too few valid edges found. No surface can be generated\n");
            return;
        }

        sort(allvals.begin(), allvals.end());
        m_upper = allvals[97*allvals.size()/100];
        m_lower = allvals[5*allvals.size()/100];

        cv::Mat grid_mer(grid_y, grid_x, CV_32FC1, 0.0);
        extract_mtf_grid(MERIDIONAL, grid_mer, blocks);
        cv::Mat grid_sag(grid_y, grid_x, CV_32FC1, 0.0);
        extract_mtf_grid(SAGITTAL, grid_sag, blocks);

        FILE* file = fopen((wdir+fname).c_str(), "wt");
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                fprintf(file, "%lf %lf %.3lf\n", 
                    x*img.cols/grid_x/pixel_size, 
                    y*img.rows/grid_y/pixel_size, 
                    grid_mer.at<float>(y,x)*pixel_size
                );
            }
            fprintf(file, "\n");
        }
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                fprintf(file, "%lf %lf %.3lf\n", 
                    x*img.cols/grid_x/pixel_size, 
                    (y*img.rows/grid_y + img.rows)/pixel_size, 
                    grid_sag.at<float>(y,x)*pixel_size
                );
            }
            fprintf(file, "\n");
        }
        fclose(file);

        const int width_in_pixels = 600;
        const int width_in_pixels_3d = 900;
        const int height_in_pixels_3d = 1200;
        FILE* gpf = fopen((wdir + std::string("grid.gnuplot")).c_str(), "wt");
        fprintf(gpf, "set yrange [] reverse\n");
        fprintf(gpf, "set palette defined (0 1 1 1, 1 0 0 1, 3 1 1 0, 4 1 0 0, 6 0 1 0)\n");
        fprintf(gpf, "set pm3d at bs depthorder interpolate 2,2\n");
        fprintf(gpf, "set cbrange [%lf:%lf]\n", m_lower*pixel_size, m_upper*pixel_size);
        fprintf(gpf, "set xlab \"column (%s)\"\n", lpmm_mode ? "mm" : "pixels");
        fprintf(gpf, "set ylab \"row (%s)\"\n",  lpmm_mode ? "mm" : "pixels");
        fprintf(gpf, "set term png size %d, %d\n", width_in_pixels, (int)lrint(width_in_pixels*2*grid_mer.rows/double(grid_mer.cols)));
        fprintf(gpf, "set output \"%sgrid_image.png\"\n", wdir.c_str());
        fprintf(gpf, "set multiplot\n");
        fprintf(gpf, "set size 1,0.5\n");   
        fprintf(gpf, "set origin 0.0,0.5\n");
        fprintf(gpf, "set title \"Meridional\"\n");
        fprintf(gpf, "plot [0:%lf][0:%lf] \"%s\" t \"MTF50 (%s)\" w image\n", 
                (img.cols)/pixel_size, (img.rows-1)/pixel_size, (wdir+fname).c_str(),
                lpmm_mode ? "lp/mm" : "c/p"
        );
        fprintf(gpf, "set origin 0.0,0.0\n");
        fprintf(gpf, "set title \"Sagittal\"\n");
        fprintf(gpf, "plot [0:%lf][0:%lf] \"%s\" u 1:($2-%lf):3 t \"MTF50 (%s)\" w image\n", 
                img.cols/pixel_size, img.rows/pixel_size, 
                (wdir+fname).c_str(),
                img.rows/pixel_size,
                lpmm_mode ? "lp/mm" : "c/p"
        );
        fprintf(gpf, "unset multiplot\n");
        fprintf(gpf, "set term png size %d, %d font \"arial,9\"\n", width_in_pixels_3d, height_in_pixels_3d);
        fprintf(gpf, "set output \"%sgrid_surface.png\"\n", wdir.c_str());
        fprintf(gpf, "unset xlab\n");
        fprintf(gpf, "unset ylab\n");
        fprintf(gpf, "set multiplot\n");
        fprintf(gpf, "set ticslevel %lf\n", m_lower);
        fprintf(gpf, "set view 25, 350\n");
        fprintf(gpf, "set title \"Meridional\"\n");
        fprintf(gpf, "set size 1,0.5\n");   
        fprintf(gpf, "set origin 0.0,0.5\n");
        fprintf(gpf, "splot [0:%lf][0:%lf] \"%s\" w d notitle\n", 
                img.cols/pixel_size, 
                (img.rows-1)/pixel_size, 
                (wdir+fname).c_str()
        );
        fprintf(gpf, "set view 25, 350\n");
        fprintf(gpf, "set title \"Sagittal\"\n");
        fprintf(gpf, "set origin 0.0,0.0\n");
        fprintf(gpf, "splot [0:%lf][0:%lf] \"%s\" u 1:($2-%lf):3 w d notitle\n", 
                img.cols/pixel_size, 
                (img.rows-1)/pixel_size, 
                (wdir+fname).c_str(), 
                img.rows/pixel_size
        );
        fprintf(gpf, "unset multiplot\n");
        fclose(gpf);
        
        char* buffer = new char[1024];
        #ifdef _WIN32
        sprintf(buffer, "\"\"%s\" \"%sgrid.gnuplot\"\"", gnuplot_binary.c_str(), wdir.c_str());
        #else
        sprintf(buffer, "\"%s\" \"%sgrid.gnuplot\"", gnuplot_binary.c_str(), wdir.c_str());
        #endif
        int rval = system(buffer);
        if (rval != 0) {
            printf("Failed to execute gnuplot (error code %d)\n", rval);
            printf("You can try to execute [%s] to render the plots manually\n", buffer);
            gnuplot_failure = true;
        } else {
            printf("Gnuplot plot completed successfully. Look for grid_image.png and grid_surface.png\n");
        }
        
        delete [] buffer;
        
    }
    
    bool gnuplot_failed(void) {
        return gnuplot_failure;
    }

  private:

    typedef enum {MERIDIONAL, SAGITTAL, NEITHER} Edge_type;

    void extract_mtf_grid(Edge_type target_edge_type, cv::Mat& grid, const vector<Block>& blocks) {
        Point centr(0,0);
        for (size_t i=0; i < blocks.size(); i++) {
            centr += blocks[i].get_centroid();
        }
        centr = centr*(1.0/double(blocks.size()));
        
        vector<double> vals;
        cv::Mat grid_binary(grid_y, grid_x, CV_8UC1, 1);
        cv::Mat grid_indices(grid_y, grid_x, CV_32SC1, 1);
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0 && blocks[i].get_quality(k) >= 0.2) {
                    Point cent = blocks[i].get_edge_centroid(k);

                    Point dir = cent - centr;
                    dir = dir * (1.0/norm(dir));

                    Point norm = blocks[i].get_normal(k);
                    double delta = dir.x*norm.x + dir.y*norm.y;

                    Edge_type edge_type = NEITHER;
                    if (target_edge_type == MERIDIONAL) {
                        if (fabs(delta) > cos(45.0/180*M_PI)) {
                            edge_type = MERIDIONAL;
                        } else {
                            edge_type = SAGITTAL;
                        }
                    } else {
                        if (fabs(delta) > cos(65.0/180*M_PI)) {
                            edge_type = MERIDIONAL;
                        } else {
                            edge_type = SAGITTAL;
                        }
                    }
                    
                    int ly = lrint(cent.y*grid_y/img.rows); 
                    int lx = lrint(cent.x*grid_x/img.cols);

                    if (val > grid.at<float>(ly,lx) && edge_type == target_edge_type) {
                        grid.at<float>(ly,lx) = val;
                        vals.push_back(val);
                        grid_binary.at<uchar>(ly,lx) = 0;
                        grid_indices.at<int32_t>(ly,lx) = vals.size();
                    }
                }
                
            }
        }

        // todo: should add NMS to prevent connected components from forming
        // in grid_binary

        #if 0
        const int erosion_size = 3;
        cv::Mat element = cv::getStructuringElement( 
            cv::MORPH_ELLIPSE,
            cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
            cv::Point( erosion_size, erosion_size ) 
        );
        
        cv::Mat dest(grid);
        cv::dilate(grid, grid, element, cv::Point(-1,-1), 1);
        cv::blur(grid, dest, cv::Size(5,5));
        grid = dest;
        #else

        cv::Mat grid_labels;
        cv::Mat distmat;
        cv::distanceTransform(grid_binary, distmat, grid_labels, CV_DIST_L2, 5);

        cv::Mat flood_labels;
        grid_labels.convertTo(flood_labels, CV_32FC1);

        // now floodfill this image, using the grid indices as guide:
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                int val = grid_indices.at<int32_t>(y,x);
                if (val > 1) {
                    cv::floodFill(flood_labels, Point(x,y), cv::Scalar::all((float)val), 0, cv::Scalar(), cv::Scalar(), 4);
                }
            }
        }
        // replace indices with actual mtf50 values
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                grid.at<float>(y,x) = vals[int(flood_labels.at<float>(y,x))-1];
            }
        }

        sort(vals.begin(), vals.end());
        double thresh95 = vals[95*vals.size()/100];
        // suppress upper 5% of values (to filter outliers)
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                if (grid.at<float>(y,x) > thresh95) {
                    grid.at<float>(y,x) = thresh95;
                }
            }
        }

        cv::Mat grid2(grid_y, grid_x, CV_32FC1, 0.0);
        cv::Mat grid3(grid_y, grid_x, CV_32FC1, 0.0);
        cv::boxFilter(grid, grid3, grid.type(), cv::Size(7,7));

        const double dist_thresh = 6.5;
        double far_fraction = 0;
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                if (distmat.at<float>(y,x) > dist_thresh) {
                    far_fraction += 1.0;
                }
            }
        }
        far_fraction /= (grid_x*grid_y);
        const double far_fraction_threshold = 0.3;
        printf("far fraction = %lf\n", far_fraction);

        // todo: replace this loop with opencv masked copy?
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                if (distmat.at<float>(y,x) > dist_thresh) {
                    if (far_fraction > far_fraction_threshold) {
                        grid.at<float>(y, x) = 0;
                    } else {
                        grid.at<float>(y, x) = grid3.at<float>(y, x);
                    }
                } else {
                    //if (grid_binary.at<uchar>(y,x) != 0) {
                        grid.at<float>(y, x) = grid3.at<float>(y, x);
                    //}   
                }
            }
        }
        
        // perform another pass to smooth the edges of the
        // surface to remove discontinuities (which render badly)
        if (far_fraction > far_fraction_threshold) {
            // edge pixels are to be replaced with Gaussian smoothed values
            cv::Mat smoothed;
            cv::GaussianBlur(grid, smoothed, cv::Size(7,7), 1.5, 1.5);

            // find edge pixels through morphology
            cv::Mat thresh(grid.rows, grid.cols, CV_8UC1);
            for (size_t y=0; y < grid_y; y++) {
                for (size_t x=0; x < grid_x; x++) {
                    thresh.at<char>(y,x) = grid.at<float>(y,x) > 0 ? 255 : 0;
                }
            }
            const int erosion_size = 3;
            cv::Mat element = cv::getStructuringElement( 
                cv::MORPH_ELLIPSE,
                cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                cv::Point( erosion_size, erosion_size ) 
            );
            cv::Mat dilated;
            cv::Mat eroded;
            cv::dilate(thresh, dilated, element, cv::Point(-1,-1), 1);
            cv::erode(thresh, eroded, element, cv::Point(-1,-1), 1);
            thresh = dilated - eroded;

            smoothed.copyTo(grid, thresh);
        }
        #endif
    }

    static double angular_diff(double a, double b) {
        return acos(cos(a)*cos(b) + sin(a)*sin(b));
    }
    
    std::string wdir;
    std::string fname;
    std::string gnuplot_binary;
    size_t grid_x;
    size_t grid_y;
    
    double img_y;
    double img_x;
    const cv::Mat& img;

    bool    lpmm_mode;
    double  pixel_size;
    
    bool gnuplot_failure;
    bool gnuplot_warning;

    double m_lower;
    double m_upper;
};

#endif
