#ifndef MTF_RENDERER_GRID_H
#define MTF_RENDERER_GRID_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "gaussfilter.h"
#include <cvwimage.h>

class Mtf_renderer_grid : public Mtf_renderer {
  public:
    Mtf_renderer_grid(const std::string& in_fname, const cv::Mat& img)
      :  fname(in_fname), img_y(img.rows), img_x(img.cols), img(img) {
      
          if (img.rows > img.cols) {
              grid_y = 150;
              grid_x = 150 * img.cols / img.rows;
          } else {
              grid_x = 150;
              grid_y = 150 * img.rows / img.cols;
          }
      
    }
    
    void render(const vector<Block>& blocks) {
    
        cv::Mat grid(grid_y, grid_x, CV_32FC1, 0.0);
        
        
        vector<double> vals;
        for (size_t i=0; i < blocks.size(); i++) {
            for (size_t k=0; k < 4; k++) {
                double val = blocks[i].get_mtf50_value(k);
                if (val > 0) {
                    Point cent = blocks[i].get_edge_centroid(k);
                    
                    int ly = lrint(cent.y*grid_y/img.rows); 
                    int lx = lrint(cent.x*grid_x/img.cols);
                    
                    if (val > grid.at<float>(ly,lx)) {
                        grid.at<float>(ly,lx) = val;
                        vals.push_back(val);
                    }
                }
                
            }
        }
        
        if (vals.size() < 20) {
            printf("Too few valid edges found. No surface can be generated\n");
            return;
        }
        
        sort(vals.begin(), vals.end());
        double thresh = vals[95*vals.size()/100];
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                if (grid.at<float>(y,x) > thresh) {
                    grid.at<float>(y,x) = thresh;
                }
            }
        }
        
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

        FILE* file = fopen(fname.c_str(), "wt");
        for (size_t y=0; y < grid_y; y++) {
            for (size_t x=0; x < grid_x; x++) {
                fprintf(file, "%d %d %.3lf\n", int(x*img.cols/grid_x), int(y*img.rows/grid_y), grid.at<float>(y,x));
            }
            fprintf(file, "\n");
        }
        fclose(file);
        
        FILE* gpf = fopen("grid.gnuplot", "wt");
        fprintf(gpf, "set size ratio %lf\n", grid.rows / double(grid.cols));
        fprintf(gpf, "set palette rgbformulae 23,28,3 negative\n");
        fprintf(gpf, "set pm3d at bs depthorder interpolate 2,2\n");
        fprintf(gpf, "set xlab \"column (pixels)\"\n");
        fprintf(gpf, "set ylab \"row (pixels)\"\n");
        fprintf(gpf, "set term png\n");
        fprintf(gpf, "set output \"grid_image.png\"\n");
        fprintf(gpf, "plot [0:%d][0:%d] \"%s\" t \"MTF50 (c/p)\" w image\n", img.cols, img.rows, fname.c_str());
        fprintf(gpf, "set output \"grid_surface.png\"\n");
        fprintf(gpf, "set view 25, 350\n");
        fprintf(gpf, "splot [0:%d][0:%d] \"%s\" t \"MTF50 (c/p)\" w d\n", img.cols, img.rows, fname.c_str());
        fclose(gpf);
        
        printf("execute \"gnuplot grid.gnuplot\" to render the plots\n");
        
    }
    
    std::string fname;
    size_t grid_x;
    size_t grid_y;
    
    double img_y;
    double img_x;
    const cv::Mat& img;
};

#endif
