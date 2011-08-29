#include <math.h>
#include <stdlib.h>
#include "render.h"

#include "cv.h"
#include "highgui.h"

#include "config.h"

#include "tclap/CmdLine.h"

#include <vector>
using std::vector;

#include "tbb/tbb.h"
using namespace tbb;

#include <iostream>
#include <string>

using std::string;
using std::stringstream;

using std::cout;
using std::endl;

#define SQR(x) ((x)*(x))


inline unsigned char reverse_gamma(double x) {
    const double C_linear = 0.0031308;
    const double S_linear = 12.9232102;
    const double SRGB_a = 0.055;
    
    if (x < C_linear) {
        return lrint(255 * x * S_linear);
    }
    return lrint(255 * ((1 + SRGB_a) * pow(x, 1.0/2.4) - SRGB_a));
}


double randdu(void) {
    return double(rand())/double(RAND_MAX);
}

double rand_norm(double m, double s) {
    double x = normal_sampler::Moro_norm_inv(randdu());
    return( m + x * s );
}
                                                                                                                                                                                                            

// functor for tbb
class Render_rows {
  public:
    Render_rows(cv::Mat& in_img, const Render_rectangle& in_r, double noise_sigma=0.005)
     : img(in_img), rect(in_r), noise(img.rows*img.cols, 0.0) {
     
        for (size_t i=0; i < size_t(img.rows*img.cols); i++) {
            noise[i] = rand_norm(0, noise_sigma);
        }
    }
     
    void operator()(const blocked_range<size_t>& r) const {
        for (size_t row=r.begin(); row != r.end(); ++row) {
        
            for (int col = 0; col < img.cols; col++) {
                int index = row * img.cols + col;
                
                double background = 0.01; // start with white background
                double rval = 0;
                rval = rect.evaluate(col, row, background);
                if (rval != transparent) {
                    background = rval;
                }
                
                background *= 0.9; // only use 90% of dynamic range
                background += 0.01; // shift black level to 5%
                background += noise[index];
                if (background < 0.0) {
                    background = 0.0;
                }
                if (background > 1.0) {
                    background = 1.0;
                }
                
                img.at<uchar>(row, col) = reverse_gamma(background);
            }
        }
    } 
     
    cv::Mat& img;
    const Render_rectangle& rect;
    vector<double> noise;
};


//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    const int width  = 300;
    const int height = 300;
    
    double sigma = 0.3;
    
    double theta = 4.0/180.0*M_PI;
    
    int rseed = time(0);
    
    stringstream ss;
    ss << genrectangle_VERSION_MAJOR << "." << genrectangle_VERSION_MINOR;
    
    TCLAP::CmdLine cmd("Generate rectangles with known MTF50 values", ' ', ss.str());
    TCLAP::ValueArg<std::string> tc_out_name("o", "output", "Output file name", false, "rect.png", "string");
    cmd.add(tc_out_name);
    TCLAP::ValueArg<double> tc_theta("a", "angle", "Orientation angle (degrees)", false, 4.0, "double");
    cmd.add(tc_theta);
    TCLAP::ValueArg<int> tc_seed("s", "seed", "Noise random seed", false, time(0), "int");
    cmd.add(tc_seed);
    TCLAP::ValueArg<double> tc_noise("n", "noise", "Noise magnitude (linear standard deviation, range [0,1])", false, 0.005, "double");
    cmd.add(tc_noise);
    TCLAP::ValueArg<double> tc_blur("b", "blur", "Blur magnitude (linear standard deviation, range [0.185, +inf))", false, 0.3, "double");
    cmd.add(tc_blur);
    
    cmd.parse(argc, argv);
    
    rseed = tc_seed.getValue();
    srand(rseed);
    
    theta = tc_theta.getValue() / 180.0 * M_PI;
    sigma = tc_blur.getValue();
    if (sigma < 0.185) {
        printf("It does not make sense to set blur below 0.185; you are on your own ...\n");
    }
    
    printf("output filename = %s, sigma = %lf, theta = %lf degrees, seed = %d, noise = %lf\n", 
        tc_out_name.getValue().c_str(), sigma, theta/M_PI*180, rseed, tc_noise.getValue()
    );
    
    const double rwidth = width / 4;
    const double rheight = height / 3;
    
    Render_rectangle rect(width*0.5 - 0.5*rwidth, height*0.5 - 0.5*rheight, rwidth, rheight, theta, sigma);
    
    cv::Mat img(height, width, CV_8UC1);
    
    Render_rows rr(img, rect, tc_noise.getValue());
    parallel_for(blocked_range<size_t>(size_t(0), height), rr); 
    
    double a = 2.0/(2*sigma*sigma);
    printf("MTF curve:  exp(%le*x*x)\n", -2*M_PI*M_PI/a);
    double s = -2*M_PI*M_PI/a;
    printf("MTF50 = %lf\n", sqrt(log(0.5)/(s)));    

    imwrite(tc_out_name.getValue(), img);
    
    return 0;
}
