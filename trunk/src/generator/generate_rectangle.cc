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

inline double gamma(double x) { // x in [0,1]
    const double S_linear = 12.9232102;
    const double C_srgb = 0.04045;
    const double SRGB_a = 0.055;
    
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    
    if (x < C_srgb) {
        return x / S_linear;
    }
    
    return  pow( (x + SRGB_a)/(1+SRGB_a), 2.4 );
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
    Render_rows(cv::Mat& in_img, const Render_rectangle& in_r, double noise_sigma=0.05, double contrast_reduction=0.05, bool gamma_correct=true)
     : img(in_img), rect(in_r), noise(img.rows*img.cols, 0.0), 
       gamma_correct(gamma_correct),contrast_reduction(contrast_reduction) {
     
        double gc_noise_sigma = noise_sigma;
        //if (gamma_correct) {
        //    gc_noise_sigma = gamma(noise_sigma);
        //} 
        for (size_t i=0; i < size_t(img.rows*img.cols); i++) {
            noise[i] = rand_norm(0, gc_noise_sigma);
        }
    }
     
    void operator()(const blocked_range<size_t>& r) const {
        // assume a dynamic range of 5% to 95% of full scale
        double object_value = contrast_reduction / 2.0;
        double background_value = 1 - object_value;
        //if (gamma_correct) {
        //    background_value = gamma(background_value);
        //    object_value = gamma(object_value);
        //    printf("gamma correction on, new values are: %lf, %lf\n", background_value, object_value);
        //}
        for (size_t row=r.begin(); row != r.end(); ++row) {
        
            for (int col = 0; col < img.cols; col++) {
                int index = row * img.cols + col;
                
                double rval = 0;
                rval = rect.evaluate(col, row, object_value, background_value);
                
                rval += noise[index];
                
                if (rval < 0.0) {
                    rval = 0.0;
                }
                if (rval > 1.0) {
                    rval = 1.0;
                }
                
                if (gamma_correct) {
                    img.at<uchar>(row, col) = reverse_gamma(rval);
                } else {
                    img.at<uchar>(row, col) = lrint(rval*255);
                }
            }
        }
    } 
     
    cv::Mat& img;
    const Render_rectangle& rect;
    vector<double> noise;
    bool gamma_correct;
    double contrast_reduction;
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
    TCLAP::ValueArg<std::string> tc_out_name("o", "output", "Output file name", false, "rect.png", "string", cmd);
    TCLAP::ValueArg<double> tc_theta("a", "angle", "Orientation angle (degrees)", false, 4.0, "double", cmd);
    TCLAP::ValueArg<int> tc_seed("s", "seed", "Noise random seed", false, time(0), "int", cmd);
    TCLAP::ValueArg<double> tc_noise("n", "noise", "Noise magnitude (linear standard deviation, range [0,1])", false, 0.01, "double", cmd);
    TCLAP::ValueArg<double> tc_blur("b", "blur", "Blur magnitude (linear standard deviation, range [0.185, +inf))", false, 0.3, "double", cmd);
    TCLAP::ValueArg<double> tc_cr("c", "contrast", "Contrast reduction [0,1]", false, 0.1, "double", cmd);
    TCLAP::SwitchArg tc_gamma("g","gamma","Generate output image with sRGB gamma", cmd, true);
    
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
    printf("\t output in sRGB gamma = %d, intensity range [%lf, %lf]\n",
        tc_gamma.getValue(), tc_cr.getValue()/2.0, 1 - tc_cr.getValue()/2.0
    );
    
    const double rwidth = width / 4;
    const double rheight = height / 3;
    
    Render_rectangle rect(width*0.5 - 0.5*rwidth, height*0.5 - 0.5*rheight, rwidth, rheight, theta, sigma);
    
    cv::Mat img(height, width, CV_8UC1);
    
    Render_rows rr(img, rect, tc_noise.getValue(), tc_cr.getValue(), tc_gamma.getValue());
    parallel_for(blocked_range<size_t>(size_t(0), height), rr); 
    
    double a = 2.0/(2*sigma*sigma);
    printf("MTF curve:  exp(%le*x*x)\n", -2*M_PI*M_PI/a);
    double s = -2*M_PI*M_PI/a;
    printf("MTF50 = %lf\n", sqrt(log(0.5)/(s)));    

    imwrite(tc_out_name.getValue(), img);
    
    return 0;
}
