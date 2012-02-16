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
#include <math.h>
#include <stdlib.h>
#include "render.h"
#include "noise_source.h"

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

const int border = 30;


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



// functor for tbb
class Render_rows {
  public:
    Render_rows(cv::Mat& in_img, const Render_rectangle& in_r, Noise_source& noise_source, 
        double contrast_reduction=0.05, bool gamma_correct=true, bool use_16bit=false,
        int buffer_border=30)
     : img(in_img), rect(in_r), noise_source(noise_source),
       gamma_correct(gamma_correct),contrast_reduction(contrast_reduction),
       use_16bit(use_16bit), buffer_border(buffer_border) {
     
        
    }


    inline void putpixel(int row, int col, double value) const {

        value = noise_source.sample(value, row*img.cols + col);

        if (value < 0.0) {
            value = 0.0;
        }
        if (value > 1.0) {
            value = 1.0;
        }

        if (use_16bit) {
            img.at<uint16_t>(row, col) = lrint(value*65535);
        } else {
            if (gamma_correct) {
                img.at<uchar>(row, col) = reverse_gamma(value);
            } else {
                img.at<uchar>(row, col) = lrint(value*255);
            }
        }
    }
     
    void operator()(const blocked_range<size_t>& r) const {
        // assume a dynamic range of 5% to 95% of full scale
        double object_value = contrast_reduction / 2.0;
        double background_value = 1 - object_value;

        for (size_t row=r.begin(); row != r.end(); ++row) {

            if (int(row) < buffer_border || int(row) > img.rows - buffer_border) {
                for (int col = 0; col < img.cols; col++) {
                    putpixel(row, col, background_value);
                }
            } else {
                for (int col = 0; col < buffer_border; col++) {
                    putpixel(row, col, background_value);
                }
                for (int col = img.cols - buffer_border; col < img.cols; col++) {
                    putpixel(row, col, background_value);
    
                }
                for (int col = buffer_border; col < img.cols-buffer_border; col++) {
                    
                    double rval = 0;
                    rval = rect.evaluate(col, row, object_value, background_value);
                    
                    putpixel(row, col, rval);
                }
            }
        }
    } 
     
    cv::Mat& img;
    const Render_rectangle& rect;
    Noise_source& noise_source;
    bool gamma_correct;
    double contrast_reduction;
    bool use_16bit;
    int buffer_border;
};


//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    int width  = 300;
    int height = 300;
    
    double sigma = 0.3;
    
    double theta = 4.0/180.0*M_PI;
    
    int rseed = time(0);
    
    stringstream ss;
    ss << genrectangle_VERSION_MAJOR << "." << genrectangle_VERSION_MINOR;
    
    TCLAP::CmdLine cmd("Generate rectangles with known MTF50 values", ' ', ss.str());
    TCLAP::ValueArg<std::string> tc_out_name("o", "output", "Output file name", false, "rect.png", "filename", cmd);
    TCLAP::ValueArg<double> tc_theta("a", "angle", "Orientation angle (degrees)", false, 4.0, "angle (degrees)", cmd);
    TCLAP::ValueArg<int> tc_seed("s", "seed", "Noise random seed", false, time(0), "seed", cmd);
    TCLAP::ValueArg<double> tc_noise("n", "noise", "Noise magnitude (linear standard deviation, range [0,1])", false, 0.01, "std. dev", cmd);
    TCLAP::ValueArg<double> tc_blur("b", "blur", "Blur magnitude (linear standard deviation, range [0.185, +inf))", false, 0.374781, "std. dev", cmd);
    TCLAP::ValueArg<double> tc_mtf("m", "mtf50", "Desired MTF50 value (range (0, 1.0])", false, 0.3, "mtf50", cmd);
    TCLAP::ValueArg<double> tc_cr("c", "contrast", "Contrast reduction [0,1]", false, 0.1, "factor", cmd);
    TCLAP::ValueArg<double> tc_dim("d", "dimension", "Dimension of the image, in pixels", false, 100, "dimension", cmd);
    TCLAP::ValueArg<double> tc_yoff("y", "yoffset", "Subpixel y offset [0,1]", false, 0, "pixels", cmd);
    TCLAP::ValueArg<double> tc_xoff("x", "xoffset", "Subpixel x offset [0,1]", false, 0, "pixels", cmd);
    TCLAP::ValueArg<double> tc_read_noise("", "read-noise", "Read noise magnitude (linear standard deviation, in electrons, range [0,+inf))", false, 3.0, "std. dev", cmd);
    TCLAP::ValueArg<double> tc_pattern_noise("", "pattern-noise", "Fixed pattern noise magnitude (linear percentage of signal, range [0,1])", false, 0.02, "percentage", cmd);
    TCLAP::ValueArg<double> tc_adc_gain("", "adc-gain", "ADC gain (linear electrons-per-DN, range (0,+inf))", false, 1.5, "electrons", cmd);
    TCLAP::ValueArg<int> tc_adc_depth("", "adc-depth", "ADC depth (number of bits, range [1,32])", false, 14, "bits", cmd);

    TCLAP::SwitchArg tc_linear("l","linear","Generate output image with linear intensities (default is sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_16("","b16","Generate linear 16-bit image", cmd, false);
    
    cmd.parse(argc, argv);
    
    rseed = tc_seed.getValue();
    srand(rseed);
    
    theta = tc_theta.getValue() / 180.0 * M_PI;
    
    double rwidth = 0;
    double rheight = 0;

    rwidth = tc_dim.getValue();
    rheight = rwidth*3/4;

    double diag = sqrt(rwidth*rwidth + rheight*rheight);
    width = lrint(diag + 32) + 2*border;
    height = width += width % 2;
    
    if (tc_mtf.isSet() && tc_blur.isSet()) {
        printf("Warning: you can not specify both blur and mtf50 values; choosing mtf50 value, and proceeding ...\n");
    }
    
    double mtf = tc_mtf.getValue();
    if (tc_mtf.isSet()) {
        sigma = sqrt( log(0.5)/(-2*M_PI*M_PI*mtf*mtf) );
    } else {
        sigma = tc_blur.getValue();
        mtf   = sqrt( log(0.5)/(-2*M_PI*M_PI*sigma*sigma) );
    }
    
    if (sigma < 0.185) {
        printf("It does not make sense to set blur below 0.185; you are on your own ...\n");
    }
    
    bool use_gamma = !tc_linear.getValue();
    bool use_16bit = tc_16.getValue();
    
    if (use_gamma && use_16bit) {
        printf("Setting both gamma and 16-bit output not supported. Choosing to disable gamma, and enable 16-bit output\n");
        use_gamma = false;
    }

    bool use_sensor_model = false;
    if (tc_read_noise.isSet() || tc_pattern_noise.isSet() ||
        tc_adc_depth.isSet()  || tc_adc_depth.isSet() ) {
        use_sensor_model = true;
    }

    if (tc_noise.isSet() && use_sensor_model) {
        printf("Noise sigma set, but full sensor noise model parameters also specified.\n"
               "Ignoring noise sigma, and going with full sensor model instead.\n");
    }

    printf("output filename = %s, sigma = %lf (or mtf50 = %lf), theta = %lf degrees, seed = %d, ", 
        tc_out_name.getValue().c_str(), sigma, mtf, theta/M_PI*180, rseed
    );
    if (use_sensor_model) {
        printf("\n\t full sensor noise model, with read noise = %.1lf electrons, fixed pattern noise fraction = %.3lf,\n"
               "\t adc gain = %.2lf e/DN, adc depth = %d bits\n",
               tc_read_noise.getValue(), tc_pattern_noise.getValue(), 
               tc_adc_gain.getValue(), tc_adc_depth.getValue()
        );
    } else {
        printf("\n\t additive Gaussian noise with sigma = %lf\n", tc_noise.getValue());
    }
    
    printf("\t output in sRGB gamma = %d, intensity range [%lf, %lf], 16-bit output:%d, dimension: %dx%d\n",
        use_gamma, tc_cr.getValue()/2.0, 1 - tc_cr.getValue()/2.0, use_16bit, width, height
    );


    Render_rectangle rect(
        width*0.5 + tc_xoff.getValue(), 
        height*0.5 + tc_yoff.getValue(),
        rwidth,
        rheight,
        M_PI/2 - theta,
        sigma
    );
    
    cv::Mat img;
    if (use_16bit) {
        img = cv::Mat(height, width, CV_16UC1);
    } else {
        img = cv::Mat(height, width, CV_8UC1);
    }

    
    Noise_source* ns = 0;

    if (use_sensor_model) {
        ns = new Sensor_model_noise(
            img.rows*img.cols,
            tc_read_noise.getValue(),
            tc_pattern_noise.getValue(),
            tc_adc_gain.getValue(),
            tc_adc_depth.getValue()
        );
    } else {
        ns = new Additive_gaussian_noise(img.rows*img.cols, tc_noise.getValue());
    }

    Render_rows rr(img, rect, *ns, tc_cr.getValue(), use_gamma, use_16bit, border);
    parallel_for(blocked_range<size_t>(size_t(0), height), rr); 
    
    double a = 1.0/(sigma*sigma);
    printf("MTF curve:  exp(%.8lg*x*x)\n", -2*M_PI*M_PI/a);
    printf("PSF : exp(-x*x/%lg)\n", 2*sigma*sigma);
    printf("MTF50 = %lf\n", mtf);    

    imwrite(tc_out_name.getValue(), img);

    delete ns;
    
    return 0;
}
