#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <string.h>

#include <tclap/CmdLine.h>

using std::string;
using std::stringstream;

#include "include/common_types.h"
#include "include/thresholding.h"
#include "include/gamma_lut.h"
#include "include/mtf_core.h"
#include "include/mtf_core_tbb_adaptor.h"
#include "include/mtf_renderer_annotate.h"
#include "include/mtf_renderer_profile.h"
#include "include/mtf_renderer_grid.h"
#include "include/mtf_renderer_print.h"
#include "config.h"

void convert_8bit_input(cv::Mat& cvimg, bool gamma_correct=true) {

    cv::Mat newmat = cv::Mat(cvimg.rows, cvimg.cols, CV_16UC1);
    cv::MatIterator_<uint16_t> it16 = newmat.begin<uint16_t>();
    cv::MatConstIterator_<uchar> it = cvimg.begin<uchar>();
    cv::MatConstIterator_<uchar> it_end = cvimg.end<uchar>();
    
    // can we determine what the current gamma is ?
    if (gamma_correct) {
        Gamma gamma;
        for(; it != it_end; ) {
            *it16 = gamma.linearize_gamma(*it);
            it++;
            it16++;
        }
    } else {
        for(; it != it_end; ) {
            *it16 = (uint16_t)(*it);
            it++;
            it16++;
        }
    }
    
    cvimg = newmat;
}

//-----------------------------------------------------------------------------
void print_version_info(void) {
    printf("MTF mapper version %d.%d\n", mtfmapper_VERSION_MAJOR, mtfmapper_VERSION_MINOR);
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv) {

    stringstream ss;
    ss << mtfmapper_VERSION_MAJOR << "." << mtfmapper_VERSION_MINOR;
    
    TCLAP::CmdLine cmd("Measure MTF50 values across edges of rectangular targets", ' ', ss.str());
    TCLAP::UnlabeledValueArg<std::string>  tc_in_name("input", 
        "Input image file name (many extensions supported)", true, "input.png", "string", cmd
    );
    TCLAP::SwitchArg tc_profile("p","profile","Generate MTF50 profile", cmd, true);
    TCLAP::SwitchArg tc_annotate("a","annotate","Annotate input image with MTF50 values", cmd, true);
    TCLAP::SwitchArg tc_surface("s","surface","Generate MTF50 surface plots", cmd, true);
    TCLAP::SwitchArg tc_linear("l","linear","Input image is linear 8-bit (default for 8-bit is assumed to be sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_print("r","raw","Print raw MTF50 values", cmd, false);
    TCLAP::ValueArg<double> tc_angle("g", "angle", "Angular filter [0,360)", false, 1000, "double", cmd);
    
    cmd.parse(argc, argv);

    cv::Mat cvimg = cv::imread(tc_in_name.getValue(),-1);
    
    if (cvimg.type() == CV_8UC3 || cvimg.type() == CV_16UC3) {
        printf("colour input image detected; converting to grayscale using 0.299R + 0.587G + 0.114B\n");
        cv::Mat dest;
        cv::cvtColor(cvimg, dest, CV_RGB2GRAY);  // force to grayscale
        cvimg = dest;
    }
    
    if (cvimg.type() == CV_8UC1) {
        printf("8-bit input image, upconverting\n");
        convert_8bit_input(cvimg, !tc_linear.getValue());        
    } else {
        printf("16-bit input image, no upconversion required\n");
    }
   
    assert(cvimg.type() == CV_16UC1);
    
    cv::Mat masked_img;
    
    printf("Thresholding image ...\n");
    int brad_S = 2*cvimg.cols/3;
    const double brad_threshold = 0.75;
    bradley_adaptive_threshold(cvimg, masked_img, brad_threshold, brad_S);
    
    printf("Computing gradients ...\n");
    Gradient gradient(cvimg, false);
    
    printf("Component labelling\n");
    Component_labeller::zap_borders(masked_img);    
    Component_labeller cl(masked_img, 100, false, 4000);
    
    // now we can destroy the thresholded image
    masked_img = cv::Mat(1,1, CV_8UC1);
    
    Mtf_core mtf_core(cl, gradient, cvimg);
    
    Mtf_core_tbb_adaptor ca(&mtf_core);
    
    printf("Parallel MTF50 calculation\n");
    parallel_for(blocked_range<size_t>(size_t(0), mtf_core.num_objects()), ca); 
    
    
    // now render the computed MTF values
    if (tc_annotate.getValue()){
        Mtf_renderer_annotate annotate(cvimg, string("annotated.png"));
        annotate.render(mtf_core.get_blocks());
    }
    
    if (tc_profile.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            printf("Warning: fewer than 10 edges found, so MTF50 profiles will not be generated. Are you using suitable input images?\n");
        } else {
            Mtf_renderer_profile profile(string("profile.txt"), string("profile_peak.txt"), cvimg);
            profile.render(mtf_core.get_blocks());
        }
    }

    if (tc_surface.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            printf("Warning: fewer than 10 edges found, so MTF50 surfaces will not be generated. Are you using suitable input images?\n");
        } else {
            Mtf_renderer_grid grid(string("grid.txt"),  cvimg);
            grid.render(mtf_core.get_blocks());
        }
    }
    
    if (tc_print.getValue()) {
        Mtf_renderer_print printer(string("raw_mtf_values.txt"), tc_angle.getValue() != 1000, tc_angle.getValue()/180.0*M_PI);
        printer.render(mtf_core.get_blocks());
    }
    
    return 0;
}
