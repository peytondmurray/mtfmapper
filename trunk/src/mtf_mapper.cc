#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <string.h>

using std::string;

#include "include/common_types.h"
#include "include/thresholding.h"
#include "include/gamma_lut.h"
#include "include/mtf_core.h"
#include "include/mtf_core_tbb_adaptor.h"
#include "include/mtf_renderer_annotate.h"
#include "include/mtf_renderer_profile.h"
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
    
    print_version_info();

    if (argc < 2) {
        printf("usage %s <filename>\n", argv[0]);
        return -1;
    }

    cv::Mat cvimg = cv::imread(string(argv[1]), 0);
    
    if (cvimg.type() == CV_8UC1) {
        convert_8bit_input(cvimg);        
    }
   
    assert(cvimg.type() == CV_16UC1);
    
    cv::Mat masked_img;
    
    printf("Thresholding image ...\n");
    int brad_S = 2*cvimg.cols/3;
    const double brad_threshold = 0.75;
    bradley_adaptive_threshold(cvimg, masked_img, brad_threshold, brad_S);
    
    printf("Computing gradients ...\n");
    Gradient gradient(cvimg, true);
    
    printf("Component labelling\n");
    Component_labeller::zap_borders(masked_img);    
    Component_labeller cl(masked_img, 100, true, 4000);
    
    // now we can destroy the thresholded image
    masked_img = cv::Mat(1,1, CV_8UC1);
    
    Mtf_core mtf_core(cl, gradient, cvimg);
    
    Mtf_core_tbb_adaptor ca(&mtf_core);
    
    printf("Parallel MTF50 calculation\n");
    parallel_for(blocked_range<size_t>(size_t(0), mtf_core.num_objects()), ca); 
    
    
    // now render the computed MTF values
    
    Mtf_renderer_annotate annotate(cvimg, string("annotated.png"));
    
    annotate.render(mtf_core.get_blocks());
    
    Mtf_renderer_profile profile(string("profile.txt"), string("profile_peak.txt"));
    
    profile.render(mtf_core.get_blocks());

    return 0;
}
