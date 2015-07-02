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
#include "include/mtf_renderer_stats.h"
#include "include/mtf_renderer_sfr.h"
#include "include/mtf_renderer_esf.h"
#include "include/mtf_tables.h"
#include "include/scanline.h"
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
            *it16 = ((uint16_t)(*it)) << 8;
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
    TCLAP::UnlabeledValueArg<std::string>  tc_in_name("<input_filename>", 
        "Input image file name (many extensions supported)", true, "input.png", "string", cmd
    );
    TCLAP::UnlabeledValueArg<std::string>  tc_wdir("<working_directory>", 
        "Working directory (for output files); \".\" is fine", true, ".", "string", cmd
    );
    TCLAP::SwitchArg tc_profile("p","profile","Generate MTF50 profile", cmd, false);
    TCLAP::SwitchArg tc_annotate("a","annotate","Annotate input image with MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_surface("s","surface","Generate MTF50 surface plots", cmd, false);
    TCLAP::SwitchArg tc_linear("l","linear","Input image is linear 8-bit (default for 8-bit is assumed to be sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_print("r","raw","Print raw MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_sfr("f","sfr","Store raw SFR curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_esf("e","esf","Store raw ESF and PSF curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_border("b","border","Add a border of 20 pixels to the image", cmd, false);
    TCLAP::SwitchArg tc_absolute("","absolute-sfr","Generate absolute SFR curve (MTF) i.s.o. relative SFR curve", cmd, false);
    TCLAP::SwitchArg tc_smooth("","nosmoothing","Disable SFR curve (MTF) smoothing", cmd, false);
    TCLAP::ValueArg<double> tc_angle("g", "angle", "Angular filter [0,360)", false, 0, "angle", cmd);
    TCLAP::ValueArg<double> tc_snap("", "snap-angle", "Snap-to angle modulus [0,90)", false, 1000, "angle", cmd);
    TCLAP::ValueArg<double> tc_thresh("t", "threshold", "Dark object threshold (0,1)", false, 0.75, "threshold", cmd);
    TCLAP::ValueArg<string> tc_gnuplot("", "gnuplot-executable", "Full path (including filename) to gnuplot executable ", false, "gnuplot", "filepath", cmd);
    TCLAP::ValueArg<double> tc_pixelsize("", "pixelsize", "Pixel size in microns. This also switches units to lp/mm", false, 1.0, "size", cmd);

    vector<string> allowed_bayer_subsets;
    allowed_bayer_subsets.push_back("red");
    allowed_bayer_subsets.push_back("green");
    allowed_bayer_subsets.push_back("blue");
    allowed_bayer_subsets.push_back("none");
    TCLAP::ValuesConstraint<string> bayer_constraints(allowed_bayer_subsets);
    TCLAP::ValueArg<std::string> tc_bayer("", "bayer", "Select Bayer subset", false, "none", &bayer_constraints );
    cmd.add(tc_bayer);
    
    cmd.parse(argc, argv);

    bool lpmm_mode = false;
    double pixel_size = 1;
    if (tc_pixelsize.isSet()) {
        printf("Info: Pixel size has been specified, measurement will be reported in lp/mm, rather than c/p\n");
        lpmm_mode = true;
        pixel_size = 1000 / tc_pixelsize.getValue();
        printf("working with=%lf pixels per mm\n", pixel_size);
    }

    if (!tc_profile.getValue() && !tc_annotate.getValue() && !tc_surface.getValue() && !tc_print.getValue() && !tc_sfr.getValue()) {
        printf("Warning: No output specified. You probably want to specify at least one of the following flags: [-r -p -a -s -f]\n");
    }

    cv::Mat cvimg;
	try {
	    cvimg = cv::imread(tc_in_name.getValue(),-1);
	} catch (const cv::Exception& ex) {
		cout << ex.what() << endl;
	}

	if (!cvimg.data) {
		printf("Fatal error: could not open input file <%s>.\nFile is missing, or not where you said it would be, or you do not have read permission.\n", tc_in_name.getValue().c_str());
		return -2;
	}
    
    if (cvimg.type() == CV_8UC3 || cvimg.type() == CV_16UC3) {
        printf("colour input image detected; converting to grayscale using 0.299R + 0.587G + 0.114B\n");
        cv::Mat dest;
        cv::cvtColor(cvimg, dest, CV_RGB2GRAY);  // force to grayscale
        cvimg = dest;
    }
    
    if (cvimg.type() == CV_8UC1) {
        printf("8-bit input image, upconverting %s\n", tc_linear.getValue() ? "with linear scaling" : "with sRGB gamma correction");
        convert_8bit_input(cvimg, !tc_linear.getValue());        
    } else {
        printf("16-bit input image, no upconversion required\n");
    }
   
    assert(cvimg.type() == CV_16UC1);


    if (tc_border.getValue()) {
        printf("The -b option has been specified, adding a 20-pixel border to the image\n");
        double max_val = 0;
        double min_val = 0;
        cv::minMaxLoc(cvimg, &min_val, &max_val);
        cv::Mat border;
        cv::copyMakeBorder(cvimg, border, 20, 20, 20, 20, cv::BORDER_CONSTANT, cv::Scalar((int)max_val));
        cvimg = border;
    }
    
    
    // process working directory
    std::string wdir(tc_wdir.getValue());
    if (wdir[wdir.length()-1]) {
        wdir = tc_wdir.getValue() + "/";
    }

	#ifdef _WIN32
	// on windows, mangle the '/' into a '\\'
	std::string wdm;
	for (size_t i=0; i < wdir.length(); i++) {
		if (wdir[i] == '/') {
			wdm.push_back('\\');
			wdm.push_back('\\');
		} else {
			wdm.push_back(wdir[i]);
		}
	}
	wdir = wdm;
	#endif
	
	// initialize the global apodization filter
	global_apodization_instance = new Apodization(512);
	// initialize the global mtf correction instance
	global_mtf_correction_instance = new Mtf_correction();

    cv::Mat masked_img;
    
    printf("Computing gradients ...\n");
    Gradient gradient(cvimg, false);
    
    printf("Thresholding image ...\n");
    int brad_S = 2*std::min(cvimg.cols, cvimg.rows)/3;
    double brad_threshold = tc_thresh.getValue();
    bradley_adaptive_threshold(cvimg, masked_img, brad_threshold, brad_S);
    
    printf("Component labelling ...\n");
    Component_labeller::zap_borders(masked_img);    
    Component_labeller cl(masked_img, 60, false, 8000);

    #if 0
    // blank out the largest block
    cv::Mat outimg = cvimg.clone();
    Boundarylist bl = cl.get_boundaries();
    for (Boundarylist::const_iterator it=bl.begin(); it != bl.end(); it++) {
        if (it->second.size() > 3000) {
            Point cent = centroid(it->second);
            printf("found a big one: %d (%lf,%lf)\n", (int)it->second.size(), cent.x, cent.y); 

            map<int, scanline> scanset;
            for (size_t i=0; i < it->second.size(); i++) {
                Point p = it->second[i];
                int iy = lrint(p.y);
                int ix = lrint(p.x);
                if (iy > 0 && iy < cvimg.rows && ix > 0 && ix < cvimg.cols) {
                    map<int, scanline>::iterator it2 = scanset.find(iy);
                    if (it2 == scanset.end()) {
                        scanline sl(ix,ix);
                        scanset.insert(make_pair(iy, sl));
                    }
                    if (ix < scanset[iy].start) {
                        scanset[iy].start = ix;
                    }
                    if (ix > scanset[iy].end) {
                        scanset[iy].end = ix;
                    }
                }
            }
            for (map<int,scanline>::const_iterator it3=scanset.begin(); it3 != scanset.end(); it3++) {
                printf("line %d, from %d to %d\n", it3->first, scanset[it3->first].start, scanset[it3->first].end);
                for (int col=scanset[it3->first].start+1; col < scanset[it3->first].end; col++) {
                    outimg.at<uint16_t>(it3->first, col) = 0;
                }
            }
        }
    }
    imwrite(string("out.tiff"), outimg);
    #endif
    
    if (cl.get_boundaries().size() == 0) {
        printf("No black objects found. Try a lower threshold value with the -t option.\n");
        return 0;
    }
    
    // now we can destroy the thresholded image
    masked_img = cv::Mat(1,1, CV_8UC1);
    
    Mtf_core mtf_core(cl, gradient, cvimg, tc_bayer.getValue());
    mtf_core.set_absolute_sfr(tc_absolute.getValue());
    mtf_core.set_sfr_smoothing(!tc_smooth.getValue());
    if (tc_snap.isSet()) {
        mtf_core.set_snap_angle(tc_snap.getValue()/180*M_PI);
    }
    
    Mtf_core_tbb_adaptor ca(&mtf_core);
    
    printf("Parallel MTF50 calculation\n");
    parallel_for(blocked_range<size_t>(size_t(0), mtf_core.num_objects()), ca); 
    
    
    // now render the computed MTF values
    if (tc_annotate.getValue()){
        Mtf_renderer_annotate annotate(cvimg, wdir + string("annotated.png"), lpmm_mode, pixel_size);
        annotate.render(mtf_core.get_blocks());
    }
    
    bool gnuplot_warning = true;
    bool few_edges_warned = false;
    
    if (tc_profile.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            printf("Warning: fewer than 10 edges found, so MTF50 surfaces/profiles will not be generated. Are you using suitable input images?\n");
            few_edges_warned = true;
        } else {
            Mtf_renderer_profile profile(
                wdir, 
                string("profile.txt"),
                string("profile_peak.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                lpmm_mode,
                pixel_size
            );
            profile.render(mtf_core.get_blocks());
            gnuplot_warning = !profile.gnuplot_failed();
        }
    }

    if (tc_surface.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            if (!few_edges_warned) {
                printf("Warning: fewer than 10 edges found, so MTF50 surfaces/profiles will not be generated. Are you using suitable input images?\n");
            }
        } else {
            Mtf_renderer_grid grid(
                wdir, 
                string("grid.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                lpmm_mode,
                pixel_size
            );
            grid.set_gnuplot_warning(gnuplot_warning);
            grid.render(mtf_core.get_blocks());
        }
    }
    
    if (tc_print.getValue()) {
        Mtf_renderer_print printer(
            wdir + string("raw_mtf_values.txt"), 
            tc_angle.isSet(), 
            tc_angle.getValue()/180.0*M_PI,
            lpmm_mode,
            pixel_size
        );
        printer.render(mtf_core.get_blocks());
    }

    if (tc_sfr.getValue() || tc_absolute.getValue()) {
        Mtf_renderer_sfr sfr_writer(
            wdir + string("raw_sfr_values.txt"), 
            lpmm_mode,
            pixel_size
        );
        sfr_writer.render(mtf_core.get_blocks());
    }
    
    if (tc_esf.getValue()) {
        Mtf_renderer_esf esf_writer(
            wdir + string("raw_esf_values.txt"), 
            wdir + string("raw_psf_values.txt")
        );
        esf_writer.render(mtf_core.get_blocks());
    }
    
    Mtf_renderer_stats stats(lpmm_mode, pixel_size);
    stats.render(mtf_core.get_blocks());
    
    return 0;
}
