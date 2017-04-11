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

#include "include/logger.h"
Logger logger;

#include "include/common_types.h"
#include "include/thresholding.h"
#include "include/gamma_lut.h"
#include "include/mtf_core.h"
#include "include/mtf_core_tbb_adaptor.h"
#include "include/mtf_renderer_annotate.h"
#include "include/mtf_renderer_profile.h"
#include "include/mtf_renderer_mfprofile.h"
#include "include/mtf_renderer_grid.h"
#include "include/mtf_renderer_print.h"
#include "include/mtf_renderer_stats.h"
#include "include/mtf_renderer_sfr.h"
#include "include/mtf_renderer_esf.h"
#include "include/mtf_renderer_edges.h"
#include "include/mtf_renderer_lensprofile.h"
#include "include/mtf_renderer_chart_orientation.h"
#include "include/mtf_renderer_focus.h"
#include "include/mtf_tables.h"
#include "include/scanline.h"
#include "include/distance_scale.h"
#include "include/autocrop.h"
#include "include/bayer.h"
#include "include/demosaic.h"
#include "include/stride_range.h"
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
    printf("MTF mapper version %d.%d.%d\n", mtfmapper_VERSION_MAJOR, mtfmapper_VERSION_MINOR, mtfmapper_VERSION_SUB);
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv) {

    stringstream ss;
    ss << mtfmapper_VERSION_MAJOR << "." << mtfmapper_VERSION_MINOR << "." << mtfmapper_VERSION_SUB;
    
    TCLAP::CmdLine cmd("Measure MTF50 values across edges of rectangular targets", ' ', ss.str());
    TCLAP::UnlabeledValueArg<std::string>  tc_in_name("<input_filename>", 
        "Input image file name (many extensions supported)", true, "input.png", "string", cmd
    );
    TCLAP::UnlabeledValueArg<std::string>  tc_wdir("<working_directory>", 
        "Working directory (for output files); \".\" is fine", true, ".", "string", cmd
    );
    TCLAP::SwitchArg tc_profile("p","profile","Generate MTF50 profile", cmd, false);
    TCLAP::SwitchArg tc_mf_profile("","mf-profile","Generate MTF50 profile (Manual focus chart)", cmd, false);
    TCLAP::SwitchArg tc_annotate("a","annotate","Annotate input image with MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_surface("s","surface","Generate MTF50 surface plots", cmd, false);
    TCLAP::SwitchArg tc_linear("l","linear","Input image is linear 8-bit (default for 8-bit is assumed to be sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_print("r","raw","Print raw MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_edges("q","edges","Print raw MTF50 values, grouped by edge location", cmd, false);
    TCLAP::SwitchArg tc_sfr("f","sfr","Store raw SFR curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_esf("e","esf","Store raw ESF and PSF curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_lensprof("","lensprofile","Render M/S lens profile plot", cmd, false);
    TCLAP::SwitchArg tc_chart_orientation("","chart-orientation","Visualize chart orientation relative to camera", cmd, false);
    TCLAP::SwitchArg tc_border("b","border","Add a border of 20 pixels to the image", cmd, false);
    TCLAP::SwitchArg tc_absolute("","absolute-sfr","Generate absolute SFR curve (MTF) i.s.o. relative SFR curve", cmd, false);
    TCLAP::SwitchArg tc_smooth("","nosmoothing","Disable SFR curve (MTF) smoothing", cmd, false);
    TCLAP::SwitchArg tc_autocrop("","autocrop","Automatically crop image to the chart area", cmd, false);
    TCLAP::SwitchArg tc_focus("","focus","Compute focus depth using special 'focus' type chart", cmd, false);
    TCLAP::SwitchArg tc_log_append("", "log-append", "Append to log file in stead of overwriting log file", cmd, false);
    TCLAP::SwitchArg tc_debug("", "debug", "Enable debug output messages", cmd, false);
    #ifdef MDEBUG
    TCLAP::SwitchArg tc_bradley("", "bradley", "Use Bradley thresholding i.s.o Sauvola thresholding", cmd, false);
    #endif
    TCLAP::ValueArg<double> tc_angle("g", "angle", "Angular filter [0,360)", false, 0, "angle", cmd);
    TCLAP::ValueArg<double> tc_snap("", "snap-angle", "Snap-to angle modulus [0,90)", false, 1000, "angle", cmd);
    TCLAP::ValueArg<double> tc_thresh("t", "threshold", "Dark object threshold (0,1), default 0.55", false, 0.55, "threshold", cmd);
    TCLAP::ValueArg<string> tc_gnuplot("", "gnuplot-executable", "Full path (including filename) to gnuplot executable ", false, "gnuplot", "filepath", cmd);
    TCLAP::ValueArg<double> tc_pixelsize("", "pixelsize", "Pixel size in microns. This also switches units to lp/mm", false, 1.0, "size", cmd);
    TCLAP::ValueArg<double> tc_lp1("", "lp1", "Lens profile resolution 1 (lp/mm or c/p)", false, 10.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp2("", "lp2", "Lens profile resolution 2 (lp/mm or c/p)", false, 30.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp3("", "lp3", "Lens profile resolution 3 (lp/mm or c/p)", false, 50.0, "lp/mm", cmd);
    TCLAP::ValueArg<string> tc_logfile("", "logfile", "Output written to <logfile> in stead of standard out", false, "", "filename", cmd);
    TCLAP::ValueArg<int> tc_gpwidth("", "gnuplot-width", "Width of images rendered by gnuplot", false, 1024, "pixels", cmd);
    TCLAP::ValueArg<double> tc_focal("", "focal-ratio", "Specify focal ratio for use in chart orientation estimation", false, -2, "ratio", cmd);
    #ifdef MDEBUG
    TCLAP::SwitchArg tc_single("","single-threaded","Force single-threaded operation", cmd, false);
    #endif

    vector<string> allowed_bayer_subsets;
    allowed_bayer_subsets.push_back("red");
    allowed_bayer_subsets.push_back("green");
    allowed_bayer_subsets.push_back("blue");
    allowed_bayer_subsets.push_back("none");
    TCLAP::ValuesConstraint<string> bayer_constraints(allowed_bayer_subsets);
    TCLAP::ValueArg<std::string> tc_bayer("", "bayer", "Select Bayer subset", false, "none", &bayer_constraints );
    cmd.add(tc_bayer);
    
    cmd.parse(argc, argv);

    if (tc_logfile.isSet()) {
        logger.redirect(tc_logfile.getValue(), tc_log_append.getValue());
    }
    if (tc_debug.getValue()) {
        logger.enable_level(Logger::DEBUG);
    }

    bool lpmm_mode = false;
    double pixel_size = 1;
    if (tc_pixelsize.isSet()) {
        logger.info("Info: Pixel size has been specified, measurement will be reported in lp/mm, rather than c/p\n");
        lpmm_mode = true;
        pixel_size = 1000 / tc_pixelsize.getValue();
        logger.info("working with=%lf pixels per mm\n", pixel_size);
    }

    if (!tc_profile.isSet() && !tc_annotate.isSet() && !tc_surface.isSet() && !tc_print.isSet() && !tc_sfr.isSet() && !tc_edges.isSet()) {
        logger.info("Warning: No output specified. You probably want to specify at least one of the following flags: [-r -p -a -s -f -q]\n");
    }

    cv::Mat cvimg;
	try {
	    cvimg = cv::imread(tc_in_name.getValue(),-1);
	} catch (const cv::Exception& ex) {
		cout << ex.what() << endl;
	}

	if (!cvimg.data) {
		logger.error("Fatal error: could not open input file <%s>.\nFile is missing, or not where you said it would be, or you do not have read permission.\n", tc_in_name.getValue().c_str());
		return -2;
	}
	
	struct STAT sb;
	if (STAT(tc_wdir.getValue().c_str(), &sb) != 0) {
	    logger.error("Fatal error: specified output directory <%s> does not exist\n", tc_wdir.getValue().c_str());
	    return -3;
	} else {
	    if (!S_ISDIR(sb.st_mode)) {
	        logger.error("Fatal error: speficied output directory <%s> is not a directory\n", tc_wdir.getValue().c_str());
	        return -3;
	    }
	}
    
    if (cvimg.type() == CV_8UC3 || cvimg.type() == CV_16UC3) {
        logger.info("colour input image detected; converting to grayscale using 0.299R + 0.587G + 0.114B\n");
        cv::Mat dest;
        cv::cvtColor(cvimg, dest, CV_RGB2GRAY);  // force to grayscale
        cvimg.release();
        cvimg = dest;
    }
    
    if (cvimg.type() == CV_8UC1) {
        logger.info("8-bit input image, upconverting %s\n", tc_linear.getValue() ? "with linear scaling" : "with sRGB gamma correction");
        convert_8bit_input(cvimg, !tc_linear.getValue());        
    } else {
        logger.info("16-bit input image, no upconversion required\n");
    }
   
    assert(cvimg.type() == CV_16UC1);

    const int border_width = 100;
    if (tc_border.getValue()) {
        logger.info("The -b option has been specified, adding a %d-pixel border to the image\n", border_width);
        double max_val = 0;
        double min_val = 0;
        cv::minMaxLoc(cvimg, &min_val, &max_val);
        cv::Mat border;
        cv::copyMakeBorder(cvimg, border, border_width, border_width, border_width, border_width, cv::BORDER_CONSTANT, cv::Scalar((int)max_val));
        cvimg = border;
    }
    
    int gnuplot_width = std::max(1024, tc_gpwidth.getValue());
    
    // process working directory
    std::string wdir(tc_wdir.getValue());
    if (wdir[wdir.length()-1]) {
        wdir = tc_wdir.getValue() + "/";
    }

    char slashchar='/';
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
	slashchar='\\';
	#endif
	
	// strip off supposed extention suffix,
	// and supposed path prefix
	int ext_idx=-1;
	int path_idx=0;
	for (int idx=tc_in_name.getValue().length()-1; idx >= 0 && path_idx == 0; idx--) {
	    if (tc_in_name.getValue()[idx] == '.' && ext_idx < 0) {
	        ext_idx = idx;
        }
	    if (tc_in_name.getValue()[idx] == slashchar && path_idx == 0) {
	        path_idx = idx;
	    }
    }
    if (ext_idx < 0) {
        ext_idx = tc_in_name.getValue().length();
    }
    std::string img_filename;
    for (int idx=path_idx; idx < ext_idx; idx++) {
        char c = tc_in_name.getValue()[idx];
        if (c == slashchar) continue;
        if (c == '_') {
            img_filename.push_back('\\');
            img_filename.push_back('\\');
        }
        img_filename.push_back(c);
    }
	
	// initialize the global mtf correction instance
	global_mtf_correction_instance = new Mtf_correction();

    cv::Mat masked_img;
    
    cv::Rect img_dimension_correction(0,0, cvimg.cols, cvimg.rows);
    
    if (tc_autocrop.getValue()) {
        Autocropper ac(cvimg);
        cvimg = ac.subset(cvimg, &img_dimension_correction);
    }
    
    cv::Mat rawimg = cvimg;
    if (tc_bayer.isSet()) {
        simple_demosaic(cvimg, rawimg, Bayer::from_string(tc_bayer.getValue()));
        //imwrite(string("prewhite.png"), rawimg);
        //imwrite(string("white.png"), cvimg);
    }
    
    size_t nthreads = std::thread::hardware_concurrency();
    ThreadPool tp (nthreads);
    
    logger.info("Thresholding image ...\n");
    int brad_S = tc_border.getValue() ? max(cvimg.cols, cvimg.rows) : max(500, min(cvimg.cols, cvimg.rows)/3);
    double brad_threshold = tc_thresh.getValue();
    #ifdef MDEBUG
        if (tc_bradley.getValue()) {
            printf("using Bradley thresholding\n");
            bradley_adaptive_threshold(cvimg, masked_img, brad_threshold, brad_S);
        } else {
            sauvola_adaptive_threshold(cvimg, masked_img, brad_threshold/0.55*0.85, brad_S, tp);
        }
    #else
        sauvola_adaptive_threshold(cvimg, masked_img, brad_threshold/0.55*0.85, brad_S, tp); // fudge the threshold to maintain backwards compatibility
    #endif
    
    logger.info("Computing gradients ...\n");
    Gradient gradient(cvimg);
    
    logger.info("Component labelling ...\n");
    Component_labeller::zap_borders(masked_img);    
    Component_labeller cl(masked_img, 60, false, 8000);

    if (cl.get_boundaries().size() == 0) {
        logger.error("Error: No black objects found. Try a lower threshold value with the -t option.\n");
        return 0;
    }
    
    // now we can destroy the thresholded image
    masked_img = cv::Mat(1,1, CV_8UC1);
    
    Mtf_core mtf_core(cl, gradient, cvimg, rawimg, tc_bayer.getValue());
    mtf_core.set_absolute_sfr(tc_absolute.getValue());
    mtf_core.set_sfr_smoothing(!tc_smooth.getValue());
    if (tc_border.getValue()) {
        logger.debug("setting border to %d\n", border_width);
        mtf_core.set_border(border_width+1);
    }
    
    if (tc_snap.isSet()) {
        mtf_core.set_snap_angle(tc_snap.getValue()/180*M_PI);
    }
    if (tc_focus.isSet() || tc_mf_profile.getValue()) {
        mtf_core.set_sliding(true);
        if (tc_mf_profile.getValue()) {
            mtf_core.set_samples_per_edge(5);
        }
    }
    if (tc_chart_orientation.isSet()) {
        mtf_core.set_find_fiducials(true);
    }
    
    Mtf_core_tbb_adaptor ca(&mtf_core);
    
    
    
    #ifdef MDEBUG
    if (tc_single.getValue()) {
        ca(Stride_range(size_t(0), mtf_core.num_objects()-1, 1));
    } else {
        logger.debug("Parallel MTF50 calculation\n");
        Stride_range::parallel_for(ca, tp, mtf_core.num_objects());
    }
    #else
    Stride_range::parallel_for(ca, tp, mtf_core.num_objects());
    #endif
    
    if (mtf_core.get_blocks().size() == 0) {
        logger.error("Error: No suitable target objects found.\n");
        return 0;
    }
    
    
    Distance_scale distance_scale;
    if (tc_mf_profile.getValue() || tc_focus.getValue() || tc_chart_orientation.getValue()) {
        distance_scale.construct(mtf_core, true, &img_dimension_correction, tc_focal.getValue());
    }
    
    // release most of the resources we no longer need
    masked_img.release();
    gradient.release();
    cl.release();
    
    // now render the computed MTF values
    if (tc_annotate.getValue()){
        Mtf_renderer_annotate annotate(cvimg, wdir + string("annotated.png"), lpmm_mode, pixel_size);
        annotate.render(mtf_core.get_blocks());
    }
    
    bool gnuplot_warning = true;
    bool few_edges_warned = false;
    
    if (tc_profile.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            logger.info("Warning: fewer than 10 edges found, so MTF50 surfaces/profiles will not be generated. Are you using suitable input images?\n");
            few_edges_warned = true;
        } else {
            Mtf_renderer_profile profile(
                img_filename,
                wdir, 
                string("profile.txt"),
                string("profile_peak.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                gnuplot_width,
                lpmm_mode,
                pixel_size
            );
            profile.render(mtf_core.get_blocks());
            gnuplot_warning = !profile.gnuplot_failed();
        }
    }
    
    if (tc_mf_profile.getValue()) {
        Mtf_renderer_mfprofile profile(
            distance_scale,
            wdir, 
            string("focus_peak.png"),
            cvimg,
            lpmm_mode,
            pixel_size
        );
        profile.render(mtf_core.get_samples());
    }

    if (tc_surface.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            if (!few_edges_warned) {
                logger.info("Warning: fewer than 10 edges found, so MTF50 surfaces/profiles will not be generated. Are you using suitable input images?\n");
            }
        } else {
            Mtf_renderer_grid grid(
                img_filename,
                wdir, 
                string("grid.txt"),
                tc_gnuplot.getValue(),
                cvimg,
                gnuplot_width,
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
    
    if (tc_edges.getValue()) {
        Mtf_renderer_edges printer(
            wdir + string("edge_mtf_values.txt"), 
            wdir + string("edge_sfr_values.txt"),
            lpmm_mode, pixel_size
        );
        printer.render(mtf_core.get_blocks());
    }
    
    if (tc_lensprof.getValue()) {
    
        vector<double> resolutions;
        // try to infer what resolutions the user wants
        if (!(tc_lp1.isSet() || tc_lp2.isSet() || tc_lp3.isSet())) {
            if (lpmm_mode) {
                // if nothing is specified explicitly, use the first two defaults
                resolutions.push_back(tc_lp1.getValue());
                resolutions.push_back(tc_lp2.getValue());
            } else {
                // otherwise just pick some arbitrary values
                resolutions.push_back(0.1);
                resolutions.push_back(0.2);
            }
        } else {
            if (tc_lp1.isSet()) {
                resolutions.push_back(tc_lp1.getValue());
            }
            if (tc_lp2.isSet()) {
                resolutions.push_back(tc_lp2.getValue());
            }
            if (tc_lp3.isSet()) {
                resolutions.push_back(tc_lp3.getValue());
            }
        }

        sort(resolutions.begin(), resolutions.end());
        
        Mtf_renderer_lensprofile printer(
            img_filename,
            wdir, 
            string("lensprofile.txt"),
            tc_gnuplot.getValue(),
            cvimg,
            resolutions,
            gnuplot_width,
            lpmm_mode,
            pixel_size
        );
        printer.render(mtf_core.get_blocks());
    }
    
    if (tc_chart_orientation.getValue()) {
        Mtf_renderer_chart_orientation co_renderer(
            img_filename,
            wdir, 
            string("chart_orientation.png"),
            cvimg,
            gnuplot_width,
            distance_scale,
            &img_dimension_correction
        );
        co_renderer.render(mtf_core.get_blocks());
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
    
    if (tc_focus.getValue()) {
        Mtf_renderer_focus profile(
            distance_scale,
            wdir, 
            string("focus_peak.png"),
            cvimg,
            lpmm_mode,
            pixel_size
        );
        profile.render(mtf_core.get_samples(), mtf_core.bayer, &mtf_core.ellipses, &img_dimension_correction);
    }
    
    Mtf_renderer_stats stats(lpmm_mode, pixel_size);
    if (tc_focus.getValue() || tc_mf_profile.getValue()) {
        stats.render(mtf_core.get_samples());
    } else {
        stats.render(mtf_core.get_blocks());
    }
    
    return 0;
}
