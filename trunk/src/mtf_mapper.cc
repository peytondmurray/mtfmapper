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
#include "include/mtf_renderer_mfprofile.h"
#include "include/mtf_renderer_grid.h"
#include "include/mtf_renderer_print.h"
#include "include/mtf_renderer_stats.h"
#include "include/mtf_renderer_sfr.h"
#include "include/mtf_renderer_esf.h"
#include "include/mtf_renderer_edges.h"
#include "include/mtf_renderer_lensprofile.h"
#include "include/mtf_tables.h"
#include "include/scanline.h"
#include "include/ellipse_decoder.h"
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
    TCLAP::SwitchArg tc_mf_profile("","mf-profile","Generate MTF50 profile (Manual focus chart)", cmd, false);
    TCLAP::SwitchArg tc_annotate("a","annotate","Annotate input image with MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_surface("s","surface","Generate MTF50 surface plots", cmd, false);
    TCLAP::SwitchArg tc_linear("l","linear","Input image is linear 8-bit (default for 8-bit is assumed to be sRGB gamma corrected)", cmd, false);
    TCLAP::SwitchArg tc_print("r","raw","Print raw MTF50 values", cmd, false);
    TCLAP::SwitchArg tc_edges("q","edges","Print raw MTF50 values, grouped by edge location", cmd, false);
    TCLAP::SwitchArg tc_sfr("f","sfr","Store raw SFR curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_esf("e","esf","Store raw ESF and PSF curves for each edge", cmd, false);
    TCLAP::SwitchArg tc_lensprof("","lensprofile","Render M/S lens profile plot", cmd, false);
    TCLAP::SwitchArg tc_border("b","border","Add a border of 20 pixels to the image", cmd, false);
    TCLAP::SwitchArg tc_absolute("","absolute-sfr","Generate absolute SFR curve (MTF) i.s.o. relative SFR curve", cmd, false);
    TCLAP::SwitchArg tc_smooth("","nosmoothing","Disable SFR curve (MTF) smoothing", cmd, false);
    TCLAP::ValueArg<double> tc_angle("g", "angle", "Angular filter [0,360)", false, 0, "angle", cmd);
    TCLAP::ValueArg<double> tc_snap("", "snap-angle", "Snap-to angle modulus [0,90)", false, 1000, "angle", cmd);
    TCLAP::ValueArg<double> tc_thresh("t", "threshold", "Dark object threshold (0,1)", false, 0.75, "threshold", cmd);
    TCLAP::ValueArg<string> tc_gnuplot("", "gnuplot-executable", "Full path (including filename) to gnuplot executable ", false, "gnuplot", "filepath", cmd);
    TCLAP::ValueArg<double> tc_pixelsize("", "pixelsize", "Pixel size in microns. This also switches units to lp/mm", false, 1.0, "size", cmd);
    TCLAP::ValueArg<double> tc_lp1("", "lp1", "Lens profile resolution 1 (lp/mm or c/p)", false, 10.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp2("", "lp2", "Lens profile resolution 2 (lp/mm or c/p)", false, 30.0, "lp/mm", cmd);
    TCLAP::ValueArg<double> tc_lp3("", "lp3", "Lens profile resolution 3 (lp/mm or c/p)", false, 50.0, "lp/mm", cmd);

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

    if (!tc_profile.isSet() && !tc_annotate.isSet() && !tc_surface.isSet() && !tc_print.isSet() && !tc_sfr.isSet() && !tc_edges.isSet()) {
        printf("Warning: No output specified. You probably want to specify at least one of the following flags: [-r -p -a -s -f -q]\n");
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
	
	struct STAT sb;
	if (STAT(tc_wdir.getValue().c_str(), &sb) != 0) {
	    printf("Fatal error: specified output directory <%s> does not exist\n", tc_wdir.getValue().c_str());
	    return -3;
	} else {
	    if (!S_ISDIR(sb.st_mode)) {
	        printf("Fatal error: speficied output directory <%s> is not a directory\n", tc_wdir.getValue().c_str());
	        return -3;
	    }
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
    
    Point zero;
    Point transverse;
    Point longitudinal;
    double chart_scale=1.0;
    vector<cv::Point3d> distance_scale;
    if (tc_mf_profile.getValue()) {
    
        // TODO: move all this out of main program
        
        // define reference frame from three solid circular dots in centre of chart
        printf("#solid ellipses: %d\n", (int)mtf_core.solid_ellipses.size());
        for (size_t j=0; j < mtf_core.solid_ellipses.size(); j++) {
            printf("(%lf %lf)\n", mtf_core.solid_ellipses[j].x, mtf_core.solid_ellipses[j].y);
        }
        printf("\n\n");
        if (mtf_core.solid_ellipses.size() >= 3) {
            printf("got all three solid ellipses\n");
            double maxdelta = 0;
            Point first;
            Point last;
            for (int i=0; i <= 1; i++) {
                for (int j=i+1; j < 3; j++) {
                    Point d = mtf_core.solid_ellipses[i] - mtf_core.solid_ellipses[j];
                    double dist = sqrt(d.x*d.x + d.y*d.y);
                    printf("considering (%lf, %lf) vs (%lf, %lf), dist=%lf ",
                        mtf_core.solid_ellipses[i].x, mtf_core.solid_ellipses[i].y,
                        mtf_core.solid_ellipses[j].x, mtf_core.solid_ellipses[j].y,
                        dist
                    );
                    if (dist > maxdelta) {
                        first = mtf_core.solid_ellipses[i];
                        last  = mtf_core.solid_ellipses[j];
                        maxdelta = dist;
                        printf("<- keeping this one\n");
                    } else printf("\n");
                }
            }
            Point base = first - last;
            chart_scale = 91.0 / sqrt(base.x*base.x + base.y*base.y);
            zero = Point(0.5*(first.x+last.x), 0.5*(first.y+last.y));
            printf("absolute centre: (%lf, %lf)\n", zero.x, zero.y);
            transverse = normalize(first - last);
            printf("transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
            // sign of transverse unknown at this stage
            longitudinal = Point(-transverse.y, transverse.x);
            printf("longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
            printf("chart scale is %lf mm per pixels\n", chart_scale);
            
            // now try to decode the ellipses
            vector< pair<double, int> > dist_by_idx;
            for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                //Ellipse_decoder ed(mtf_core.ellipses[i], mtf_core.cl, transverse);
                Ellipse_decoder ed(mtf_core.ellipses[i], mtf_core.img, transverse);
                mtf_core.ellipses[i].set_code(ed.code);
                if (!mtf_core.ellipses[i].solid) {
                    Point d(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
                    dist_by_idx.push_back(make_pair(sqrt(d.x*d.x + d.y*d.y), i));
                }
            }
            sort(dist_by_idx.begin(), dist_by_idx.end());
            int up = dist_by_idx[0].second;
            int down = dist_by_idx[1].second;
            if (mtf_core.ellipses[up].code > 3 || mtf_core.ellipses[down].code > 3) {
                for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                    mtf_core.ellipses[i].set_code(Ellipse_decoder::reverse(mtf_core.ellipses[i].code));
                }
                transverse = -transverse;
                printf("flipping transverse direction\n");
            }
            if (mtf_core.ellipses[up].code > 3 || mtf_core.ellipses[down].code > 3) {
                printf("Warning: Ellipse targets closest to zero appear to have incorrect codes %d, %d\n",
                    mtf_core.ellipses[up].code, mtf_core.ellipses[down].code
                );
            } else {
                if (mtf_core.ellipses[up].code != 1) {
                    std::swap(up, down);
                }
                if (mtf_core.ellipses[up].code == 1) {
                    printf("up ellipse identified\n");
                } else {
                    printf("Warning: could not identify 'up' ellipse (codes are up=%d, down=%d)\n", 
                        mtf_core.ellipses[up].code, mtf_core.ellipses[down].code
                    );
                }
                Point dv(mtf_core.ellipses[up].centroid_x - zero.x, mtf_core.ellipses[up].centroid_y - zero.y);
                double dot = longitudinal.x * dv.x + longitudinal.y * dv.y;
                if (dot < 0) {
                    longitudinal = -longitudinal;
                }
                printf("Oriented transverse vector: (%lf, %lf)\n", transverse.x, transverse.y);
                printf("Oriented longitudinal vector: (%lf, %lf)\n", longitudinal.x, longitudinal.y);
                
                // project ellipses onto longitudinal vector
                dist_by_idx.clear();
                for (size_t i=0; i < mtf_core.ellipses.size(); i++) {
                    Point delta(mtf_core.ellipses[i].centroid_x - zero.x, mtf_core.ellipses[i].centroid_y - zero.y);
                    double dot = longitudinal.x * delta.x + longitudinal.y * delta.y;
                    dist_by_idx.push_back(make_pair(dot, i));
                }
                sort(dist_by_idx.begin(), dist_by_idx.end());
                bool centre_done = false;
                for (size_t i=0; i < dist_by_idx.size(); i++) {
                    int j=dist_by_idx[i].second;
                    printf("(%lf, %lf), dist=%lf -> code=%d\n",
                        mtf_core.ellipses[j].centroid_x, mtf_core.ellipses[j].centroid_y,
                        dist_by_idx[i].first,
                        mtf_core.ellipses[j].code
                    );
                    if (fabs(dist_by_idx[i].first) < 5) {
                        if (!centre_done) { 
                            distance_scale.push_back(cv::Point3d(0,0,0));
                            centre_done = true;
                        }
                    } else {
                        distance_scale.push_back(
                            cv::Point3d(
                                dist_by_idx[i].first, 
                                (((dist_by_idx[i].first < 0) ? -2 : 0) + mtf_core.ellipses[j].code) * (1.5*5) * 
                                 (dist_by_idx[i].first < 0 ? -1 : 1),
                                mtf_core.ellipses[j].minor_axis/mtf_core.ellipses[j].major_axis
                            )
                        );
                    }
                }
            }
            
            // construct a distance scale
            
        } else {
            printf("Warning: Manual focus profile chart mode requested, but central dots not found\n");
            const vector<Block>& blocks = mtf_core.get_blocks();
            // find largest block
            vector< pair<double, int> > by_size;
            for (size_t i=0; i < blocks.size(); i++) {
                by_size.push_back(make_pair(blocks[i].get_area(), i));
            }
            sort(by_size.begin(), by_size.end());
            double delta_1 = by_size[by_size.size()-1].first / by_size[by_size.size()-2].first;
            double delta_2 = by_size[by_size.size()-2].first / by_size[by_size.size()-3].first;
            if (delta_1/delta_2 > 100) {
                const Block& lblock = blocks[by_size.back().second];
                // we have a clear largest block. now determine its orientation
                vector<double> xp;
                vector<double> yp;
                for (size_t i=0; i < blocks.size(); i++) {
                    xp.push_back(blocks[i].centroid.x);
                    yp.push_back(blocks[i].centroid.y);
                }
                sort(xp.begin(), xp.end());
                sort(yp.begin(), yp.end());
                double idx_x = (find(xp.begin(), xp.end(), lblock.centroid.x) - xp.begin())/double(xp.size());
                double idx_y = (find(yp.begin(), yp.end(), lblock.centroid.y) - yp.begin())/double(yp.size());
                printf("xfrac=%lf, yfrac=%lf\n", idx_x, idx_y);
                
                Point median(xp[xp.size()/2], yp[yp.size()/2]); 
                
                // orientations relative to chart, not image
                int top_i=lblock.get_edge_index(Block::TOP);
                int bot_i=lblock.get_edge_index(Block::BOTTOM);
                int left_i=lblock.get_edge_index(Block::LEFT);
                int right_i=lblock.get_edge_index(Block::RIGHT);
                if (fabs(idx_x - 0.5) < fabs(idx_y - 0.5)) {
                    // outer rows arranged in columns
                    if (fabs(lblock.get_edge_centroid(top_i).y - median.y) <
                        fabs(lblock.get_edge_centroid(bot_i).y - median.y) ) {
                        
                        // chart is upside down
                        std::swap(top_i, bot_i);
                        std::swap(left_i, right_i);
                    }
                } else {
                    // outer rows arranged in rows
                    std::swap(bot_i, right_i);
                    std::swap(top_i, left_i);
                    
                    if (fabs(lblock.get_edge_centroid(top_i).x - median.x) <
                        fabs(lblock.get_edge_centroid(bot_i).x - median.x) ) {
                        
                        // chart is upside down
                        std::swap(top_i, bot_i);
                        std::swap(left_i, right_i);
                    }
                }
                transverse = normalize(lblock.get_edge_centroid(right_i) - lblock.get_edge_centroid(left_i));
                longitudinal = normalize(lblock.get_edge_centroid(bot_i) - lblock.get_edge_centroid(top_i));
                zero = lblock.get_edge_centroid(bot_i);
                chart_scale = 0.15;
                printf("Warning: choosing (potentially) poor chart scale of %lf mm/pixel\n", chart_scale);
            } else { 
                printf("Warning: Could not identify largest block, choosing poor defaults\n");
                chart_scale = 0.15;
                zero = Point(932, 710);
                transverse = Point(1,0);
                longitudinal = Point(0,1);
            }
            printf("zero: %lf %lf\n", zero.x, zero.y);
            printf("transverse: %lf %lf\nlongitudinal: %lf %lf\n", transverse.x, transverse.y, longitudinal.x, longitudinal.y);
            
        }
    }
    
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
    
    if (tc_mf_profile.getValue()) {
        Mtf_renderer_mfprofile profile(
            zero,
            transverse,
            longitudinal,
            chart_scale, 
            wdir, 
            string("focus_peak.png"),
            cvimg,
            distance_scale,
            lpmm_mode,
            pixel_size
        );
        profile.render(mtf_core.get_blocks());
    }

    if (tc_surface.getValue()) {
        if (mtf_core.get_blocks().size() < 10) {
            if (!few_edges_warned) {
                printf("Warning: fewer than 10 edges found, so MTF50 surfaces/profiles will not be generated. Are you using suitable input images?\n");
            }
        } else {
            Mtf_renderer_grid grid(
                img_filename,
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
