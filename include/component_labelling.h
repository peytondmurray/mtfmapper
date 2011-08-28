#ifndef COMPONENT_LABELLING_H
#define COMPONENT_LABELLING_H

#include "common_types.h"
#include <assert.h>
#include <string.h>

// A class for extracting the boundaries of a binary image.
//
// It expects black objects on white backgrounds.
//
// Implements the method in (F. Chang, C.-J. Chen, C.-J. Jen, A linear-time
// component labeling algorithm using contour tracing technique, Computer
// Vision and Image Understanding, 93:206-220, 2004)

//==============================================================================
class Component_labeller {
public:
    Component_labeller(void);
    Component_labeller(const Component_labeller& b);
    Component_labeller(const cv::Mat& in_img,
        int min_boundary_length = 10, bool snapshot = false, 
        int max_boundary_length = 5000);

    ~Component_labeller(void);

    Component_labeller& operator=(const Component_labeller& b);

    void configure(const cv::Mat& in_img,
        int min_boundary_length = 10, 
        int max_boundary_length = 5000,
        bool snapshot = false);


    const Boundarylist& get_boundaries(void) const {
        assert(configured);
        return _boundaries;
    }

    inline int operator[](int index) const {
        return _labels[index];
    }

    inline int operator()(int x, int y) const {
        if (!(x >= 0 && x < _width && y >= 0 && y < _height)) {
            printf("trying to access %d, %d (size is %d,%d)\n", x, y, _width, _height);
        }
        //assert(x >= 0 && x < _width && y >= 0 && y < _height);
        return _labels[y * _width + x];
    }

    inline int get_width(void) const {
        return _width;
    }

    inline int get_height(void) const {
        return _height;
    }
    
    // set borders to the background color, since objects that
    // touch the image borders cause havoc with the main algorithm
    static void zap_borders(cv::Mat& masked_img, int fill=255) {
        for (int y=0; y < masked_img.rows; y++) {
            masked_img.at<uchar>(y,0) = fill;
            masked_img.at<uchar>(y,1) = fill;
            masked_img.at<uchar>(y,masked_img.cols-1) = fill;
            masked_img.at<uchar>(y,masked_img.cols-2) = fill;
        }
        for (int x=0; x < masked_img.cols; x++) {
            masked_img.at<uchar>(0,x) = fill;
            masked_img.at<uchar>(1,x) = fill;
            masked_img.at<uchar>(masked_img.rows-1,x) = fill;
            masked_img.at<uchar>(masked_img.rows-2,x) = fill;
        }
    }

private:
    typedef enum {
        INTERNAL = 0,
        EXTERNAL = 1,
        EXTERNAL_FIRST = 2,
        INTERNAL_FIRST = 3
    } mode_type;

    void _find_components(void);

    void _contour_tracing(int x, int y, int label, mode_type mode);

    void _tracer(int x, int y, int& nx, int& ny,
        int& from, mode_type mode, bool mark_white = true);

    void _draw_snapshot(void);

    int _width;
    int _height;

    unsigned char* _pix;
    int* _labels;
    Boundarylist _boundaries;

    int _min_boundary_length;
    int _max_boundary_length;

    int C;

    bool configured;
};

#endif // COMPONENT_LABELLING_H


