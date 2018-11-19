#include "include/gradient.h"
#include "include/common_types.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//------------------------------------------------------------------------------
Gradient::Gradient(const cv::Mat& in_img, bool snapshot, bool thin)
 : _width(in_img.cols), _height(in_img.rows), _thin(thin)
{

    float* smoothed;
    _blur(in_img, &smoothed);

    _compute_gradients(smoothed, in_img.rows, in_img.cols, &_gradient_x, &_gradient_y);

    cv::Mat grad(in_img.rows, in_img.cols, CV_8UC1);

    _gradient_m = new float[in_img.cols * in_img.rows];
    for (int y=0; y < in_img.rows; y++) {
        for (int x=0; x < in_img.cols; x++) {
            int pos = x + y * in_img.cols;
            _gradient_m[pos] = sqrt(SQR(_gradient_x[pos]) + SQR(_gradient_y[pos]));
            grad.data[pos] = std::min(255, int(_gradient_m[pos]/65536.0*255));
        }
    }

    if (snapshot) {
        imwrite(std::string("gradient.png"), grad);
    }

    delete [] smoothed;
}

//------------------------------------------------------------------------------
Gradient::~Gradient(void) {
    delete [] _gradient_x;
    delete [] _gradient_y;
    delete [] _gradient_m;
}


//------------------------------------------------------------------------------
void Gradient::_blur(const cv::Mat& image, float **smoothedim) {
   const int& rows = image.rows;
   const int& cols = image.cols;

   int   windowsize;     /* Dimension of the gaussian kernel. */
   int   center;         /* Half of the windowsize. */
   double* tempim;       /* Buffer for separable filter gaussian smoothing. */

   // fixed gaussian kernel, as used by circle detector
   double gkernel[5] = {0.1353, 0.6065, 1, 0.6065, 0.1353}; // sigma = 1
   windowsize = 5;
   center = windowsize / 2;

   // Allocate a temporary buffer image and the smoothed image.
   if ((tempim = new double[rows*cols]) == NULL){
      fprintf(stderr, "Error allocating the buffer image.\n");
      exit(1);
   }
   memset(tempim, 0, rows*cols*sizeof(double));
   if ( ((*smoothedim) = new float[rows*cols]) == NULL){
      fprintf(stderr, "Error allocating the smoothed image.\n");
      exit(1);
   }
   memset(*smoothedim, 0, rows*cols*sizeof(float));

   // Blur in the x - direction.
   for(int r=0;r<rows;r++){
      for(int c=0;c<cols;c++){
         double dot = 0.0;
         double sum = 0.0;
         for(int cc=(-center);cc<=center;cc++){
            if(((c+cc) >= 0) && ((c+cc) < cols)){
               dot += image.at<uint16_t>(r, c+cc) * gkernel[center+cc];
               sum += gkernel[center+cc];
            }
         }
         tempim[r*cols+c] = dot/sum;
      }
   }

   // Blur in the y - direction.
   for(int c=0;c<cols;c++){
      for(int r=0;r<rows;r++){
         double sum = 0.0;
         double dot = 0.0;
         for(int rr=(-center);rr<=center;rr++){
            if(((r+rr) >= 0) && ((r+rr) < rows)){
               dot += tempim[(r+rr)*cols+c] * gkernel[center+rr];
               sum += gkernel[center+rr];
            }
         }
         (*smoothedim)[r*cols+c] = (float)dot/sum;
      }
   }


   delete [] tempim;
}

//------------------------------------------------------------------------------
void Gradient::_compute_gradients(const float* smoothed_im, int rows, int cols, float** grad_x, float** grad_y) {

    if ((*grad_x = new float[rows*cols]) == NULL){
       fprintf(stderr, "Error allocating the x-gradient buffer image.\n");
       exit(1);
    }
    memset(*grad_x, 0, sizeof(float)*rows*cols);
    if ((*grad_y = new float[rows*cols]) == NULL){
       fprintf(stderr, "Error allocating the y-gradient buffer image.\n");
       exit(1);
    }
    memset(*grad_y, 0, sizeof(float)*rows*cols);

    for (int y=0; y < rows; y++) {
        (*grad_y)[y*cols] = 0;
        (*grad_y)[cols-1 + y*cols] = 0;
    }
    for (int x=0; x < cols; x++) {
        (*grad_x)[x] = 0;
        (*grad_x)[x + (rows-1)*cols] = 0;
    }


    if (_thin) {
        for(int r=0;r<rows;r++){
            for(int c=1;c<cols-1;c++){
                (*grad_x)[r*cols+c] = smoothed_im[r*cols + c ] - smoothed_im[r*cols + c - 1];
            }
        }
        for(int r=1;r<rows-1;r++){
            for(int c=0;c<cols;c++){
                (*grad_y)[r*cols+c] = smoothed_im[r*cols + c ] - smoothed_im[r*cols + c - cols];
            }
        }
    } else {
        for(int r=0;r<rows;r++){
            for(int c=1;c<cols-1;c++){
                (*grad_x)[r*cols+c] = smoothed_im[r*cols + c + 1] - smoothed_im[r*cols + c - 1];
            }
        }
        for(int r=1;r<rows-1;r++){
            for(int c=0;c<cols;c++){
                (*grad_y)[r*cols+c] = smoothed_im[r*cols + c + cols] - smoothed_im[r*cols + c - cols];
            }
        }
    }


}

