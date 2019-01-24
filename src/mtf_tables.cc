#include "include/mtf_tables.h"
#include "include/afft.h"
#include "include/fastexp.h"

#include <complex>


typedef std::complex<double> cplex;

inline double sinc(double x) {
    return x == 0 ? 1 : sin(x)/x;
}

Mtf_correction::Mtf_correction(void) 
: w(NYQUIST_FREQ*4, 0.0), sdev(13.0), nbins(16) {
    update_correction_tables();
}

void Mtf_correction::update_correction_tables(void) {
    
    constexpr int c_fft_size = 512*16;
    vector<double> fft_buffer(c_fft_size);
    AFFT<c_fft_size> afft;
    
    //FILE* fc = fopen("kernel.txt", "wt");
    
    for (int i=0; i < c_fft_size; i++) {
        double x = (i - c_fft_size/2)/(16*8.0);
        fft_buffer[i] = evaluate(x);
        //fprintf(fc, "%le %le\n", x, fft_buffer[i]);
    }
    
    //fclose(fc);
    
    afft.realfft(fft_buffer.data());
    
    //fc = fopen("correction.txt", "wt");
    
    vector<double> correction(w.size());
    double n0 = fabs(fft_buffer[0]);
    for (int i=0; i < NYQUIST_FREQ*4; i++) {
        double dc_x = 2*M_PI*i/double(NYQUIST_FREQ*2);
        correction[i] = sqrt(SQR(fft_buffer[i]) + SQR(fft_buffer[c_fft_size - i])) / n0;
        //fprintf(fc, "%le %le\n", dc_x, correction[i]);
    }

    //fclose(fc);
    
    w[0] = 1.0;
    for (int i=1; i < NYQUIST_FREQ*4; i++) {
        double dc_x = 2*M_PI*i/double(NYQUIST_FREQ*2);
        w[i] = (sinc(dc_x/8.0)) * correction[i];
    }
}

double Mtf_correction::evaluate(double x) {
    if (fabs(x) > nbins*0.125) return 0;
    return fastexp(-sdev*x*x);
}
