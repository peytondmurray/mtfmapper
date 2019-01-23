#include "include/mtf_tables.h"

Mtf_correction::Mtf_correction(void) 
: w(NYQUIST_FREQ*4, 0.0) {
    w[0] = 1.0;
    for (int i=1; i < NYQUIST_FREQ*4; i++) {
        double dc_x = 2*M_PI*i/double(NYQUIST_FREQ*2);
        w[i] = sinc(dc_x/32.0);
    }
}
