#ifndef LOESS_FIT_H
#define LOESS_FIT_H

#include <vector>
using std::vector;

#include "include/common_types.h"

class Ordered_point {
public:
    Ordered_point(double in_first=0, double in_second=0) : first(in_first), second(in_second) {}
    bool operator< (const Ordered_point& b) const {
        return first < b.first;
    }
    
    double first;
    double second;
};

double loess_core(vector<Ordered_point>& ordered, size_t start_idx, size_t end_idx,
    double mid,  Point& sol, int mode=1);
    
void loess_fit(vector< Ordered_point  >& ordered, double* fft_in_buffer, const int fft_size, bool deriv=true);

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#endif // LOESS_FIT_H
