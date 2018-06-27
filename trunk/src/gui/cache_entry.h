#ifndef CACHE_ENTRY_H
#define CACHE_ENTRY_H

#include <opencv2/opencv.hpp>

class Cache_entry {
  public:
    Cache_entry(void): img(cv::Mat()), seq_num(0) {}
    Cache_entry(cv::Mat img): img(img), seq_num(counter()++) {}
    
    cv::Mat fetch(void) {
        seq_num = counter()++;
        return img;
    }
    
    uint32_t seq(void) const { 
        return seq_num;
    }
    
    uint64_t size(void) const {
        return uint64_t(img.total())*img.elemSize();
    }
    
  private:
    cv::Mat img;
    uint32_t seq_num;
    
    static uint32_t& counter(void) {
        static uint32_t global_counter = 1;
        return global_counter;
    }
};

#endif
