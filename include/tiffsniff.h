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
#ifndef TIFFSNIFF_H
#define TIFFSNIFF_H

#include <stdio.h>
#include <string>
using std::string;

#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include "include/display_profile.h"

typedef enum {
    EXIF=1,
    ICC=2,
    NONE=15
} jpeg_app_t;


class Tiffsniff {
  public:
    Tiffsniff(const string& fname, bool is_8bit = false);
    ~Tiffsniff(void);
    bool profile_found(void) const { return has_profile; }
    Display_profile profile(void);
    
  private:
    typedef enum {
        sRGB,
        adobeRGB,
        UNKNOWN,
        CUSTOM
    } profile_t;
    
    typedef enum {
        ICC,
        TIFF,
        EXIF,
        EXIF_INTEROP
    } ifd_t;
   
    void parse_tiff(off_t offset) throw(int);
    void read_ifd(off_t offset, off_t base_offset = 0, ifd_t ifd_type = ifd_t::TIFF) throw(int);
    void read_icc_profile(off_t offset) throw(int);
    void read_trc_entry(off_t offset, uint32_t size) throw(int);
    vector<double> read_xyztype_entry(off_t offset, uint32_t size) throw(int);
    void read_curv_trc(off_t offset, uint32_t size) throw(int);
    void read_para_trc(off_t offset) throw(int);
    vector< pair<jpeg_app_t, off_t> > scan_jpeg_app_blocks(void) throw(int);
    double read_exif_gamma(off_t offset) throw(int);
    
    uint32_t read_uint32(void);
    uint16_t read_uint16(void);
    
    FILE* fin;
    bool big_endian = false;
    bool has_profile = false;
    bool assumed_sRGB = false;
    
    int exif_cs = -1;
    bool exif_gamma_found = false;
    bool exif_interop_r03 = false;
    
    profile_t inferred_profile = profile_t::UNKNOWN;
    off_t file_size;
    
    vector<double> gparm; // should also have the option of a vector<> ?
    vector< pair<uint16_t, uint16_t> > gtable;
    vector<double> luminance_weights;
};

typedef struct {
    uint16_t tag_id;
    uint16_t data_type;
    uint32_t data_count;
    uint32_t data_offset;
} tiff_field;

typedef struct {
    uint32_t tag_signature;
    uint32_t data_offset;
    uint32_t element_size;
    
	// we have to explicitly perform all the reads (fgets) first, or the optimizer may
	// re-arrange the arguments to the binary operators, thus performing the reads in
	// the wrong order
    static double read_fixed8_8(FILE* fin) {
		unsigned char b0 = fgetc(fin) & 0xff;
		unsigned char b1 = fgetc(fin) & 0xff;
        return double(b0) + double(b1)/256.0;
    }
    
    static uint16_t read_uint16(FILE* fin) {
		unsigned char b0 = fgetc(fin) & 0xff;
		unsigned char b1 = fgetc(fin) & 0xff;
        return (uint16_t(b0) << 8) | uint16_t(b1);
    }
    
    static double read_fixed15_16(FILE* fin) {
        char sb0 = fgetc(fin) & 0xff;
		unsigned char b1 = fgetc(fin) & 0xff;
		unsigned char b2 = fgetc(fin) & 0xff;
		unsigned char b3 = fgetc(fin) & 0xff;
        return double((int16_t(sb0) << 8) | b1) + 
            double(((uint16_t(b2) << 8) | b3))/65536.0;
    }
    
    static uint32_t read_uint32(FILE* fin) {
		unsigned char b0 = fgetc(fin) & 0xff;
		unsigned char b1 = fgetc(fin) & 0xff;
		unsigned char b2 = fgetc(fin) & 0xff;
		unsigned char b3 = fgetc(fin) & 0xff;
        return (uint32_t(b0) << 24) | (uint32_t(b1) << 16) |
            (uint32_t(b2) << 8) | uint32_t(b3);
    }
    
} icc_tag;

#endif
