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

#include "include/tiffsniff.h"
#include "include/logger.h"
#include "include/common_types.h"
#include <string.h>

Tiffsniff::Tiffsniff(const string& fname) : gparm(7, 0.0), luminance_weights(3, 0.0) {
    // linear gamma, by default
    gparm[0] = 1.0;
    gparm[1] = 1.0;
    // use sRGB weights, by default
    luminance_weights[0] = 0.2126;
    luminance_weights[1] = 0.7152;
    luminance_weights[2] = 0.0722;

    fin = fopen(fname.c_str(), "rb");
    if (fin) {
        if (fseek(fin, 0, SEEK_END) == 0) {
            file_size = ftell(fin);
        }
        rewind(fin);
        
        try {
        
            unsigned char magic[4];
            size_t nread = fread(magic, 1, 4, fin);
        
            if ((magic[0] == 0x4d && magic[1] == 0x4d) ||
                (magic[0] == 0x49 && magic[1] == 0x49)) { // possible TIFF file
                
                parse_tiff(0);
            }
            
            if (magic[0] == 0xff && magic[1] == 0xd8) { // possible JPEG file
                fseek(fin, -2, SEEK_END);
                if (fgetc(fin) == 0xff && fgetc(fin) == 0xd9) {
                    logger.debug("Valid JPEG file detected.\n");
                    has_profile = true;
                    // now build a list of APP blocks
                    auto blocks = scan_jpeg_app_blocks();
                    
                    if (blocks.size() == 0) {
                        // this is a raw JPEG file
                        assumed_sRGB = true;
                    } else {
                        bool icc_found = false;
                        // if we have an ICC block, use that
                        for (size_t i=0; i < blocks.size() && !icc_found; i++) {
                            if (blocks[i].first == jpeg_app_t::ICC) {
                                icc_found = true;
                                read_icc_profile(blocks[i].second);
                            }
                        }
                        
                        bool exif_found = false;
                        if (!icc_found) { // no ICC block, settle for EXIF
                            for (size_t i=0; i < blocks.size() && !exif_found; i++) {
                                if (blocks[i].first == jpeg_app_t::EXIF) {
                                    exif_found = true;
                                    parse_tiff(blocks[i].second);
                                    
                                    if (confirmed_sRGB) {
                                        logger.debug("sRGB colour space confirmed by EXIF data\n");
                                    } 
                                }
                            }
                        }
                        
                        assumed_sRGB = true;
                    }
                }
            }
        } catch (int) {
            logger.error("Error while parsing input image, unable to extract colourspace information.\n");
        }
    }
    
    if (confirmed_sRGB || assumed_sRGB) {
        gparm[0] = 2.4;
        gparm[1] = 1.0/1.055;
        gparm[2] = 0.055/1.055;
        gparm[3] = 1.0/12.92;
        gparm[4] = 0.04045;
        
        luminance_weights[0] = 0.2126;
        luminance_weights[1] = 0.7152;
        luminance_weights[2] = 0.0722;
    }
}

Tiffsniff::~Tiffsniff(void) {
    if (fin) {
        fclose(fin);
    }
}

// we have to explicitly perform all the reads (fgets) first, or the optimizer may
// re-arrange the arguments to the binary operators, thus performing the reads in
// the wrong order
uint32_t Tiffsniff::read_uint32(void) {
    if (big_endian) {
        return icc_tag::read_uint32(fin);
    } 
	unsigned char b0 = fgetc(fin) & 0xff;
	unsigned char b1 = fgetc(fin) & 0xff;
	unsigned char b2 = fgetc(fin) & 0xff;
	unsigned char b3 = fgetc(fin) & 0xff;
    return uint32_t(b0) | (uint32_t(b1) << 8) |
        (uint32_t(b2) << 16) | (uint32_t(b3) << 24);
}

uint16_t Tiffsniff::read_uint16(void) {
    if (big_endian) {
        return icc_tag::read_uint16(fin);
    } 
	unsigned char b0 = fgetc(fin) & 0xff;
	unsigned char b1 = fgetc(fin) & 0xff;
    return uint32_t(b0) | (uint32_t(b1) << 8);
}

Display_profile Tiffsniff::profile(void) {
    if (gtable.size() == 0) {
        return Display_profile(gparm, luminance_weights);
    }
    return Display_profile(gtable, luminance_weights);
}

void Tiffsniff::parse_tiff(off_t offset) throw(int) {
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        const char be_id[4] = {0x4D, 0x4D, 0x00, 0x2A};
        const char le_id[4] = {0x49, 0x49, 0x2A, 0x00};
        
        unsigned char magic[4];
        size_t nread = fread(magic, 1, 4, fin);
        
        if (nread == 4) {
            // TIFF
            bool is_valid_tiff = false;
            if (memcmp(le_id, magic, 4) == 0) {
                if (offset == 0) {
                    logger.debug("Little endian TIFF detected\n");
                } else {
                    logger.debug("EXIF / Little endian TIFF detected\n");
                }
                is_valid_tiff = true;
            }
            if (memcmp(be_id, magic, 4) == 0) {
                if (offset == 0) {
                    logger.debug("Big endian TIFF detected\n");
                } else {
                    logger.debug("EXIF / Big endian TIFF detected\n");
                }
                is_valid_tiff = true;
                big_endian = true;
            }
            
            if (is_valid_tiff) {
                uint32_t ifd_offset = read_uint32();
                read_ifd(ifd_offset + offset, offset);
                has_profile = true;
            } 
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_ifd(off_t offset, off_t base_offset) throw(int) {
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        uint16_t cur_ifd_entries = read_uint16();
        
        // check for weird values; realistic tiff files will not have thousands of entries
        if (cur_ifd_entries == 0 || cur_ifd_entries > 32765 || feof(fin)) {
            throw -1;
        }
        
        tiff_field field;
        for (uint16_t i=0; i < cur_ifd_entries; i++) {
            field.tag_id = read_uint16();
            field.data_type = read_uint16();
            field.data_count = read_uint32();
            field.data_offset = read_uint32();;
            
            if (i < (cur_ifd_entries-1) && feof(fin)) {
                throw -1;
            }
            
            if (field.tag_id == 0x8773) { // ICC profile 
                long fpos = ftell(fin);
                read_icc_profile(field.data_offset);
                fseek(fin, fpos, SEEK_SET);
            }
            
            if (field.tag_id == 0x8769) { // EXIF field
                long fpos = ftell(fin);
                read_ifd(base_offset + field.data_offset, base_offset);
                fseek(fin, fpos, SEEK_SET);
            }
            
            if (field.tag_id == 0xa001) { // EXIF colourspace tag
                confirmed_sRGB = (field.data_offset / 65536) == 1;
            }
        }
        
        // read next IDF offset
        uint32_t next_offset = read_uint32();
        if (next_offset) {
            if (next_offset > file_size || feof(fin)) {
                throw -1;
            }
            read_ifd(base_offset + next_offset, base_offset);
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_icc_profile(off_t offset) throw(int) {
    char icc_header[128];
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        if (fread(icc_header, 1, 128, fin) != 128) {
            throw -1;
        }
        logger.debug("Device: %c%c%c%c\n", icc_header[12], icc_header[13], icc_header[14], icc_header[15]);
        logger.debug("Colour space: %c%c%c%c\n", icc_header[16], icc_header[17], icc_header[18], icc_header[19]);
        
        // looks like ICC profiles are stored in big endian format
        uint32_t icc_entries = icc_tag::read_uint32(fin);
        if (icc_entries == 0 || icc_entries > 32765 || feof(fin)) {
            throw -1;
        }
        
        bool found_trc = false;
        
        icc_tag tag;
        for (uint32_t i=0; i < icc_entries; i++) {
            tag.tag_signature = icc_tag::read_uint32(fin);
            tag.data_offset = icc_tag::read_uint32(fin);
            tag.element_size = icc_tag::read_uint32(fin);
            
            const off_t icc_tag_offset = offset + tag.data_offset;
            if (icc_tag_offset > file_size || feof(fin)) {
                throw -1;
            }
            
            long fpos = ftell(fin); // before we jump to an individual entry
            
            // Just grab the first TRC get find, since they are probably all the same anyway
            if (!found_trc && (tag.tag_signature == 0x67545243 || tag.tag_signature == 0x6B545243 ||
                tag.tag_signature == 0x62545243 || tag.tag_signature == 0x72545243)) {
                read_trc_entry(icc_tag_offset, tag.element_size);
                found_trc = true;
            }
            
            if (tag.tag_signature == 0x7258595A) {
                vector<double> col = read_matrix_column_entry(icc_tag_offset, tag.element_size);
                luminance_weights[0] = col[1];
            }
            
            if (tag.tag_signature == 0x6758595a) {
                vector<double> col = read_matrix_column_entry(icc_tag_offset, tag.element_size);
                luminance_weights[1] = col[1];
            }
            
            if (tag.tag_signature == 0x6258595A) {
                vector<double> col = read_matrix_column_entry(icc_tag_offset, tag.element_size);
                luminance_weights[2] = col[1];
            }
            
            if (tag.tag_signature == 0x77747074) {
                read_xyztype_entry(icc_tag_offset, tag.element_size);
            }
            
            
            
            fseek(fin, fpos, SEEK_SET); // return to ICC IFD
        }
        if (!found_trc) {
            logger.error("Warning: image contains ICC profile, but no TRC curve was found.\n");
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_trc_entry(off_t offset, uint32_t size) throw(int) {
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        char trc_type[4];
        int nread = fread(&trc_type, 1, 4, fin);
        
        if (strncmp(trc_type, "curv", 4) == 0) {
            read_curv_trc(offset);
        } else {
            if (strncmp(trc_type, "para", 4) == 0) {
                read_para_trc(offset);
            }
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_curv_trc(off_t offset) throw(int) {
    int seekerr = fseek(fin, offset + 8, SEEK_SET);
    if (!seekerr) {
        uint32_t ecount = icc_tag::read_uint32(fin);
        switch(ecount) {
        case 0: // nothing to do here
            break;
        case 1:
            gparm[0] = icc_tag::read_fixed8_8(fin);
            logger.debug("ICC 'curv' gamma is %lf\n", gparm[0]);
            break;
        default: // with two or more entries, use the table
            gtable = vector< pair<uint16_t, uint16_t> >(ecount);
            for (int i=0; i < ecount; i++) {
                gtable[i].first = i*65535/(ecount-1);
                gtable[i].second = icc_tag::read_uint16(fin);
            }
            logger.debug("ICC 'curv' gamma table with %ld entries\n", gtable.size());
            break;
        }
    } else {
        throw -1;
    }
}

void Tiffsniff::read_para_trc(off_t offset) throw(int) {
    int seekerr = fseek(fin, offset + 8, SEEK_SET);
    if (!seekerr) {
        uint16_t ftype = icc_tag::read_uint16(fin);
        icc_tag::read_uint16(fin); // dump the reserved bytes
        
        // See Table 70 of ICC 1:2001-12 (ICC profile v4)
        switch(ftype) {
        case 0:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = 1.0;
            break;
        case 1:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = 0;
            gparm[4] = -gparm[2]/gparm[1];
            break;
        case 2:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = 0;
            gparm[4] = -gparm[2]/gparm[1];
            gparm[6] = icc_tag::read_fixed15_16(fin);
            break;
        case 3: // sRGB, supposedly
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = icc_tag::read_fixed15_16(fin);
            gparm[4] = icc_tag::read_fixed15_16(fin);
            break;
        case 4:
            gparm[0] = icc_tag::read_fixed15_16(fin);
            gparm[1] = icc_tag::read_fixed15_16(fin);
            gparm[2] = icc_tag::read_fixed15_16(fin);
            gparm[3] = icc_tag::read_fixed15_16(fin);
            gparm[4] = icc_tag::read_fixed15_16(fin);
            gparm[5] = icc_tag::read_fixed15_16(fin);
            gparm[6] = icc_tag::read_fixed15_16(fin);
            break;
        default: // treat this as linear
            logger.debug("unknown ICC parametricCurveType %d\n", ftype);
            break;
        }
    } else {
        throw -1;
    }
}

vector<double> Tiffsniff::read_matrix_column_entry(off_t offset, uint32_t size) throw(int) {
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        uint32_t sig = icc_tag::read_uint32(fin);
        icc_tag::read_uint32(fin); // drop reserved field
        
        double x = icc_tag::read_fixed15_16(fin);
        double y = icc_tag::read_fixed15_16(fin);
        double z = icc_tag::read_fixed15_16(fin);
        
        return vector<double>{x,y,z};
    } else {
        throw -1;
    }
    return vector<double>{0,0,0};
}

vector<double> Tiffsniff::read_xyztype_entry(off_t offset, uint32_t size) throw(int) {
    int seekerr = fseek(fin, offset, SEEK_SET);
    if (!seekerr) {
        uint32_t sig = icc_tag::read_uint32(fin);
        icc_tag::read_uint32(fin); // drop reserved field
        
        int n = (size - 8)/12;
        
        double x = icc_tag::read_fixed15_16(fin);
        double y = icc_tag::read_fixed15_16(fin);
        double z = icc_tag::read_fixed15_16(fin);
        
        logger.debug("mediaWhitepoint : x=%lf, y=%lf, z=%lf\n", x, y, z);
        
        return vector<double>{x,y,z};
    } else {
        throw -1;
    }
    return vector<double>{0,0,0};
}

vector< pair<jpeg_app_t, off_t> > Tiffsniff::scan_jpeg_app_blocks(void) throw(int) {
    vector< pair<jpeg_app_t, off_t> > blocks;
    fseek(fin, 2, SEEK_SET);
    bool done = false;
    
    while (!feof(fin) && !done) {
        uint16_t app_id = icc_tag::read_uint16(fin);
        uint16_t bsize  = icc_tag::read_uint16(fin); // this size excludes the size bytes
        
        if ((app_id & 0xff00) != 0xff00) {
            throw -1;
        }
        
        done = (app_id & 0xffe0) != 0xffe0;
        if (!done) {
            int app_n = (app_id & 0xff) - 0xe0;
            off_t fpos = ftell(fin);
            
            unsigned char sig[128];
            int sn = 0;
            while (sn < 127 && (sig[sn] = fgetc(fin)) != 0) sn++;
            
            if (strncasecmp((char*)sig, "ICC_PROFILE", 11) == 0) {
                blocks.push_back(make_pair(jpeg_app_t::ICC, ftell(fin) + 2)); // skip over chunk numbers
            }
            if (strncasecmp((char*)sig, "EXIF", 4) == 0) {
                int i;
                for (i=0; i < 4 && fgetc(fin) == 0; i++);
                blocks.push_back(make_pair(jpeg_app_t::EXIF, ftell(fin) - 1));
            }
            
            if (fseek(fin, fpos + bsize - 2, SEEK_SET) != 0) {
                throw -1;
            }
        }
    }
    return blocks;
}

