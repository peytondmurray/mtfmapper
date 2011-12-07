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

#include "exiv2_property.h"
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>

using std::string;
using std::ifstream;
using std::cout;
using std::endl;


Exiv2_property::Exiv2_property(QString bin_name, QString ifname, QString tfname)
: exiv2_binary(bin_name),ifname(ifname), tfname(tfname)
{

    // ping exiv2 to get camera make
    QString make = extract_property(QString("Exif.Image.Make"));

    if (make.startsWith(QString("NIKON"), Qt::CaseInsensitive)) mode = NIKON;
    else if (make.startsWith(QString("Canon"), Qt::CaseInsensitive)) mode = CANON;
    else if (make.length() == 0) mode = NONE;
    else mode = OTHER;

    p_af_tune = extract_af_tune();
    p_comment = extract_comment();
    p_focus_distance = extract_focus_distance();
    p_focal_length = extract_focal_length();
    p_aperture = extract_aperture();
}

char*   Exiv2_property::eat_non_whitespace(char* cp) {
    while (*cp != ' ' && *cp != '\t') {
        cp++;
    }
    return cp;
}

char*   Exiv2_property::eat_whitespace(char* cp) {
    while (*cp == ' ' || *cp == '\t') {
        cp++;
    }
    return cp;
}

QString Exiv2_property::extract_property(QString propname) {
    char* buffer = new char[4096];
    sprintf(buffer, "\"\"%s\" \"%s\" -g %s > %s\"",
		exiv2_binary.toLocal8Bit().constData(),
        ifname.toLocal8Bit().constData(),
        propname.toLocal8Bit().constData(),
        tfname.toLocal8Bit().constData()
    );

    int rval2 = system(buffer);

    if (rval2 < 0) {
        delete [] buffer;
        return QString("");
    }

    ifstream ifs(tfname.toLocal8Bit().constData());
    if (!ifs.fail()) {
        ifs.getline(buffer, 4096);
        buffer[4095] = 0; // ensure strlen will stop

        if (strlen(buffer) == 0) {
            return QString("N/A");
        }

        // seek over three whitespace regions
        char* cp = buffer;
        for (int i=0; i < 3; i++) {
            cp = eat_non_whitespace(cp);
            cp = eat_whitespace(cp);
        }

        QString rval = QString(cp);
        delete [] buffer;
        // now cp points to the last record, which is the one we want
        return rval;

    } else {
        printf("failed to open %s\n", tfname.toLocal8Bit().constData());
    }

    delete [] buffer;
    return QString("not found");
}


QString Exiv2_property::extract_af_tune(void) {
    if (mode != NONE) {
        switch (mode) {
        case NIKON:
            {
                QString aft = extract_property(QString("Exif.NikonAFT.AFFineTuneAdj"));
                bool ok;
                int value = aft.toInt(&ok);
                if (value > 20) {
                    value -= 256;
                }
                return QString("%1").arg(value);
            }
            
            break;
        case CANON:
            return extract_property(QString("Exif.NikonAFT.AFFineTuneAdj"));
            break;
        case OTHER:
        case NONE:
        default:
            return QString("N/A");
        }
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_comment(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.UserComment"));
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_focus_distance(void) {
    if (mode != NONE) {
        switch (mode) {
        case NIKON:
            return extract_property(QString("Exif.NikonLd3.FocusDistance"));
            break;
        case CANON:
            return extract_property(QString("Exif.CanonSi.SubjectDistance")); // right tag?
            break;
        case OTHER:
        case NONE:
        default:
            return QString("N/A");
        }
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_focal_length(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.FocalLength"));
    } else {
        return QString("N/A");
    }
}

QString Exiv2_property::extract_aperture(void) {
    if (mode != NONE) {
        return extract_property(QString("Exif.Photo.FNumber"));
    } else {
        return QString("N/A");
    }
}


