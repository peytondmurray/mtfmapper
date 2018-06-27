#ifndef SFR_MARKER_H
#define SFR_MARKER_H

class Sfr_marker {
  public:
    Sfr_marker(const QPoint& p, int dot_no) : p(p), dot_no(dot_no) {};
    Sfr_marker(int ix=0, int iy=0, int dot_no=0) : p(ix, iy), dot_no(dot_no) {};
    
    QPoint p;
    int dot_no;
};

#endif
