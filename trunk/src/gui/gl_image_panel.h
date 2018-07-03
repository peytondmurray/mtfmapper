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
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <vector>
using std::vector;
#include <map>
using std::pair;
using std::make_pair;
#include <string>
using std::string;

#include <cmath>
#include <algorithm>

#include <opencv2/opencv.hpp>

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>

#include "sfr_marker.h"
#include "cache_entry.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram);
QT_FORWARD_DECLARE_CLASS(QOpenGLTexture)

class GL_image_panel : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit GL_image_panel(QWidget *parent = 0);
    ~GL_image_panel();

    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;
    QSize img_size() const { return imgsize; }
    
    void move(int nx, int ny);
    QPoint zoom(int step, int x, int y);
    
    double zoom_scale(void) const { return scale_factor; }
    QPoint pos(void) const { return QPoint(tlx, tly); }
    QPoint locate(QPoint pos);
    void click_marker(QPoint pos, bool add=false);
    void load_image(const QString& fname);
    void load_image(QImage& qimg);
    void set_cache_size(uint64_t size);
    void clear_dots(void) { dot_list.clear(); }
    void set_default_image(QImage* qimg) { default_image = qimg; }

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int width, int height) override;

private:
    void make_dots();
    void draw_dot(double x, double y, double r, double g, double b);
    void load_image(cv::Mat cvimg);
    void reset_bias(void);
    void trim_cache(void);

    QColor clearColor;
    QPoint lastPos;
    vector<QOpenGLTexture*> textures;
    QOpenGLShaderProgram* program;
    QOpenGLBuffer vbo;
    
    QOpenGLShaderProgram* dots_program;
    QOpenGLBuffer dots_vbo;
    
    std::map<std::string, Cache_entry> image_cache;
    uint64_t image_cache_size = 0;
    uint64_t max_image_cache_size = uint64_t(1) << 30;
    std::map<std::string, vector<Sfr_marker> > dot_list;
    
    QSize imgsize = QSize(0,0);
    
    double tlx;
    double tly;
    double scale_factor = 1.0;
    double bias_x = 0;
    double bias_y = 0;
    
    QImage* default_image;
    QPoint dot_pos;
    std::string current_fname;
};

#endif
