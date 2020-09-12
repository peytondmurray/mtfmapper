/*
Copyright 2020 Frans van den Bergh. All rights reserved.

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
#ifndef GL_IMAGE_PANEL_EDGES_H
#define GL_IMAGE_PANEL_EDGES_H

#include "gl_image_panel.h"

class GL_image_panel_edges : public GL_image_panel {
  public:
    explicit GL_image_panel_edges(QWidget *parent = 0);
    
    virtual void click_marker(QPoint pos, bool add=false) override;
    virtual void clear_overlay(void) override;
    virtual void initialize_overlay(void) override;
    virtual void paint_overlay(void) override;
    
    void line_endpoint(QPoint pos);
    
  private:
    void make_dot_vbo();
    void draw_dot(double x, double y, double r, double g, double b);
    void draw_line(QPointF start, QPointF end, double r, double g, double b);
    void draw_roi(QPointF start, QPointF end, double r, double g, double b);
    void draw_box(QPointF start, QPointF end, float box_width, double r, double g, double b, double t = 1.0);
    
    QOpenGLShaderProgram* line_program;
    QOpenGLShaderProgram* box_program;
    QOpenGLShaderProgram* dots_program;
    
    QOpenGLShaderProgram* edges_program;
    QOpenGLBuffer edges_vbo;
    
    QOpenGLBuffer dot_vbo;
    vector<QPoint> dot_list;
    
    QPoint line_endp;
    
    //vector<Edge_marker> edge_list;
};

#endif

