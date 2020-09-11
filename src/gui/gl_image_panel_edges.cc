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
#include "include/logger.h"
#include "gl_image_panel_edges.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include <cmath>

GL_image_panel_edges::GL_image_panel_edges(QWidget *parent)
    : GL_image_panel(parent) {
      
    logger.debug("GL_image_panel_edges ctor: OpenGL version: %d.%d, samples=%d\n", format().majorVersion(), format().minorVersion(), format().samples());
    set_cache_state(false);
}


void GL_image_panel_edges::initialize_overlay(void) {
    make_dot_vbo();
    
    // shader program for general transformed objects
    QOpenGLShader* obj_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    const char* obj_vsrc =
        "#version 120\n"
        "attribute highp vec4 vertex;\n"
        "attribute mediump vec4 texCoord;\n"
        "varying mediump vec4 texc;\n"
        "uniform mediump mat4 projectionMatrix;\n"
        "uniform mediump mat4 modelMatrix;\n"
        "uniform mediump mat4 viewMatrix;\n"
        "void main(void)\n"
        "{\n"
        "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * vertex;\n"
        "    texc = texCoord;\n"
        "}\n";
    obj_vshader->compileSourceCode(obj_vsrc);
    
    
    QOpenGLShader* uniform_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char* uniform_fsrc =
        "#version 120\n"
        "uniform vec4 dcolour;\n"
        "uniform float uform;\n"
        "varying mediump vec4 texc;\n"
        "void main(void)\n"
        "{\n"
        "    gl_FragColor = vec4(dcolour.xyz, 1 - uform*abs(texc.t - 0.5));\n"
        "}\n";
    if (!uniform_fshader->compileSourceCode(uniform_fsrc)) {
        printf("uniform_fshader log: %s\n", uniform_fshader->log().toLocal8Bit().constData());
    }
    
    obj_program = new QOpenGLShaderProgram;
    obj_program->addShader(obj_vshader);
    obj_program->addShader(uniform_fshader);
    obj_program->bindAttributeLocation("vertex", prog_vert_att);
    obj_program->bindAttributeLocation("texCoord", prog_texcoord_att);
    if (!obj_program->link()) {
        printf("obj_program log: %s\n", obj_program->log().toLocal8Bit().constData());
    }

    obj_program->bind();
}

void GL_image_panel_edges::paint_overlay(void) {
    obj_program->bind();
    obj_program->setUniformValue("viewMatrix", view);
    obj_program->setUniformValue("projectionMatrix", projection);
    
    dot_vbo.bind();
    obj_program->enableAttributeArray(prog_vert_att);
    obj_program->enableAttributeArray(prog_texcoord_att);
    obj_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    obj_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    constexpr double colours[3][3] = {
        {0.129, 0.627, 0.875},
        {0.604, 0.792, 0.329},
        {0.965, 0.651, 0.149}
    };
    
    for (auto p: dot_list) {
        draw_dot(p.x(), p.y(), colours[0][0], colours[0][1], colours[0][2]);
    }
    
    if (dot_list.size() > 0) {
        draw_dot(line_endp.x(), line_endp.y(), colours[1][0], colours[1][1], colours[1][2]);
        
        //draw_line(dot_list.back(), line_endp, colours[2][0], colours[2][1], colours[2][2]);
        draw_roi(dot_list.back(), line_endp, colours[2][0], colours[2][1], colours[2][2]);
    } 
}


void GL_image_panel_edges::draw_dot(double x, double y, double r, double g, double b) {
    
    double lsf = 1.0;
    if (scale_factor < 1) {
        lsf = 0.5/scale_factor + 0.75;
    }
    if (scale_factor == 1.0) {
        lsf = 1.25;
    }
  
    QMatrix4x4 model = QMatrix4x4();
    model = QMatrix4x4();
    model.translate(x, y, 0.1);
    model.scale(lsf);
    obj_program->setUniformValue("modelMatrix", model);
    obj_program->setUniformValue("dcolour", r, g, b, 1.0);
    obj_program->setUniformValue("uform", float(0.0));
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
}

void GL_image_panel_edges::draw_line(QPointF start, QPointF end, double r, double g, double b) {
    
    double lsf = 1.0;
    if (scale_factor < 1) {
        lsf = 0.5/scale_factor + 0.75;
    }
    if (scale_factor == 1.0) {
        lsf = 1.25;
    }
    
    QPointF delta = end - start;
    double line_len = sqrt(QPointF::dotProduct(delta, delta));
  
    QMatrix4x4 model = QMatrix4x4();
    
    model = QMatrix4x4();
    model.translate(start.x() + delta.x()*0.5, start.y() + delta.y()*0.5, 0.1);
    model.scale(lsf);
    
    model.rotate(atan2(delta.y(), delta.x())/M_PI*180.0, 0, 0, 1.0);
    
    // TODO: not quite 100% happy with line thickness scaling with image size
    model.scale(line_len/20.0/lsf, std::min(1.0/scale_factor, 2.0)/10.0);
    
    obj_program->setUniformValue("modelMatrix", model);
    obj_program->setUniformValue("dcolour", r, g, b, 1.0);
    obj_program->setUniformValue("uform", float(1.0));
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
}

void GL_image_panel_edges::draw_roi(QPointF start, QPointF end, double r, double g, double b) {
    
    draw_line(start, end, r, g, b);
    
    QPointF delta = end - start;
    QPointF n(-delta.y(), delta.x());
    n /= sqrt(n.x()*n.x() + n.y()*n.y());
    
    constexpr float roi_hw = 16.0;
    
    draw_line(start - roi_hw*n, end - roi_hw*n, r, g, b);
    draw_line(start + roi_hw*n, end + roi_hw*n, r, g, b);
    
    draw_line(start - roi_hw*n, start + roi_hw*n, r, g, b);
    draw_line(end - roi_hw*n, end + roi_hw*n, r, g, b);
}

void GL_image_panel_edges::make_dot_vbo(void) {
    QVector<GLfloat> vd;
    
    constexpr double rad = 10;
    vd.append(-rad);
    vd.append(-rad);
    vd.append(0.1);
    vd.append(0);
    vd.append(0);
    
    vd.append(rad);
    vd.append(-rad);
    vd.append(0.1);
    vd.append(1);
    vd.append(0);
    
    vd.append(rad);
    vd.append(rad);
    vd.append(0.1);
    vd.append(1);
    vd.append(1);
    
    vd.append(-rad);
    vd.append(rad);
    vd.append(0.1);
    vd.append(0);
    vd.append(1);
    

    dot_vbo.create();
    dot_vbo.bind();
    dot_vbo.allocate(vd.constData(), vd.count() * sizeof(GLfloat));
}

void GL_image_panel_edges::click_marker(QPoint pos, bool add) {
    int ix = pos.x() - imgsize.width()/2;
    int iy = pos.y() - imgsize.height()/2;
    
    printf("adding click marker (%d, %d)  / (%d, %d)\n", pos.x(), pos.y(), ix, iy);
    
    dot_list.push_back(QPoint(ix, iy));
    line_endp = QPoint(ix, iy);
    update();
    
}

void GL_image_panel_edges::line_endpoint(QPoint pos) {
    int ix = pos.x() - imgsize.width()/2;
    int iy = pos.y() - imgsize.height()/2;
    
    line_endp = QPoint(ix, iy);
    update();
}

void GL_image_panel_edges::clear_overlay(void) {
    dot_list.clear();
}
