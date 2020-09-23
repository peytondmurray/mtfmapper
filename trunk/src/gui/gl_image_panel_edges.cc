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

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

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
    
    
    QOpenGLShader* edgefade_shader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char* edgefade_fsrc =
        "#version 120\n"
        "uniform vec4 dcolour;\n"
        "uniform float uform;\n"
        "varying mediump vec4 texc;\n"
        "void main(void)\n"
        "{\n"
        "    gl_FragColor = vec4(dcolour.xyz, 1 - uform*abs(texc.t - 0.5));\n"
        "}\n";
    if (!edgefade_shader->compileSourceCode(edgefade_fsrc)) {
        printf("edgefade_shader log: %s\n", edgefade_shader->log().toLocal8Bit().constData());
    }
    
    line_program = new QOpenGLShaderProgram;
    line_program->addShader(obj_vshader);
    line_program->addShader(edgefade_shader);
    line_program->bindAttributeLocation("vertex", prog_vert_att);
    line_program->bindAttributeLocation("texCoord", prog_texcoord_att);
    if (!line_program->link()) {
        printf("line_program log: %s\n", line_program->log().toLocal8Bit().constData());
    }
    
    QOpenGLShader* flat_shader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char* flat_fsrc =
        "#version 120\n"
        "uniform vec4 dcolour;\n"
        "void main(void)\n"
        "{\n"
        "    gl_FragColor = dcolour;\n"
        "}\n";
    if (!flat_shader->compileSourceCode(flat_fsrc)) {
        printf("flat_shader log: %s\n", flat_shader->log().toLocal8Bit().constData());
    }
    
    box_program = new QOpenGLShaderProgram;
    box_program->addShader(obj_vshader);
    box_program->addShader(flat_shader);
    box_program->bindAttributeLocation("vertex", prog_vert_att);
    box_program->bindAttributeLocation("texCoord", prog_texcoord_att);
    if (!box_program->link()) {
        printf("box_program log: %s\n", box_program->log().toLocal8Bit().constData());
    }
    
    QOpenGLShader *dots_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char *dots_fsrc =
        "#version 120\n"
        "uniform vec4 dcolour;\n"
        "uniform sampler2D texture;\n"
        "varying mediump vec4 texc;\n"
        "void main(void)\n"
        "{\n"
        "    vec2 uv = texc.st - vec2(0.5, 0.5);\n"
        "    float dist = sqrt(dot(uv, uv));\n"
        "    float t = smoothstep(0.20, 0.30, dist) - smoothstep(0.30, 0.40, dist);\n"
        "    gl_FragColor = vec4(dcolour.xyz, t);\n"
        "}\n";
    dots_fshader->compileSourceCode(dots_fsrc);
    
    dots_program = new QOpenGLShaderProgram;
    dots_program->addShader(obj_vshader);
    dots_program->addShader(dots_fshader);
    dots_program->bindAttributeLocation("vertex", prog_vert_att);
    dots_program->bindAttributeLocation("texCoord", prog_texcoord_att);
    dots_program->link();

    line_program->bind();
}

void GL_image_panel_edges::paint_overlay(void) {
    img_centre = QPointF(imgsize.width()/2, imgsize.height()/2);

    box_program->bind();
    box_program->setUniformValue("viewMatrix", view);
    box_program->setUniformValue("projectionMatrix", projection);
    
    dot_vbo.bind();
    box_program->enableAttributeArray(prog_vert_att);
    box_program->enableAttributeArray(prog_texcoord_att);
    box_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    box_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    for (auto& r : rois) {
        draw_box(r.get(0) - img_centre, r.get(1) - img_centre, 56, 0.2, 0.75, 0.75, 0.3);
    }
    
    line_program->bind();
    line_program->setUniformValue("viewMatrix", view);
    line_program->setUniformValue("projectionMatrix", projection);
    
    dot_vbo.bind();
    line_program->enableAttributeArray(prog_vert_att);
    line_program->enableAttributeArray(prog_texcoord_att);
    line_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    line_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    constexpr double colours[3][3] = {
        {0.2, 0.75, 0.75},
        {0.75, 0.35, 0.35},
        {0.965, 0.651, 0.149}
    };
    
    
    for (auto& r : rois) {
        int col_idx = &r == current_roi ? (closebox.is_valid() ? 1 : 0) : 0;
        draw_line(r.get(0) - img_centre, r.get(1) - img_centre, colours[col_idx][0], colours[col_idx][1], colours[col_idx][2]);
    }
    
    dots_program->bind();
    dots_program->setUniformValue("viewMatrix", view);
    dots_program->setUniformValue("projectionMatrix", projection);
    
    dot_vbo.bind();
    dots_program->enableAttributeArray(prog_vert_att);
    dots_program->enableAttributeArray(prog_texcoord_att);
    dots_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    dots_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    for (auto& r : rois) {
        int col_idx = &r == current_roi ? (closebox.is_valid() ? 1 : 2) : 0;
        draw_dot(r.get(0).x() - img_centre.x(), r.get(0).y() - img_centre.y(), colours[col_idx][0], colours[col_idx][1], colours[col_idx][2]);
        draw_dot(r.get(1).x() - img_centre.x(), r.get(1).y() - img_centre.y(), colours[col_idx][0], colours[col_idx][1], colours[col_idx][2]);
    }
    
    // finally, draw the close box if there is one
    if (closebox.is_valid()) {
        box_program->bind();
        box_program->setUniformValue("viewMatrix", view);
        box_program->setUniformValue("projectionMatrix", projection);
        
        dot_vbo.bind();
        box_program->enableAttributeArray(prog_vert_att);
        box_program->enableAttributeArray(prog_texcoord_att);
        box_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
        box_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
        
        double lsf = 1.0;
        if (scale_factor < 1) {
            lsf = 0.5/scale_factor + 0.75;
        }
        if (scale_factor == 1.0) {
            lsf = 1.25;
        }
        
        draw_box(
            closebox.get_pos() - lsf*10.0 * closebox.get_dir() - img_centre, 
            closebox.get_pos() + lsf*10.0 * closebox.get_dir() - img_centre, 
            lsf*20,
            0.75, 0.45, 0.45, 0.5
        );
        
        line_program->bind();
        line_program->setUniformValue("viewMatrix", view);
        line_program->setUniformValue("projectionMatrix", projection);
        
        dot_vbo.bind();
        line_program->enableAttributeArray(prog_vert_att);
        line_program->enableAttributeArray(prog_texcoord_att);
        line_program->setAttributeBuffer(prog_vert_att, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
        line_program->setAttributeBuffer(prog_texcoord_att, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
        draw_close_symbol(closebox.get_pos() - img_centre, closebox.get_dir(), 0.75, 0.2, 0.2);
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
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", r, g, b, 1.0);
    dots_program->setUniformValue("uform", float(0.0));
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
}

void GL_image_panel_edges::draw_line(QPointF start, QPointF end, double r, double g, double b, double width) {
    
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
    model.scale(line_len/20.0/lsf, std::min(width/scale_factor, 2.0)/10.0);
    
    line_program->setUniformValue("modelMatrix", model);
    line_program->setUniformValue("dcolour", r, g, b, 1.0);
    line_program->setUniformValue("uform", float(1.0));
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
}

void GL_image_panel_edges::draw_box(QPointF start, QPointF end, float box_width, double r, double g, double b, double t) {
    
    QPointF delta = end - start;
    double line_len = sqrt(QPointF::dotProduct(delta, delta));
  
    QMatrix4x4 model = QMatrix4x4();
    
    model = QMatrix4x4();
    model.translate(start.x() + delta.x()*0.5, start.y() + delta.y()*0.5, 0.1);
    
    model.rotate(atan2(delta.y(), delta.x())/M_PI*180.0, 0, 0, 1.0);
    
    model.scale(line_len/20.0, box_width/20.0);
    
    box_program->setUniformValue("modelMatrix", model);
    box_program->setUniformValue("dcolour", r, g, b, t);
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

void GL_image_panel_edges::draw_close_symbol(QPointF pos, QPointF dir, double r, double g, double b) {
    
    QPointF n(-dir.y(), dir.x());
    
    double lsf = 1.0;
    if (scale_factor < 1) {
        lsf = 0.5/scale_factor + 0.75;
    }
    if (scale_factor == 1.0) {
        lsf = 1.25;
    }
    
    QPointF dx = 8.0*dir*lsf;
    QPointF dy = 8.0*n*lsf;
    
    double lw = std::max(1.4, lsf);
    
    draw_line(pos - dy - dx, pos - dy + dx, r, g, b, lw);
    draw_line(pos + dy - dx, pos + dy + dx, r, g, b, lw);
    draw_line(pos - dy - dx, pos + dy - dx, r, g, b, lw);
    draw_line(pos - dy + dx, pos + dy + dx, r, g, b, lw);
    
    constexpr float rf = 5.0/8.0;
    draw_line(pos - rf*dy - rf*dx, pos + rf*dy + rf*dx, r, g, b, 1.25*lw);
    draw_line(pos - rf*dy + rf*dx, pos + rf*dy - rf*dx, r, g, b, 1.25*lw);
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

void GL_image_panel_edges::click_marker(QPoint pos, [[maybe_unused]] bool add) {
    QPointF icoords(pos.x(), pos.y());
    
    rois.push_back(GL_roi(icoords, icoords));
    current_roi = &rois.back();
    
    update();
    
}

void GL_image_panel_edges::line_endpoint(QPoint pos) {
    QPointF icoords(pos.x(), pos.y());
    
    if (current_roi && current_roi_handle_idx >= 0 && current_roi_handle_idx <= 1) {
        current_roi->get(current_roi_handle_idx) = icoords;
    }
    
    update();
}

void GL_image_panel_edges::clear_overlay(void) {
    rois.clear();
    current_roi = nullptr;
    closebox.make_invalid();
    state = NONE;
}

void GL_image_panel_edges::mousePressEvent(QMouseEvent* event) { 
    click = event->pos();
    
    bool any_handle = false;
    
    QPointF img_coords = locate(click);
    
    if ( (state == NONE || state == ROI_SELECTED) &&
        !(closebox.is_valid() && closebox.selected(img_coords)) ) {
        
        for (auto& r: rois) {
            int handle_idx = r.handle_selected(img_coords);
            
            if (handle_idx >= 0) {
                current_roi = &r;
                current_roi_handle_idx = handle_idx;
                closebox.make_invalid();
                state = DRAGGING_HANDLE;
                any_handle = true;
                
            }
        }
    }
    
    if (!any_handle) {
        event->ignore();
    }
    update();
}

static double sqr(double x) {
    return x*x;
}

void GL_image_panel_edges::check_roi_boxes_and_handles(QPointF img_coords) {
    bool any_roi_hit = false;
    
    for (auto& r: rois) {
        int box_hit = r.box_selected(img_coords, 10.0, 28.0);
        
        if (box_hit) {
            if (current_roi != &r) {
                current_roi = &r;
                state = ROI_SELECTED;
                closebox = GL_closebox(r.get(0), r.get(1), imgsize.width(), imgsize.height());
                closebox.make_valid();
                any_roi_hit = true;
            } else {
                // unselect if we hit the same ROI again
                current_roi = nullptr;
                state = NONE;
            }
        } else {
    
            int handle_idx = r.handle_selected(img_coords);
            
            if (handle_idx >= 0) {
                current_roi = &r;
                current_roi_handle_idx = handle_idx;
                state = MOVING;
            }
        }
    }
    if (!any_roi_hit) {
        closebox.make_invalid();
    }
}

void GL_image_panel_edges::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        double d = sqrt( sqr(event->pos().x() - click.x()) + sqr(event->pos().y() - click.y()) );
        QPoint img_coords = locate(click);
        QPoint current_img_coords = locate(event->pos());
        
        if (d <= 2.3) { // mouse "release" close enough to mouse "press" to consider it a click (rather than drag)
            
            switch (state) {
            case NONE:
                [[fallthrough]];
            case DRAGGING_HANDLE:
                state = AFTER_FIRST_CLICK;
                
                check_roi_boxes_and_handles(img_coords);
                
                if (state == AFTER_FIRST_CLICK) {
                    click_marker(img_coords, false);
                }
                break;
            case AFTER_FIRST_CLICK:
                // do nothing, remain in this state
                break;
            case MOVING:
                state = NONE;
                current_roi = nullptr;
                break;
            case ROI_SELECTED:
                state = NONE;
                
                // first check if we clicked the close box
                if (closebox.is_valid() && closebox.selected(img_coords)) {
                    for (auto it=rois.begin(); it != rois.end();) {
                        if (&(*it) == current_roi) {
                            it = rois.erase(it);
                            current_roi = nullptr;
                            closebox.make_invalid();
                        } else {
                            ++it;
                        }
                    }
                } else {
                    check_roi_boxes_and_handles(img_coords);
                }
                
                if (state == NONE) {
                    current_roi = nullptr;
                }
                break;
            }
            
            release_pt = img_coords;
            
            update();
            
        } else {
            if (state == DRAGGING_HANDLE) {
                line_endpoint(current_img_coords);
                current_roi = nullptr;
                state = NONE;
            }
        }
    }
    event->ignore();
    update();
}

void GL_image_panel_edges::mouseMoveEvent(QMouseEvent* event) {
    
    QPoint move_pt = locate(event->pos());
    
    switch (state) {
    case NONE:
        // silent move
        break;
    case AFTER_FIRST_CLICK:
        state = MOVING;
        break;
    case MOVING:
        state = MOVING;
        break;
    case DRAGGING_HANDLE:
        break;
    case ROI_SELECTED:
        // not possible to reach this code
        break;
    }
    
    if (state == MOVING || state == DRAGGING_HANDLE) {
        if (current_roi_handle_idx < 0) {
            current_roi_handle_idx = 1;
        }
        line_endpoint(move_pt);
    }
    
    // once we start dragging a handle, fudge the original
    // click coordinates to avoid "picking up" the handle
    // if you happen to end the drag manouvre close to 
    // the starting click
    if (state == DRAGGING_HANDLE) {
        click = QPoint(-(1 << 28), -(1 << 28));
    }
    
    update();
    event->ignore();
}

bool GL_image_panel_edges::save_rois(const QString& fname) {
    if (rois.size() == 0) {
        return false;
    }
    
    FILE* fout = fopen(fname.toLocal8Bit().constData(), "wt");
    if (!fout) {
        return false;
    }
    
    for (auto& r : rois) {
        fprintf(fout, "%.3lf %.3lf %.3lf %.3lf\n",
            r.get(0).x(), r.get(0).y(),
            r.get(1).x(), r.get(1).y()
        );
    }
    
    fclose(fout);
    
    return true;
}

