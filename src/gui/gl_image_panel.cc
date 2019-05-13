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
#include "gl_image_panel.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include <cmath>

GL_image_panel::GL_image_panel(QWidget *parent)
    : QOpenGLWidget(parent),
      program(0) {
      
    logger.debug("GL_image_panel ctor: OpenGL version: %d.%d, samples=%d\n", format().majorVersion(), format().minorVersion(), format().samples());
}

GL_image_panel::~GL_image_panel() {
    makeCurrent();
    vbo.destroy();
    for (int i=0; i < (int)textures.size(); i++) {
        delete textures[i];
    }
    delete program;
    doneCurrent();
}

QSize GL_image_panel::minimumSizeHint() const {
    return QSize(50, 50);
}

QSize GL_image_panel::sizeHint() const {
    return QSize(800, 800);
}

void GL_image_panel::initializeGL() {

    initializeOpenGLFunctions();
    logger.debug("initializeGL: OpenGL version: %d.%d, samples=%d\n", format().majorVersion(), format().minorVersion(), format().samples());
    
    load_image(*default_image);
    make_dots();

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_MULTISAMPLE);

    #define PROGRAM_VERTEX_ATTRIBUTE 0
    #define PROGRAM_TEXCOORD_ATTRIBUTE 1
    
    // shader program for main texture image
    QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    const char *vsrc =
        "#version 120\n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec4 texCoord;\n"
        "varying highp vec4 texc;\n"
        "uniform mediump mat4 projectionMatrix;\n"
        "uniform mediump mat4 modelMatrix;\n"
        "void main(void)\n"
        "{\n"
        "    gl_Position = projectionMatrix * modelMatrix * vertex;\n"
        "    texc = texCoord;\n"
        "}\n";
    vshader->compileSourceCode(vsrc);

    QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    const char *fsrc =
        "#version 120\n"
        "uniform sampler2D texture;\n"
        "varying highp vec4 texc;\n"
        "void main(void)\n"
        "{\n"
        "    gl_FragColor = texture2D(texture, texc.st);\n"
        "}\n";
    fshader->compileSourceCode(fsrc);

    program = new QOpenGLShaderProgram;
    program->addShader(vshader);
    program->addShader(fshader);
    program->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
    program->bindAttributeLocation("texCoord", PROGRAM_TEXCOORD_ATTRIBUTE);
    program->link();

    program->bind();
    program->setUniformValue("texture", 0);
    
    // shader program for click-dots
    QOpenGLShader *dots_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    const char *dots_vsrc =
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
    dots_vshader->compileSourceCode(dots_vsrc);
    
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
        "    float t = smoothstep(0.3, 0.5, dist);\n"
        "    gl_FragColor = vec4(dcolour.xyz, 1 - t);\n"
        "}\n";
    dots_fshader->compileSourceCode(dots_fsrc);
    
    dots_program = new QOpenGLShaderProgram;
    dots_program->addShader(dots_vshader);
    dots_program->addShader(dots_fshader);
    dots_program->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
    dots_program->bindAttributeLocation("texCoord", PROGRAM_TEXCOORD_ATTRIBUTE);
    dots_program->link();

    dots_program->bind();
    dots_program->setUniformValue("texture", 0);
}

void GL_image_panel::paintGL() {
    
    int w = size().width();
    int h = size().height();
    
    reset_scroll_range();
    
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT);
    
    view = QMatrix4x4();
    view.translate(-vp.centre.x + int(w/2), -vp.centre.y + int(h/2));
    view.scale(scale_factor);
    
    QMatrix4x4 projection;
    projection.ortho(0, w, h, 0, -1, 1);
    
    program->bind();
    program->setUniformValue("modelMatrix", view);
    program->setUniformValue("projectionMatrix", projection);
    
    vbo.bind();
    program->enableAttributeArray(PROGRAM_VERTEX_ATTRIBUTE);
    program->enableAttributeArray(PROGRAM_TEXCOORD_ATTRIBUTE);
    program->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    program->setAttributeBuffer(PROGRAM_TEXCOORD_ATTRIBUTE, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));

    
    for (int i = 0; i < (int)textures.size(); i++) {
        textures[i]->bind();
        glDrawArrays(GL_TRIANGLE_FAN, i * 4, 4);
    }
    
    QMatrix4x4 model;
    
    
    dots_program->bind();
    dots_program->setUniformValue("viewMatrix", view);
    dots_program->setUniformValue("projectionMatrix", projection);
    
    dots_vbo.bind();
    dots_program->enableAttributeArray(PROGRAM_VERTEX_ATTRIBUTE);
    dots_program->enableAttributeArray(PROGRAM_TEXCOORD_ATTRIBUTE);
    dots_program->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
    dots_program->setAttributeBuffer(PROGRAM_TEXCOORD_ATTRIBUTE, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
    
    constexpr double colours[3][3] = {
        {0.129, 0.627, 0.875},
        {0.604, 0.792, 0.329},
        {0.965, 0.651, 0.149}
    };
    
    if (dot_list[current_fname].size() > 0) {
        for (const auto& s: dot_list[current_fname]) {
            if (s.dot_no >= 1 && s.dot_no <= 3) {
                int j = s.dot_no - 1;
                draw_dot(s.p.x(), s.p.y(), colours[j][0], colours[j][1], colours[j][2]);    
            }
        }
    }
}

void GL_image_panel::draw_dot(double x, double y, double r, double g, double b) {
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
    model.scale(1.2*lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", 1, 1, 1, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
    model = QMatrix4x4();
    model.translate(x, y, 0.2);
    model.scale(1.05*lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", 0, 0, 0, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    
    
    model = QMatrix4x4();
    model.translate(x, y, 0.3);
    model.scale(lsf);
    dots_program->setUniformValue("modelMatrix", model);
    dots_program->setUniformValue("dcolour", r, g, b, 1.0);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}


void GL_image_panel::resizeGL(int /*width*/, int /*height*/) {
    
    double w = size().width();
    double h = size().height();
    if (imgsize.width() > w || imgsize.height() > h) {
        scale_factor = std::min(w / imgsize.width(), h / imgsize.height());
    }
    
    reset_scroll_range();
    
    update();
}

static int next_pow2(uint32_t x) {
    uint32_t k=1;
    while (k < 31 && (uint32_t(1) << k) < x) {
        k++;
    }
    return uint32_t(1) << k;
}

void GL_image_panel::load_image(const QString& fname) {
    if (fname.toStdString().compare(current_fname) == 0) {
        return;
    }
    
    cv::Mat cvimg;
    if (image_cache.find(fname.toStdString()) != image_cache.end()) {
        cvimg = image_cache[fname.toStdString()].fetch();
    } else { 
        cvimg = cv::imread(fname.toStdString(), cv::IMREAD_COLOR);
        
        // manual conversion, cvtColor seems a bit slow
        uint8_t* sptr = cvimg.data;
        uint8_t* sentinel = sptr + size_t(cvimg.rows)*size_t(cvimg.cols*3);
        while (sptr < sentinel) {
            std::swap(*sptr, *(sptr+2));
            sptr += 3;
        }
        
        image_cache[fname.toStdString()] = Cache_entry(cvimg);
        image_cache_size += image_cache[fname.toStdString()].size();
        
        trim_cache();
    }
    
    current_fname = fname.toStdString();
    
    load_image(cvimg);
}

void GL_image_panel::load_image(QImage& qimg) {

    cv::Mat cvimg(qimg.height(), qimg.width(), CV_8UC3, qimg.bits());
    current_fname = "mtf_mapper_logo";
    load_image(cvimg);
}

void GL_image_panel::load_image(cv::Mat cvimg) {
    bool keep_zoom = false;

    vbo.destroy();
    for (int i=0; i < (int)textures.size(); i++) {
        delete textures[i];
    }
    textures.clear();
    
    int hw_texw = 0;
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &hw_texw);
    
    int texw = std::min(hw_texw, 2048); // we do not want the last block to become too large
    
    if (cvimg.cols == imgsize.width() && cvimg.rows == imgsize.height()) {
        keep_zoom = true;
    }
    
    imgsize = QSize(cvimg.cols, cvimg.rows);
    
    // now trim down the texture block size if the image is small
    texw = std::min(texw, std::max(next_pow2(imgsize.width()), next_pow2(imgsize.height())));
    
    QVector<GLfloat> vertData;
    
    int cblocks = (int)ceil(imgsize.width()/double(texw));
    int rblocks = (int)ceil(imgsize.height()/double(texw));
    
    int x_off = imgsize.width() / 2;
    int y_off = imgsize.height() / 2;
    
    for (int r=0; r < rblocks; r++) {
        int rstart = r*texw;
        int rend = std::min((r+1)*texw, imgsize.height());
        for (int c=0; c < cblocks; c++) {
            int cstart = c*texw;
            int cend = std::min((c+1)*texw, imgsize.width());
            
            // obtain pow2 dimensions for this block
            int tex_width = next_pow2(cend-cstart);
            int tex_height = next_pow2(rend-rstart);
            
            vertData.append(cstart - x_off);
            vertData.append(rstart - y_off);
            vertData.append(0);
            vertData.append(0);
            vertData.append(0);
            
            vertData.append(cend - x_off);
            vertData.append(rstart - y_off);
            vertData.append(0);
            vertData.append(double(cend - cstart)/tex_width);
            vertData.append(0);
            
            vertData.append(cend - x_off);
            vertData.append(rend - y_off);
            vertData.append(0);
            vertData.append(double(cend - cstart)/tex_width);
            vertData.append(double(rend - rstart)/tex_height);
            
            vertData.append(cstart - x_off);
            vertData.append(rend - y_off);
            vertData.append(0);
            vertData.append(0);
            vertData.append(double(rend - rstart)/tex_height);
            
            // we copy parts of the cv::Mat image into a QImage
            // this way we can have large images (>10k dims were problematic QT 5.7), and
            // at the same time have auto-mipmaps (which we only get from QImage)
            if ((cend - cstart) == tex_width && (rend - rstart) == tex_height) {
                size_t offset = cstart * 3 + rstart * cvimg.elemSize()*cvimg.cols;
                QImage tex_block(cvimg.data + offset, tex_width, tex_height, cvimg.elemSize()*cvimg.cols, QImage::Format_RGB888);
                textures.push_back(new QOpenGLTexture(tex_block));
            } else {
                // only create a new image (power-of-2) if we have a partial block
                QImage tex_block(tex_width, tex_height, QImage::Format_RGB888);
                cv::Mat srcimg(cvimg, cv::Rect(cstart, rstart, cend-cstart, rend-rstart)); // the ROI we want
                cv::Mat teximg(tex_height, tex_width, CV_8UC3, tex_block.bits());
                teximg = cv::Scalar(255, 255, 255); // fill the image to avoid borders (during mipmap building)
                srcimg.copyTo(teximg(cv::Rect(0, 0, cend-cstart, rend-rstart)));
                textures.push_back(new QOpenGLTexture(tex_block));
            }
            
            // but this looks better (other than the artifacts, of course)
            textures.back()->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
            textures.back()->setMagnificationFilter(QOpenGLTexture::Linear);
            textures.back()->setWrapMode(QOpenGLTexture::ClampToEdge);
        }
    }
    
    double w = size().width();
    double h = size().height();
    if (!keep_zoom) {
        scale_factor = std::min(1.0, std::min(w / imgsize.width(), h / imgsize.height()));
        
        // reset centre position (in image)
        vp.centre = cv::Point2d(0,0);
    }
    
    reset_scroll_range();
    
    vbo.create();
    vbo.bind();
    vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
}

void GL_image_panel::make_dots(void) {
    QVector<GLfloat> vd;
    
    const double rad = 10;
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
    

    dots_vbo.create();
    dots_vbo.bind();
    dots_vbo.allocate(vd.constData(), vd.count() * sizeof(GLfloat));
}

void GL_image_panel::move(int nx, int ny) {
    vp.centre = cv::Point2d(nx, ny);
}

QPoint GL_image_panel::zoom(int step, int mx, int my) {
    double iw = imgsize.width();
    double ih = imgsize.height();
    int w = size().width();
    int h = size().height();
    
    // scale to next power of two, but skip to the next larger/smaller scale
    // if the relative change in scale is less than 25% when near the minimum scale
    double min_scale_factor = std::min(1.0, std::min(w/iw, h/ih));
    double next_scale_factor = scale_factor;
    if (step < 0) {
        double p = floor(log(scale_factor*0.5)/log(2.0));
        next_scale_factor = std::min(2.0, std::max(min_scale_factor, pow(2.0, p)));
        if (fabs(next_scale_factor - min_scale_factor)/min_scale_factor < 0.25) {
            next_scale_factor = min_scale_factor;
        }
    } else {
        if (step >= 0) {
            double p = floor(log(scale_factor*2)/log(2.0));
            next_scale_factor = std::max(min_scale_factor, std::min(2.0, pow(2.0, p)));
            if (fabs(next_scale_factor - min_scale_factor)/min_scale_factor < 0.25) {
                next_scale_factor = std::max(min_scale_factor, std::min(2.0, pow(2.0, p+1)));
            }
        }
    }
    scale_factor = next_scale_factor;
    
    QPointF m(mx, my);
    QPointF mp = view.inverted().map(m);
    
    QMatrix4x4 nview;
    nview.translate(-vp.centre.x + w/2, -vp.centre.y + h/2);
    nview.scale(scale_factor);
    
    QPointF mpi = nview.map(mp);
    QPointF delta = m - mpi;
    
    vp.centre.x -= delta.x();
    vp.centre.y -= delta.y();
    
    reset_scroll_range();
    
    return QPoint(vp.centre.x, vp.centre.y);
}

QPoint GL_image_panel::locate(QPoint pos) {
    double iw = imgsize.width();
    double ih = imgsize.height();
    
    QPointF p = view.inverted().map(QPointF(pos));
    double ix = p.x() + iw/2;
    double iy = p.y() + ih/2;
    
    dot_pos = QPoint(ix, iy);
    
    return dot_pos;
}


void GL_image_panel::click_marker(QPoint pos, bool add) {
    
    int ix = pos.x() - imgsize.width()/2;
    int iy = pos.y() - imgsize.height()/2;
    
    int dot_no = 0;
    for (auto it=dot_list.begin(); it != dot_list.end(); it++) {
        dot_no += it->second.size();
    }
    
    if (add) {
        if (dot_no < 3) {
            dot_list[current_fname].push_back(Sfr_marker(ix, iy, dot_no + 1));
        } else {
            // just rebuild the whole list, skipping the one we are replacing
            std::map<std::string, vector<Sfr_marker> > new_cache;
            for (const auto& f: dot_list) {
                for (const auto& s: f.second) {
                    if (s.dot_no != dot_no) {
                        new_cache[f.first].push_back(s);
                    }
                }
            }
            dot_list = new_cache;
            dot_list[current_fname].push_back(Sfr_marker(ix, iy, dot_no));
        }
    } else {
        dot_list.clear();
        dot_list[current_fname].push_back(Sfr_marker(ix, iy, 1));
    }
}

// A bit of a hack, but this is required to ensure that small images remain
// centered. Resetting this is necessary after scale changes, or window
// resizing, or loading a new image
void GL_image_panel::reset_scroll_range(void) {
    int w = size().width();
    int h = size().height();
    double iw = imgsize.width();
    double ih = imgsize.height();
    
    vp.set_x_range(-iw/2*scale_factor + w/2, iw/2*scale_factor - w/2);
    vp.set_y_range(-ih/2*scale_factor + h/2, ih/2*scale_factor - h/2);
}

void GL_image_panel::set_cache_size(uint64_t size) { 
    max_image_cache_size = size;
    trim_cache();
}

void GL_image_panel::trim_cache(void) {
    if (image_cache_size > max_image_cache_size) {
        // evict entries if necessary (we cannot evict pre-emptively because 
        // image dimensions are unknown until after imread())
        vector< pair<uint32_t, string> > cache_age;
        for (const auto& e: image_cache) {
            cache_age.push_back(make_pair(e.second.seq(), e.first));
        }
        sort(cache_age.begin(), cache_age.end());
        
        
        size_t ci = 0;
        while (image_cache_size > max_image_cache_size && ci < cache_age.size()) {
            auto cit = image_cache.find(cache_age[ci].second);
            image_cache_size -= cit->second.size();
            image_cache.erase(cit);
            ci++;
        }
    }
}

