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
#include "gl_image_viewer.h"
#include "mtfmapper_app.h"
#include <QtCore/QtCore>
#include <QtWidgets/QtWidgets>

GL_image_viewer::GL_image_viewer(mtfmapper_app* parent)
: QAbstractScrollArea(parent), parent_app(parent) {

    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    
    verticalScrollBar()->setRange(0, 1000);
    horizontalScrollBar()->setRange(0, 1000);
}



bool GL_image_viewer::viewportEvent(QEvent* e)  {
  
    if (e->type() == QEvent::Paint) {
    
        // TODO: must find a better place (time?) to update scrollbars
        // can only do this once the image has been loaded, though
        if (must_update_bars) {
            QSize areaSize = viewport()->size();
            
            verticalScrollBar()->setPageStep(areaSize.height());
            horizontalScrollBar()->setPageStep(areaSize.width());
            
            Image_viewport vp = widget->get_viewport();
            verticalScrollBar()->setRange(vp.y_range.x, vp.y_range.y);
            horizontalScrollBar()->setRange(vp.x_range.x, vp.x_range.y);
            
            must_update_bars = false;
        }
        return false;
    }
    
    if (e->type() == QEvent::Resize) {
        must_update_bars = true;
        return false;
    }
    
    return QAbstractScrollArea::viewportEvent(e);
}

void GL_image_viewer::scrollContentsBy(int /*dx*/, int /*dy*/) {
    int hvalue = horizontalScrollBar()->value();
    int vvalue = verticalScrollBar()->value();
    
    widget->move(hvalue, vvalue);
    viewport()->update();
}

void GL_image_viewer::wheelEvent(QWheelEvent* e) {
    const int scroll_size = 1;
    switch (e->modifiers()) {
    case Qt::ControlModifier:
      // directly ask viewport to zoom ...
      zoom_action(e->angleDelta().y(), e->x(), e->y());
      break;  
    case Qt::ShiftModifier:
      //scroll x
      horizontalScrollBar()->setValue(horizontalScrollBar()->value() - scroll_size*e->angleDelta().y());
      break;
    case Qt::NoModifier:
    default:
      // scroll y
      verticalScrollBar()->setValue(verticalScrollBar()->value() - scroll_size*e->angleDelta().y());
      break;
    }
    
    scrollContentsBy(0, 0);
}

void GL_image_viewer::set_GL_widget(GL_image_panel* w) { 
    widget = w;
}

void GL_image_viewer::keyPressEvent(QKeyEvent* event) {
    
    if (last_mouse_pos == QPoint(0, 0)) { // just in case this is possible
        last_mouse_pos = QPoint(width()/2, height()/2);
    }
    
    if (event->key() == Qt::Key_Plus) {
        zoom_action(1.0, last_mouse_pos.x(), last_mouse_pos.y());
    }
    if (event->key() == Qt::Key_Minus) {
        zoom_action(-1.0, last_mouse_pos.x(), last_mouse_pos.y());
    }
    
}

void GL_image_viewer::mousePressEvent(QMouseEvent* event) {
    last_mouse_pos = event->pos();
    if (event->button() == Qt::LeftButton) {
        panning = true;
        pan = event->pos();
        setCursor(Qt::ClosedHandCursor);
        click = event->pos();
        event->accept();
        return;
    }
    if (event->button() == Qt::RightButton) {
        zooming = true;
        zoom_pos = event->pos();
        zoom_pos_temp = zoom_pos;
        setCursor(Qt::SizeVerCursor);
        event->accept();
        return;
    }
    event->ignore();
}

static double sqr(double x) {
    return x*x;
}

void GL_image_viewer::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::RightButton) {
        zooming = false;
        setCursor(Qt::ArrowCursor);
        event->accept();
        return;
    }
    if (event->button() == Qt::LeftButton) {
        double d = sqrt( sqr(event->pos().x() - click.x()) + sqr(event->pos().y() - click.y()) );
        if (is_clickable && d < 10) { // mouse "release" close enough to mouse "press" to consider it a click (rather than drag)
            
            QPoint img_coords = widget->locate(click);
            bool valid = parent_app->edge_selected(img_coords.x(), img_coords.y(), event->modifiers().testFlag(Qt::ControlModifier), event->modifiers().testFlag(Qt::ShiftModifier));
            
            if (valid) {
                widget->click_marker(img_coords, event->modifiers() == Qt::ShiftModifier);
                widget->update();
            }
        } 
        panning = false;
        setCursor(Qt::ArrowCursor);
        event->accept();
        return;
    }
    event->ignore();
}

void GL_image_viewer::mouseMoveEvent(QMouseEvent* event) {
    if (panning) {
        horizontalScrollBar()->setValue(horizontalScrollBar()->value() - (event->x() - pan.x()));
        verticalScrollBar()->setValue(verticalScrollBar()->value() - (event->y() - pan.y()));
        pan = event->pos();
        event->accept();
        return;
    }
    if (zooming) {
        double dist = zoom_pos_temp.y() - event->y();
        
        if (fabs(dist) > 20) {
            zoom_pos_temp = event->pos();
            zoom_action(dist, zoom_pos.x(), zoom_pos.y());
        }
    }
    event->ignore();
}

void GL_image_viewer::zoom_action(double direction, int zx, int zy) {
    QPoint np = widget->zoom(direction, zx, zy);
    
    QSize areaSize = viewport()->size();  
    verticalScrollBar()->setPageStep(areaSize.height());
    horizontalScrollBar()->setPageStep(areaSize.width());
    
    Image_viewport vp = widget->get_viewport();
    verticalScrollBar()->setRange(vp.y_range.x, vp.y_range.y);
    horizontalScrollBar()->setRange(vp.x_range.x, vp.x_range.y);
    
    horizontalScrollBar()->setValue(np.x());
    verticalScrollBar()->setValue(np.y());
}

void GL_image_viewer::load_image(const QString& fname) {
    widget->load_image(fname);
    must_update_bars = true;
    widget->update();
}

void GL_image_viewer::load_image(QImage* qimg) {
    widget->load_image(*qimg);
    must_update_bars = true;
    widget->update();
}

void GL_image_viewer::set_clickable(bool b) {
    is_clickable = b;
}

void GL_image_viewer::clear_dots(void) {
    widget->clear_dots();
    widget->update();
}