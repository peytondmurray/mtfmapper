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
            QSize widgetSize(
                widget->zoom_scale()*widget->img_size().width(), 
                widget->zoom_scale()*widget->img_size().height()
            );
            
            verticalScrollBar()->setPageStep(areaSize.height());
            horizontalScrollBar()->setPageStep(areaSize.width());
            verticalScrollBar()->setRange(-0.5*(widgetSize.height()), 0.5*(widgetSize.height()) - areaSize.height());
            horizontalScrollBar()->setRange(-0.5*(widgetSize.width()), 0.5*(widgetSize.width()) - areaSize.width());
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
      {
          QPoint np = widget->zoom(e->angleDelta().y(), e->x(), e->y());
          
          // update scroll bars with new scale factor
          QSize areaSize = viewport()->size();
          QSize widgetSize(
              widget->zoom_scale()*widget->img_size().width(), 
              widget->zoom_scale()*widget->img_size().height()
          );
          verticalScrollBar()->setRange(-0.5*(widgetSize.height()), 0.5*(widgetSize.height()) - areaSize.height());
          horizontalScrollBar()->setRange(-0.5*(widgetSize.width()), 0.5*(widgetSize.width()) - areaSize.width());
          
          horizontalScrollBar()->setValue(np.x());
          verticalScrollBar()->setValue(np.y());
      }
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

void GL_image_viewer::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::RightButton) {
        panning = true;
        pan = event->pos();
        setCursor(Qt::ClosedHandCursor);
        event->accept();
        return;
    }
    if (event->button() == Qt::LeftButton) {
        panning = true;
        pan = event->pos();
        setCursor(Qt::ClosedHandCursor);
        click = event->pos();
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
        panning = false;
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
    event->ignore();
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