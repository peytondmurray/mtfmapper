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
#include "imgviewer.h"
#include "mtfmapper_app.h"


Imgviewer::Imgviewer(QGraphicsScene* scene, mtfmapper_app* zoom_parent, QWidget* parent)
 : QGraphicsView(scene, parent), scene(scene), zoom_parent(zoom_parent), clickable(false) {
 
}
    
void Imgviewer::wheelEvent(QWheelEvent* event) {
    if (event->modifiers().testFlag(Qt::ControlModifier))  {
        if (event->delta() > 0) {
            zoom_parent->zoom_in();
        } else {
            zoom_parent->zoom_out();
        }
    } else {
        QGraphicsView::wheelEvent(event);
    }
}

void Imgviewer::enterEvent(QEvent* event) {
    QGraphicsView::enterEvent(event);
    if (clickable) {
        viewport()->setCursor(Qt::ArrowCursor);
    }
}

void Imgviewer::mousePressEvent(QMouseEvent* event) {
    click_down_pos = event->pos();
    QGraphicsView::mousePressEvent(event); // pass along the event in case someone else needs it
    viewport()->setCursor(Qt::ClosedHandCursor);
}

static inline double sqr(double x) {
    return x*x;
}

void Imgviewer::mouseReleaseEvent(QMouseEvent* event) {
    double d = sqrt( sqr(event->pos().x() - click_down_pos.x()) + sqr(event->pos().y() - click_down_pos.y()) );
    if (d < 10) { // mouse "release" close enough to mouse "press" to consider it a click (rather than drag)
        QPointF mpos = mapToScene(event->pos().x(), event->pos().y());
        zoom_parent->edge_selected(mpos.x(), mpos.y(), event->modifiers().testFlag(Qt::ControlModifier), event->modifiers().testFlag(Qt::ShiftModifier));
    }
    QGraphicsView::mouseReleaseEvent(event); // pass along the event in case someone else needs it
    if (clickable) {
        viewport()->setCursor(Qt::ArrowCursor);
    } else {
        viewport()->setCursor(Qt::OpenHandCursor);
    }
}
