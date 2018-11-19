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
#include "img_frame.h"
#include "img_frame.moc"

#include "common.h"

#include <QKeyEvent>
#include <QWidget>
#include <stdio.h>

Img_frame::Img_frame(QWidget* parent) : parent(parent), ctrl_down(false) {
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

const int ScrollStep = 10;

void Img_frame::keyPressEvent(QKeyEvent* event) {
    switch (event->key()) {
    case Qt::Key_Control:
        ctrl_down = true;
        break;
    case Qt::Key_Equal:
        emit zoom_to_100();
        break;
    case Qt::Key_Plus:
        zoom(1);
        break;
    case Qt::Key_Minus:
        zoom(-1);
        break;
    default:
        QWidget::keyPressEvent(event);
    }
}

void Img_frame::keyReleaseEvent(QKeyEvent* event) {
    switch (event->key()) {
    case Qt::Key_Control:
        ctrl_down = false;
        break;
    default:
        QWidget::keyReleaseEvent(event);
    }
}

void Img_frame::wheelEvent(QWheelEvent* event) {
    if (ctrl_down) {
        if (event->delta() > 0) {
            emit zoom_in();
        } else {
            emit zoom_out();
        }
    }
}

void Img_frame::zoom(int dir) {
    if (dir < 0) {
        emit zoom_out();
    } else {
        emit zoom_in();
    }
}

