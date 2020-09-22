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
#ifndef EDGE_SELECT_DIALOG_H
#define EDGE_SELECT_DIALOG_H

#include <QDialog>

#include "gl_image_viewer.h"
#include "gl_image_panel_edges.h"


class QPushButton;
class QTextEdit;

class Edge_select_dialog : public QDialog {
  Q_OBJECT
  
  public:
    Edge_select_dialog(QWidget* parent);
    void load_image(QString img_name);
    GL_image_viewer* get_viewer(void) { return img_viewer; }
    GL_image_panel* get_panel(void) { return img_panel; }
    
    void set_size_hint(QSize size) { hinted_size = size; }
    QSize sizeHint(void) const override { return hinted_size; }
    
    void set_roi_file(const QString& fname) { roi_file = fname; }

  private:
    GL_image_viewer*  img_viewer;
    GL_image_panel_edges*   img_panel;
    QPushButton*      cancel_button;
    QPushButton*      accept_button;
    QWidget* parent = nullptr;
    QImage* icon_image;
    QSize hinted_size = QSize(0, 0);
    
    QString roi_file;
    
  public slots:
    void open();
    void export_roi();
    
};

#endif

