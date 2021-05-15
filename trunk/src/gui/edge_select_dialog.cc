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
#include <QtWidgets> 
#include "edge_select_dialog.h"

Edge_select_dialog::Edge_select_dialog(QWidget* parent)
 : QDialog(parent, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint), parent(parent) {

    cancel_button = new QPushButton("Cancel");
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
    accept_button = new QPushButton("Accept");
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    help_button = new QPushButton("Help");
    help_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

    img_viewer = new GL_image_viewer(parent);
    img_viewer->set_clickable(true);
    img_viewer->set_resize_on_load(true);
    
    img_panel = new GL_image_panel_edges(img_viewer);
    img_panel->setMouseTracking(true);
    
    img_viewer->setViewport(img_panel);   // TODO: could combine these
    img_viewer->set_GL_widget(img_panel);
    
    icon_image = std::shared_ptr<QImage>(new QImage(QSize(256, 256), QImage::Format_RGB888));
    icon_image->fill(Qt::red);
    img_panel->set_default_image(icon_image.get());

    help_dialog = new Manual_roi_help_dialog(this);
    help_dialog->setModal(true);

    QGridLayout* vlayout = new QGridLayout(this);
    vlayout->addWidget(img_viewer, 0, 0, 1, 5);
    vlayout->addWidget(help_button, 1, 0);
    vlayout->addWidget(accept_button, 1, 3);
    vlayout->addWidget(cancel_button, 1, 4);
    
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(accept_button, SIGNAL(clicked()), this, SLOT( export_roi() ));
    connect(help_button, SIGNAL(clicked()), this, SLOT( show_help() ));
    
    setModal(true);
    setLayout(vlayout);
    setWindowTitle("Select one or more edge ROIs");
}

bool Edge_select_dialog::load_image(QString img_name) {
    bool load_success = img_viewer->load_image(img_name);
    
    if (!load_success) {
        return false;
    }
    
    img_panel->clear_overlay();
    update();
    // force a repaint on the GL panel once the event loop runs
    QTimer::singleShot(0, [this]() { 
        img_panel->update();
        QCoreApplication::processEvents(); 
    });
    return true;
}

void Edge_select_dialog::open(void) {
    show();
}

void Edge_select_dialog::export_roi(void) {
    bool roi_success = img_panel->save_rois(roi_file);
    if (roi_success) {
        accept();
    } else {
        close();
    }
}

void Edge_select_dialog::show_help(void) {
    help_dialog->setModal(false);
    help_dialog->show();
}

