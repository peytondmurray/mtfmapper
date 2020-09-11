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
#include "gl_viewer_functor_edgeselect.h"
#include "gl_viewer_functor_annotated.h"

Edge_select_dialog::Edge_select_dialog(QWidget* parent)
 : QDialog(parent), parent(parent) {

    dismiss_button = new QPushButton("Dismiss");
    dismiss_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
    proceed_button = new QPushButton("Proceed");
    proceed_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    edge_functor = std::unique_ptr<GL_viewer_functor>(new GL_viewer_functor_edgeselect);
    
    img_viewer = new GL_image_viewer(parent, edge_functor.get());
    img_viewer->set_clickable(true);
    
    img_panel = new GL_image_panel_edges(img_viewer);
    img_panel->setMouseTracking(true);
    ((GL_viewer_functor_edgeselect*)edge_functor.get())->set_panel((GL_image_panel_edges*)img_panel); // TODO: there has to be a more elegant way ...
    
    img_viewer->setViewport(img_panel);   // TODO: could combine these
    img_viewer->set_GL_widget(img_panel);
    
    icon_image = new QImage(QSize(256, 256), QImage::Format_RGB888);
    icon_image->fill(Qt::red);
    img_panel->set_default_image(icon_image);
    
    QGridLayout* vlayout = new QGridLayout(this);
    vlayout->addWidget(img_viewer, 0, 0, 1, 3);
    vlayout->addWidget(proceed_button, 1, 0);
    vlayout->addWidget(dismiss_button, 1, 1);
    
    connect(dismiss_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(proceed_button, SIGNAL(clicked()), this, SLOT( accept() ));
    
    setLayout(vlayout);
    setWindowTitle("Select an edge ROI");
}

void Edge_select_dialog::load_image(QString img_name) {

    img_panel->load_image(img_name);
    img_panel->clear_overlay();
    update();
    
}

void Edge_select_dialog::open(void) {
    show();
}


