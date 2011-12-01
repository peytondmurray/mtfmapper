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
#include <QtGui> 
#include "help_dialog.h"
#include "help_dialog.moc"

#include "common.h"
#include "config.h"

#define HELP_DIR_FORMAT(x) "##x"

Help_dialog::Help_dialog(QWidget* parent ATTRIBUTE_UNUSED) {

    dismiss_button = new QPushButton("Dismiss");
    dismiss_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    body = new QTextBrowser(parent);
    #if WIN32_
    if (mtfmapper_HTML_HELP_DIR.compare("./") == 0) {
        body_text = QCoreApplication::applicationDirPath() + QString("\\doc\\html\\mtf_mapper.html");
    } else {
        body_text = QString(mtfmapper_HTML_HELP_DIR.c_str());
    }
    #else
    body_text = QString(mtfmapper_HTML_HELP_DIR.c_str());
    #endif
    body->setSource(QUrl(body_text));
    body->setMinimumWidth(800);
    body->setMinimumHeight(600);
    body->setReadOnly(true);

    QGroupBox* vGroupBox = new QGroupBox(tr("About"));
    QGridLayout* vlayout = new QGridLayout;

    vlayout->addWidget(body, 0, 0, 1, 3);
    vlayout->addWidget(dismiss_button, 1, 1);
    vGroupBox->setLayout(vlayout);
    
    connect(dismiss_button, SIGNAL(clicked()), this, SLOT( close() ));
    
    setLayout(vlayout);
    setCentralWidget(body);
}


void Help_dialog::open(void) {
    body->setSource(QUrl(body_text));
    body->show();
    show();
}


