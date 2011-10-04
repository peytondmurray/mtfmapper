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
#include "settings_dialog.h"
#include "settings_dialog.moc"

#include <iostream>
using std::cout;
using std::endl;

const QString setting_threshold = "setting_threshold";
const QString setting_threshold_default = "0.75";
const QString setting_linear_gamma = "setting_gamma";
const Qt::CheckState setting_linear_gamma_default = Qt::Unchecked;
const QString setting_annotation = "setting_annotation";
const Qt::CheckState setting_annotation_default = Qt::Checked;
const QString setting_profile = "setting_profile";
const Qt::CheckState setting_profile_default = Qt::Checked;
const QString setting_grid = "setting_grid";
const Qt::CheckState setting_grid_default = Qt::Checked;

Settings_dialog::Settings_dialog(QWidget *parent)
 : settings("mtfmapper", "mtfmapper")
{
    arguments_label = new QLabel(tr("Arguments:"));
    arguments_line  = new QLineEdit;
    threshold_label = new QLabel(tr("Threshold:"));
    threshold_line  = new QLineEdit;
    
    accept_button = new QPushButton(tr("&Accept"));
    cancel_button = new QPushButton(tr("&Cancel"));
    
    cb_linear_gamma = new QCheckBox("Linear gamma (8 bit)");
    cb_annotation   = new QCheckBox("Annotation");
    cb_profile      = new QCheckBox("Profile");
    cb_grid         = new QCheckBox("Grid");
    
    threshold_line->setText(settings.value(setting_threshold, setting_threshold_default).toString());
    cb_linear_gamma->setCheckState(
        (Qt::CheckState)settings.value(setting_linear_gamma, setting_linear_gamma_default).toInt()
    );
    cb_annotation->setCheckState(
        (Qt::CheckState)settings.value(setting_annotation, setting_annotation_default).toInt()
    );
    cb_profile->setCheckState(
        (Qt::CheckState)settings.value(setting_profile, setting_profile_default).toInt()
    );
    cb_grid->setCheckState(
        (Qt::CheckState)settings.value(setting_grid, setting_grid_default).toInt()
    );
    
    
    QGroupBox* v2GroupBox = new QGroupBox(tr("Flags"));
    QVBoxLayout *cb_layout = new QVBoxLayout;
    cb_layout->addWidget(cb_linear_gamma);
    cb_layout->addWidget(cb_annotation);
    cb_layout->addWidget(cb_profile);
    cb_layout->addWidget(cb_grid);
    v2GroupBox->setLayout(cb_layout);
    
    QGroupBox* vGroupBox = new QGroupBox(tr("Settings"));
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(v2GroupBox, 0, 0, 1, 2);
    vlayout->addWidget(threshold_label, 1, 0);
    vlayout->addWidget(threshold_line, 1, 1);
    vlayout->addWidget(arguments_label, 2, 0);
    vlayout->addWidget(arguments_line, 2, 1);
    vlayout->addWidget(accept_button, 3, 0);
    vlayout->addWidget(cancel_button, 3, 1);
    vGroupBox->setLayout(vlayout);
    
    connect(accept_button, SIGNAL(clicked()), this, SLOT( save_and_close() ));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    
    setLayout(vlayout);
}

void Settings_dialog::open() {

    show();
}

void Settings_dialog::send_argument_string(void) {
    QString args = QString("-t %1").arg(threshold_line->text());
    if (cb_linear_gamma->checkState()) {
        args = args + QString(" -l");
    }
    
    if (!cb_annotation->checkState()) {
        args = args + QString(" -a");
    }
    
    if (!cb_profile->checkState()) {
        args = args + QString(" -p");
    }
    
    if (!cb_grid->checkState()) {
        args = args + QString(" -s");
    }
    
    emit argument_string(args);
}

void Settings_dialog::save_and_close() {
    printf("saving settings\n");
    settings.setValue(setting_threshold, threshold_line->text());
    settings.setValue(setting_linear_gamma, cb_linear_gamma->checkState());
    settings.setValue(setting_annotation, cb_annotation->checkState());
    settings.setValue(setting_profile, cb_profile->checkState());
    settings.setValue(setting_grid, cb_grid->checkState());
    
    send_argument_string();
    
    close();
}

