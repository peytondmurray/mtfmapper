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
const QString setting_gnuplot = "setting_gnuplot";
const QString setting_exiv = "setting_exiv";
#ifdef _WIN32
static QString setting_gnuplot_default = "gnuplot.exe";
static QString setting_exiv_default = "exiv2.exe";
#else
const QString setting_gnuplot_default = "/usr/bin/gnuplot";
const QString setting_exiv_default = "/usr/bin/exiv2";
#endif

Settings_dialog::Settings_dialog(QWidget *parent)
 : settings("mtfmapper", "mtfmapper")
{
    arguments_label = new QLabel(tr("Arguments:"));
    arguments_line  = new QLineEdit;
    threshold_label = new QLabel(tr("Threshold:"));
    threshold_line  = new QLineEdit;

    gnuplot_label  = new QLabel(tr("gnuplot executable:"));
    gnuplot_line   = new QLineEdit;
    //gnuplot_line->resize(200,20);
    gnuplot_button = new QPushButton(tr("Browse"));

    exiv_label  = new QLabel(tr("exiv2 executable:"));
    exiv_line   = new QLineEdit;
    //exiv_line->resize(200,20);
    exiv_button = new QPushButton(tr("Browse"));

    accept_button = new QPushButton(tr("&Accept"));
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    cancel_button = new QPushButton(tr("&Cancel"));
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
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

    #ifdef _WIN32
    setting_gnuplot_default = QCoreApplication::applicationDirPath() + QString("\\gnuplot\\gnuplot.exe");
    setting_exiv_default = QCoreApplication::applicationDirPath() + QString("\\exiv2\\exiv2.exe");
    #endif

    gnuplot_line->setText(settings.value(setting_gnuplot, setting_gnuplot_default).toString());
    exiv_line->setText(settings.value(setting_exiv, setting_exiv_default).toString());
    
    QGroupBox* v2GroupBox = new QGroupBox(tr("Flags"));
    QVBoxLayout *cb_layout = new QVBoxLayout;
    cb_layout->addWidget(cb_linear_gamma);
    cb_layout->addWidget(cb_annotation);
    cb_layout->addWidget(cb_profile);
    cb_layout->addWidget(cb_grid);
    v2GroupBox->setLayout(cb_layout);

    QGroupBox* v3GroupBox = new QGroupBox(tr("Helpers"));
    QGridLayout *helper_layout = new QGridLayout;
    helper_layout->addWidget(gnuplot_label, 0, 0);
    helper_layout->addWidget(gnuplot_line, 1, 0);
    helper_layout->addWidget(gnuplot_button, 1, 1);
    helper_layout->addWidget(exiv_label, 2, 0);
    helper_layout->addWidget(exiv_line, 3, 0);
    helper_layout->addWidget(exiv_button, 3, 1);
    v3GroupBox->setLayout(helper_layout);

    QGroupBox* advanced = new QGroupBox(tr("Advanced"));
    QGridLayout* adv_layout = new QGridLayout;
    adv_layout->addWidget(threshold_label, 0, 0);
    adv_layout->addWidget(threshold_line, 0, 1);
    adv_layout->addWidget(arguments_label, 1, 0);
    adv_layout->addWidget(arguments_line, 1, 1);
    advanced->setLayout(adv_layout);

    
    QGroupBox* vGroupBox = new QGroupBox(tr("Settings"));
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(v2GroupBox, 0, 0, 1, 2);
    vlayout->addWidget(v3GroupBox, 1, 0, 1, 2);
    vlayout->addWidget(advanced, 2, 0, 1, 2);
    vlayout->addWidget(accept_button, 3, 0);
    vlayout->addWidget(cancel_button, 3, 1);
    vGroupBox->setLayout(vlayout);
    
    connect(accept_button, SIGNAL(clicked()), this, SLOT( save_and_close() ));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(gnuplot_button, SIGNAL(clicked()), this, SLOT( browse_for_gnuplot() ));
    connect(exiv_button, SIGNAL(clicked()), this, SLOT( browse_for_exiv() ));
    
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
    
    if (cb_annotation->checkState()) {
        args = args + QString(" -a");
    }
    
    if (cb_profile->checkState()) {
        args = args + QString(" -p");
    }
    
    if (cb_grid->checkState()) {
        args = args + QString(" -s");
    }
    
    emit argument_string(args);
}

void Settings_dialog::save_and_close() {
    check_gnuplot_binary();
    check_exiv2_binary();
    settings.setValue(setting_threshold, threshold_line->text());
    settings.setValue(setting_linear_gamma, cb_linear_gamma->checkState());
    settings.setValue(setting_annotation, cb_annotation->checkState());
    settings.setValue(setting_profile, cb_profile->checkState());
    settings.setValue(setting_grid, cb_grid->checkState());
    settings.setValue(setting_gnuplot, gnuplot_line->text());
    settings.setValue(setting_exiv, exiv_line->text());
    
    send_argument_string();
    
    close();
}

void Settings_dialog::browse_for_gnuplot(void) {
    QString gnuplot = QFileDialog::getOpenFileName(
        this,
        "Locate gnuplot binary",
        QString("/usr/bin/gnuplot"),
        QString::null
    );

    if (gnuplot != QString::null) {
        gnuplot_line->setText(gnuplot);
    }

    check_gnuplot_binary();
}


void Settings_dialog::check_gnuplot_binary(void) {
    bool gnuplot_exists = QFile::exists(get_gnuplot_binary());
    if (!gnuplot_exists) {
        QMessageBox::warning(
            this, 
            QString("gnuplot helper"), 
            QString("gnuplot helper executable not found. Please reconfigure.")
        );
    }
}

void Settings_dialog::browse_for_exiv(void) {
    QString exiv = QFileDialog::getOpenFileName(
        this,
        "Locate exiv2 binary",
        QString("/usr/bin/exiv2"),
        QString::null
    );

    if (exiv != QString::null) {
        exiv_line->setText(exiv);
    }

    check_exiv2_binary();
}

void Settings_dialog::check_exiv2_binary(void) {    
    bool exiv_exists = QFile::exists(get_exiv2_binary());
    if (!exiv_exists) {
        QMessageBox::warning(
            this, 
            QString("Exiv2 helper"), 
            QString("Exiv2 helper executable not found. Please reconfigure.")
        );
    }
}

QString Settings_dialog::get_gnuplot_binary(void) const {
    return gnuplot_line->text();
}

QString Settings_dialog::get_exiv2_binary(void) const {
    return exiv_line->text();
}
