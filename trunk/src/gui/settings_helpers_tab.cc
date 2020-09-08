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

#include "settings_helpers_tab.h"

#include "common.h"

const QString Settings_helpers_tab::setting_gnuplot = "setting_gnuplot";
const QString Settings_helpers_tab::setting_exiv = "setting_exiv";
const QString Settings_helpers_tab::setting_dcraw = "setting_dcraw";
#ifdef _WIN32
QString Settings_helpers_tab::setting_gnuplot_default = "gnuplot.exe";
QString Settings_helpers_tab::setting_exiv_default = "exiv2.exe";
QString Settings_helpers_tab::setting_dcraw_default = "dcraw.exe";
#else
QString Settings_helpers_tab::setting_gnuplot_default = "/usr/bin/gnuplot";
QString Settings_helpers_tab::setting_exiv_default = "/usr/bin/exiv2";
QString Settings_helpers_tab::setting_dcraw_default = "/usr/bin/dcraw";
#endif

Settings_helpers_tab::Settings_helpers_tab(QWidget *parent ATTRIBUTE_UNUSED) {

    gnuplot_label  = new QLabel(tr("gnuplot executable:"), this);
    gnuplot_line   = new QLineEdit(this);
    gnuplot_button = new QPushButton(tr("Browse"), this);

    exiv_label  = new QLabel(tr("exiv2 executable:"), this);
    exiv_line   = new QLineEdit(this);
    exiv_button = new QPushButton(tr("Browse"), this);

    dcraw_label  = new QLabel(tr("dcraw executable:"), this);
    dcraw_line   = new QLineEdit(this);
    dcraw_button = new QPushButton(tr("Browse"), this);

    #ifdef _WIN32
    setting_gnuplot_default = QCoreApplication::applicationDirPath() + QString("/gnuplot/gnuplot.exe");
    setting_exiv_default = QCoreApplication::applicationDirPath() + QString("/exiv2/exiv2.exe");
    setting_dcraw_default = QCoreApplication::applicationDirPath() + QString("/dcraw/dcraw.exe");
    #endif

    QSettings settings("mtfmapper", "mtfmapper");
    gnuplot_line->setText(settings.value(setting_gnuplot, setting_gnuplot_default).toString());
    exiv_line->setText(settings.value(setting_exiv, setting_exiv_default).toString());
    dcraw_line->setText(settings.value(setting_dcraw, setting_dcraw_default).toString());
    
    QGridLayout *helper_layout = new QGridLayout;
    helper_layout->addWidget(gnuplot_label, 0, 0);
    helper_layout->addWidget(gnuplot_line, 0, 1);
    helper_layout->addWidget(gnuplot_button, 0, 2);
    helper_layout->addWidget(exiv_label, 1, 0);
    helper_layout->addWidget(exiv_line, 1, 1);
    helper_layout->addWidget(exiv_button, 1, 2);
    helper_layout->addWidget(dcraw_label, 2, 0);
    helper_layout->addWidget(dcraw_line, 2, 1);
    helper_layout->addWidget(dcraw_button, 2, 2);
    helper_layout->setColumnStretch(1, 1);

    QGroupBox* gb = new QGroupBox("Helper programs");
    gb->setLayout(helper_layout);

    QVBoxLayout* tab_v_layout = new QVBoxLayout;
    tab_v_layout->addWidget(gb);
    tab_v_layout->addStretch(1);

    setLayout(tab_v_layout);
    
    connect(gnuplot_button, SIGNAL(clicked()), this, SLOT( browse_for_gnuplot() ));
    connect(exiv_button, SIGNAL(clicked()), this, SLOT( browse_for_exiv() ));
    connect(dcraw_button, SIGNAL(clicked()), this, SLOT( browse_for_dcraw() ));
}



void Settings_helpers_tab::browse_for_gnuplot(void) {
    QString gnuplot = QFileDialog::getOpenFileName(
        this,
        "Locate gnuplot binary",
        QString("/usr/bin/gnuplot"),
        QString()
    );

    if (gnuplot != QString()) {
        gnuplot_line->setText(gnuplot);
    }

    check_gnuplot_binary();
}

void Settings_helpers_tab::check_gnuplot_binary(void) {
    bool gnuplot_exists = QFile::exists(get_gnuplot_binary());
    if (!gnuplot_exists) {
        QMessageBox::warning(
            this, 
            QString("gnuplot helper"), 
            QString("gnuplot helper executable not found. Please reconfigure.")
        );
    }
}

void Settings_helpers_tab::browse_for_exiv(void) {
    QString exiv = QFileDialog::getOpenFileName(
        this,
        "Locate exiv2 binary",
        QString("/usr/bin/exiv2"),
        QString()
    );

    if (exiv != QString()) {
        exiv_line->setText(exiv);
    }

    check_exiv2_binary();
}

void Settings_helpers_tab::check_exiv2_binary(void) {
    bool exiv_exists = QFile::exists(get_exiv2_binary());
    if (!exiv_exists) {
        QMessageBox::warning(
            this, 
            QString("Exiv2 helper"), 
            QString("Exiv2 helper executable not found. Please reconfigure.")
        );
    }
}

void Settings_helpers_tab::check_dcraw_binary(void) {
    bool dcraw_exists = QFile::exists(get_dcraw_binary());
    if (!dcraw_exists) {
        QMessageBox::warning(
            this, 
            QString("dcraw helper"), 
            QString("dcraw helper executable not found. Please reconfigure.")
        );
    }
}


void Settings_helpers_tab::browse_for_dcraw(void) {
    QString dcraw = QFileDialog::getOpenFileName(
        this,
        "Locate dcraw binary",
        QString("/usr/bin/dcraw"),
        QString()
    );

    if (dcraw != QString()) {
        dcraw_line->setText(dcraw);
    }

    check_exiv2_binary();
}


QString Settings_helpers_tab::get_gnuplot_binary(void) const {
    return gnuplot_line->text();
}

QString Settings_helpers_tab::get_exiv2_binary(void) const {
    return exiv_line->text();
}

QString Settings_helpers_tab::get_dcraw_binary(void) const {
    return dcraw_line->text();
}

