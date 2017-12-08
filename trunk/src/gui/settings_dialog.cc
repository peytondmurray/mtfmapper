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
#include "settings_dialog.h"
#include "settings_dialog.moc"
#include "nonempty_validator.h"

#include "common.h"

#include <iostream>
using std::cout;
using std::endl;

const QString setting_threshold = "setting_threshold_sauvola";
const QString setting_threshold_default = "0.55";
const QString setting_pixsize = "setting_pixelsize";
const QString setting_pixsize_default = "4.78";
const QString setting_bayer = "setting_bayer";
const QString setting_linear_gamma = "setting_gamma";
const Qt::CheckState setting_linear_gamma_default = Qt::Unchecked;
const QString setting_annotation = "setting_annotation";
const Qt::CheckState setting_annotation_default = Qt::Checked;
const QString setting_profile = "setting_profile";
const Qt::CheckState setting_profile_default = Qt::Checked;
const QString setting_grid = "setting_grid";
const Qt::CheckState setting_grid_default = Qt::Checked;
const QString setting_focus = "setting_focus";
const Qt::CheckState setting_focus_default = Qt::Unchecked;
const QString setting_lensprofile = "setting_lensprofile";
const Qt::CheckState setting_lensprofile_default = Qt::Unchecked;
const QString setting_orientation = "setting_orientation";
const Qt::CheckState setting_orientation_default = Qt::Unchecked;
const QString setting_autocrop = "setting_autocrop";
const Qt::CheckState setting_autocrop_default = Qt::Unchecked;
const QString setting_lpmm = "setting_lpmm";
const Qt::CheckState setting_lpmm_default = Qt::Unchecked;
const QString setting_gnuplot_scaled = "setting_gnuplot_scaled";
const Qt::CheckState setting_gnuplot_scaled_default = Qt::Checked;
const QString setting_zscale = "setting_zscale";
const int setting_zscale_default = 0;
const QString setting_ea_f = "setting_equiangular_f";
const QString setting_ea_f_default = "16.0";
const QString setting_sg_f = "setting_stereographic_f";
const QString setting_sg_f_default = "8.0";
const QString setting_gnuplot = "setting_gnuplot";
const QString setting_exiv = "setting_exiv";
const QString setting_dcraw = "setting_dcraw";
const QString setting_lens = "setting_lens_correction_type";
#ifdef _WIN32
static QString setting_gnuplot_default = "gnuplot.exe";
static QString setting_exiv_default = "exiv2.exe";
static QString setting_dcraw_default = "dcraw.exe";
#else
const QString setting_gnuplot_default = "/usr/bin/gnuplot";
const QString setting_exiv_default = "/usr/bin/exiv2";
const QString setting_dcraw_default = "/usr/bin/dcraw";
#endif

Settings_dialog::Settings_dialog(QWidget *parent ATTRIBUTE_UNUSED)
 : settings("mtfmapper", "mtfmapper"), gnuplot_img_width(1024)
{
    arguments_label = new QLabel(tr("Arguments:"), this);
    arguments_line  = new QLineEdit(this);
    
    QDoubleValidator* dv_thresh = new Nonempty_DoubleValidator(0.001, 1.0, 3, 0.55, this);
    threshold_label = new QLabel(tr("Threshold:"), this);
    threshold_line  = new QLineEdit(this);
    threshold_line->setValidator(dv_thresh);
    
    QDoubleValidator* dv_pixsize = new Nonempty_DoubleValidator(0.001, 999.0, 3, 4.0, this);
    pixsize_label   = new QLabel(tr("Pixel size:"), this);
    pixsize_line    = new QLineEdit(this);
    pixsize_line->setValidator(dv_pixsize);

    gnuplot_label  = new QLabel(tr("gnuplot executable:"), this);
    gnuplot_line   = new QLineEdit(this);
    gnuplot_button = new QPushButton(tr("Browse"), this);

    exiv_label  = new QLabel(tr("exiv2 executable:"), this);
    exiv_line   = new QLineEdit(this);
    exiv_button = new QPushButton(tr("Browse"), this);

    dcraw_label  = new QLabel(tr("dcraw executable:"), this);
    dcraw_line   = new QLineEdit(this);
    dcraw_button = new QPushButton(tr("Browse"), this);

    accept_button = new QPushButton(tr("&Accept"), this);
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    cancel_button = new QPushButton(tr("&Cancel"), this);
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
    cb_linear_gamma = new QCheckBox("Linear gamma (8 bit)", this);
    cb_annotation   = new QCheckBox("Annotation", this);
    cb_profile      = new QCheckBox("Profile", this);
    cb_grid         = new QCheckBox("Grid", this);
    cb_focus        = new QCheckBox("Focus position", this);
    cb_lensprofile  = new QCheckBox("Lens profile", this);
    cb_orientation  = new QCheckBox("Chart orientation", this);
    cb_autocrop     = new QCheckBox("Autocrop", this);
    cb_lpmm         = new QCheckBox("Line pairs/mm units", this);
    cb_gnuplot_scaled = new QCheckBox("Scale plots to window", this);
    
    rb_colour_none  = new QRadioButton("none");
    rb_colour_red   = new QRadioButton("red");
    rb_colour_green = new QRadioButton("green");
    rb_colour_blue  = new QRadioButton("blue");
    
    rb_lens_pw_quad  = new QRadioButton("piecewise-quadratic");
    rb_lens_quad  = new QRadioButton("quadratic");
    rb_lens_none  = new QRadioButton("none");
    rb_lens_radial = new QRadioButton("radial");
    rb_lens_equiangular = new QRadioButton("equiangular");
    rb_lens_stereo = new QRadioButton("stereographic");
    
    QDoubleValidator* dv_f = new Nonempty_DoubleValidator(0.01, 999.0, 2, 16.0, this);
    ea_f_label = new QLabel(tr("focal length"), this);
    ea_f_line = new QLineEdit(this);
    ea_f_line->setMaxLength(5);
    ea_f_line->setText(settings.value(setting_ea_f, setting_ea_f_default).toString());
    ea_f_line->setValidator(dv_f);
    sg_f_label = new QLabel(tr("focal length"), this);
    sg_f_line = new QLineEdit(this);
    sg_f_line->setMaxLength(5);
    sg_f_line->setText(settings.value(setting_sg_f, setting_sg_f_default).toString());
    sg_f_line->setValidator(dv_f);
    
    threshold_line->setText(settings.value(setting_threshold, setting_threshold_default).toString());
    pixsize_line->setText(settings.value(setting_pixsize, setting_pixsize_default).toString());
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
    cb_focus->setCheckState(
        (Qt::CheckState)settings.value(setting_focus, setting_focus_default).toInt()
    );
    cb_lensprofile->setCheckState(
        (Qt::CheckState)settings.value(setting_lensprofile, setting_lensprofile_default).toInt()
    );
    cb_orientation->setCheckState(
        (Qt::CheckState)settings.value(setting_orientation, setting_orientation_default).toInt()
    );
    cb_autocrop->setCheckState(
        (Qt::CheckState)settings.value(setting_autocrop, setting_autocrop_default).toInt()
    );
    cb_lpmm->setCheckState(
        (Qt::CheckState)settings.value(setting_lpmm, setting_lpmm_default).toInt()
    );
    cb_gnuplot_scaled->setCheckState(
        (Qt::CheckState)settings.value(setting_gnuplot_scaled, setting_gnuplot_scaled_default).toInt()
    );
    
    switch(settings.value(setting_bayer, 0).toInt()) {
        case 0: rb_colour_none->setChecked(true); rb_colour_red->setChecked(false); rb_colour_green->setChecked(false); rb_colour_blue->setChecked(false); break;
        case 1: rb_colour_none->setChecked(false); rb_colour_red->setChecked(true); rb_colour_green->setChecked(false); rb_colour_blue->setChecked(false); break;
        case 2: rb_colour_none->setChecked(false); rb_colour_red->setChecked(false); rb_colour_green->setChecked(true); rb_colour_blue->setChecked(false); break;
        case 3: rb_colour_none->setChecked(false); rb_colour_red->setChecked(false); rb_colour_green->setChecked(false); rb_colour_blue->setChecked(true); break;
    }
    
    rb_lens_pw_quad->setChecked(false);
    rb_lens_quad->setChecked(false);
    rb_lens_none->setChecked(false);
    rb_lens_radial->setChecked(false);
    rb_lens_equiangular->setChecked(false);
    rb_lens_stereo->setChecked(false);
    
    switch(settings.value(setting_lens, 0).toInt()) {
        case 0: rb_lens_pw_quad->setChecked(true); break;
        case 1: rb_lens_quad->setChecked(true); break;
        case 2: rb_lens_none->setChecked(true); break;
        case 3: rb_lens_radial->setChecked(true); break;
        case 4: rb_lens_equiangular->setChecked(true); break;
        case 5: rb_lens_stereo->setChecked(true); break;
    }

    #ifdef _WIN32
    setting_gnuplot_default = QCoreApplication::applicationDirPath() + QString("/gnuplot/gnuplot.exe");
    setting_exiv_default = QCoreApplication::applicationDirPath() + QString("/exiv2/exiv2.exe");
    setting_dcraw_default = QCoreApplication::applicationDirPath() + QString("/dcraw/dcraw.exe");
    #endif

    gnuplot_line->setText(settings.value(setting_gnuplot, setting_gnuplot_default).toString());
    exiv_line->setText(settings.value(setting_exiv, setting_exiv_default).toString());
    dcraw_line->setText(settings.value(setting_dcraw, setting_dcraw_default).toString());
    
    zscale_label = new QLabel("3D plot z-axis relative scale factor", this);
    zscale_slider = new QSlider(Qt::Horizontal, this);
    zscale_slider->setFocusPolicy(Qt::StrongFocus);
    zscale_slider->setTickPosition(QSlider::TicksAbove);
    zscale_slider->setMinimum(0);
    zscale_slider->setMaximum(20);
    zscale_slider->setTickInterval(5);
    zscale_slider->setSingleStep(1);
    zscale_slider->setValue(settings.value(setting_zscale, setting_zscale_default).toInt());

    QGroupBox* voGroupBox = new QGroupBox(tr("Output types"), this);
    QVBoxLayout *vo_layout = new QVBoxLayout(this);
    vo_layout->addWidget(cb_annotation);
    vo_layout->addWidget(cb_profile);
    vo_layout->addWidget(cb_grid);
    vo_layout->addWidget(cb_focus);
    vo_layout->addWidget(cb_lensprofile);
    vo_layout->addWidget(cb_orientation);
    voGroupBox->setLayout(vo_layout);
    
    QGroupBox* v2GroupBox = new QGroupBox(tr("Flags"), this);
    QVBoxLayout *cb_layout = new QVBoxLayout(this);
    cb_layout->addWidget(cb_linear_gamma);
    cb_layout->addWidget(cb_autocrop);
    cb_layout->addWidget(cb_lpmm);
    cb_layout->addWidget(cb_gnuplot_scaled);
    v2GroupBox->setLayout(cb_layout);
    
    QGroupBox* v4GroupBox = new QGroupBox(tr("Bayer channel"), this);
    QVBoxLayout *rb_layout = new QVBoxLayout(this);
    rb_layout->addWidget(rb_colour_none);
    rb_layout->addWidget(rb_colour_red);
    rb_layout->addWidget(rb_colour_green);
    rb_layout->addWidget(rb_colour_blue);
    v4GroupBox->setLayout(rb_layout);

    QGroupBox* v3GroupBox = new QGroupBox(tr("Helpers"), this);
    QGridLayout *helper_layout = new QGridLayout(this);
    helper_layout->addWidget(gnuplot_label, 0, 0);
    helper_layout->addWidget(gnuplot_line, 1, 0);
    helper_layout->addWidget(gnuplot_button, 1, 1);
    helper_layout->addWidget(exiv_label, 2, 0);
    helper_layout->addWidget(exiv_line, 3, 0);
    helper_layout->addWidget(exiv_button, 3, 1);
    helper_layout->addWidget(dcraw_label, 4, 0);
    helper_layout->addWidget(dcraw_line, 5, 0);
    helper_layout->addWidget(dcraw_button, 5, 1);
    v3GroupBox->setLayout(helper_layout);
    
    QGroupBox* lens = new QGroupBox(tr("Lens distortion correction"), this);
    QGridLayout* lens_layout = new QGridLayout(this);
    lens_layout->addWidget(rb_lens_pw_quad, 0, 0);
    lens_layout->addWidget(rb_lens_quad, 1, 0);
    lens_layout->addWidget(rb_lens_none, 2, 0);
    lens_layout->addWidget(rb_lens_radial, 3, 0);
    lens_layout->addWidget(rb_lens_equiangular, 4, 0);
    lens_layout->addWidget(rb_lens_stereo, 5, 0);
    lens_layout->addWidget(ea_f_label, 4, 1);
    lens_layout->addWidget(ea_f_line, 4, 2);
    lens_layout->addWidget(sg_f_label, 5, 1);
    lens_layout->addWidget(sg_f_line, 5, 2);
    lens->setLayout(lens_layout);

    QGroupBox* advanced = new QGroupBox(tr("Advanced"), this);
    QGridLayout* adv_layout = new QGridLayout(this);
    adv_layout->addWidget(threshold_label, 0, 0);
    adv_layout->addWidget(threshold_line, 0, 1);
    adv_layout->addWidget(pixsize_label, 1, 0);
    adv_layout->addWidget(pixsize_line, 1, 1);
    adv_layout->addWidget(zscale_label, 2, 0, 1, 2);
    adv_layout->addWidget(zscale_slider, 3, 0, 1, 2);
    adv_layout->addWidget(arguments_label, 4, 0);
    adv_layout->addWidget(arguments_line, 4, 1);
    advanced->setLayout(adv_layout);

    
    QGroupBox* vGroupBox = new QGroupBox(tr("Settings"), this);
    QGridLayout* vlayout = new QGridLayout(this);
    vlayout->addWidget(voGroupBox, 0, 0, 1, 2);
    vlayout->addWidget(v2GroupBox, 1, 0, 1, 2);
    vlayout->addWidget(v4GroupBox, 2, 0, 1, 2);
    vlayout->addWidget(v3GroupBox, 0, 3, 1, 2);
    vlayout->addWidget(lens, 1, 3, 1, 2);
    vlayout->addWidget(advanced, 2, 3, 1, 2);
    vlayout->addWidget(accept_button, 3, 0);
    vlayout->addWidget(cancel_button, 3, 1);
    vGroupBox->setLayout(vlayout);
    
    connect(accept_button, SIGNAL(clicked()), this, SLOT( save_and_close() ));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(gnuplot_button, SIGNAL(clicked()), this, SLOT( browse_for_gnuplot() ));
    connect(exiv_button, SIGNAL(clicked()), this, SLOT( browse_for_exiv() ));
    connect(dcraw_button, SIGNAL(clicked()), this, SLOT( browse_for_dcraw() ));
    connect(rb_lens_equiangular, SIGNAL(clicked()), this, SLOT( equiangular_toggled() ));
    connect(rb_lens_stereo, SIGNAL(clicked()), this, SLOT( stereographic_toggled() ));
    
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
        args = args + QString(" -q");
    }
    
    if (cb_profile->checkState()) {
        args = args + QString(" -p");
    }
    
    if (cb_grid->checkState()) {
        args = args + QString(" -s");
    }
    
    if (cb_focus->checkState()) {
        args = args + QString(" --focus");
    }
    
    if (cb_lensprofile->checkState()) {
        args = args + QString(" --lensprofile");
    }
    
    if (cb_orientation->checkState()) {
        args = args + QString(" --chart-orientation");
    }
    
    if (cb_autocrop->checkState()) {
        args = args + QString(" --autocrop");
    }

    if (cb_lpmm->checkState()) {
        args = args + QString(" --pixelsize " + pixsize_line->text());
    }
    
    if (rb_colour_red->isChecked()) {
        args = args + QString(" --bayer red");
    }
    if (rb_colour_green->isChecked()) {
        args = args + QString(" --bayer green");
    }
    if (rb_colour_blue->isChecked()) {
        args = args + QString(" --bayer blue");
    }
    
    
    if (rb_lens_pw_quad->isChecked()) {
        args = args + QString(" --esf-sampler piecewise-quadratic");
    }
    if (rb_lens_quad->isChecked()) {
        args = args + QString(" --esf-sampler quadratic");
    }
    if (rb_lens_none->isChecked()) {
        args = args + QString(" --esf-sampler line");
    }
    if (rb_lens_radial->isChecked()) {
        args = args + QString(" --optimize-distortion");
    }
    if (rb_lens_equiangular->isChecked()) {
        args = args + QString(" --equiangular %1").arg(ea_f_line->text());
    }
    if (rb_lens_stereo->isChecked()) {
        args = args + QString(" --stereographic %1").arg(sg_f_line->text());
    }
    
    if (cb_gnuplot_scaled->checkState()) {
        args = args + QString(" --gnuplot-width %1").arg(gnuplot_img_width);
    }
    args = args + QString(" --zscale %1").arg(zscale_slider->value()/20.0);
    
    args = args + QString(" %1").arg(arguments_line->text());
    
    emit argument_string(args);
}

void Settings_dialog::save_and_close() {
    check_gnuplot_binary();
    check_exiv2_binary();
    settings.setValue(setting_threshold, threshold_line->text());
    settings.setValue(setting_pixsize, pixsize_line->text());
    settings.setValue(setting_linear_gamma, cb_linear_gamma->checkState());
    settings.setValue(setting_annotation, cb_annotation->checkState());
    settings.setValue(setting_profile, cb_profile->checkState());
    settings.setValue(setting_lpmm, cb_lpmm->checkState());
    settings.setValue(setting_grid, cb_grid->checkState());
    settings.setValue(setting_focus, cb_focus->checkState());
    settings.setValue(setting_lensprofile, cb_lensprofile->checkState());
    settings.setValue(setting_orientation, cb_orientation->checkState());
    settings.setValue(setting_autocrop, cb_autocrop->checkState());
    settings.setValue(setting_gnuplot_scaled, cb_gnuplot_scaled->checkState());
    settings.setValue(setting_gnuplot, gnuplot_line->text());
    settings.setValue(setting_exiv, exiv_line->text());
    settings.setValue(setting_dcraw, dcraw_line->text());
    settings.setValue(setting_zscale, zscale_slider->value());
    
    if (rb_colour_none->isChecked()) {
        settings.setValue(setting_bayer, 0);
    }
    if (rb_colour_red->isChecked()) {
        settings.setValue(setting_bayer, 1);
    }
    if (rb_colour_green->isChecked()) {
        settings.setValue(setting_bayer, 2 );
    }
    if (rb_colour_blue->isChecked()) {
        settings.setValue(setting_bayer, 3);
    }

    if (rb_lens_pw_quad->isChecked()) {
        settings.setValue(setting_lens, 0);
    }
    if (rb_lens_quad->isChecked()) {
        settings.setValue(setting_lens, 1);
    }    
    if (rb_lens_none->isChecked()) {
        settings.setValue(setting_lens, 2);
    }
    if (rb_lens_radial->isChecked()) {
        settings.setValue(setting_lens, 3);
    }
    if (rb_lens_equiangular->isChecked()) {
        settings.setValue(setting_lens, 4);
        cb_lpmm->setCheckState(Qt::Checked);
        settings.setValue(setting_lpmm, cb_lpmm->checkState());
    }
    if (rb_lens_stereo->isChecked()) {
        settings.setValue(setting_lens, 5);
        cb_lpmm->setCheckState(Qt::Checked);
        settings.setValue(setting_lpmm, cb_lpmm->checkState());
    }
    
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

void Settings_dialog::check_dcraw_binary(void) {
    bool dcraw_exists = QFile::exists(get_dcraw_binary());
    if (!dcraw_exists) {
        QMessageBox::warning(
            this, 
            QString("dcraw helper"), 
            QString("dcraw helper executable not found. Please reconfigure.")
        );
    }
}


void Settings_dialog::browse_for_dcraw(void) {
    QString dcraw = QFileDialog::getOpenFileName(
        this,
        "Locate dcraw binary",
        QString("/usr/bin/dcraw"),
        QString::null
    );

    if (dcraw != QString::null) {
        dcraw_line->setText(dcraw);
    }

    check_exiv2_binary();
}

void Settings_dialog::set_gnuplot_img_width(int w) {
    gnuplot_img_width = w < 1024 ? 1024 : w;
}

QString Settings_dialog::get_gnuplot_binary(void) const {
    return gnuplot_line->text();
}

QString Settings_dialog::get_exiv2_binary(void) const {
    return exiv_line->text();
}

QString Settings_dialog::get_dcraw_binary(void) const {
    return dcraw_line->text();
}

void Settings_dialog::equiangular_toggled() {
    if (rb_lens_equiangular->isChecked() && cb_lpmm->checkState() == Qt::Unchecked) {
        cb_lpmm->setCheckState(Qt::Checked);
    }
}

void Settings_dialog::stereographic_toggled() {
    if (rb_lens_stereo->isChecked() && cb_lpmm->checkState() == Qt::Unchecked) {
        cb_lpmm->setCheckState(Qt::Checked);
    }
}