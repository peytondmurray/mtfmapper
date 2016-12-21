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
#ifndef SETTINGS_DIALOG_H
#define SETTINGS_DIALOG_H

#include <QDialog>
#include <QSettings>

class QLabel;
class QLineEdit;
class QPushButton;
class QCheckBox;
class QRadioButton;

class Settings_dialog : public QDialog 
{
  Q_OBJECT
  
  public:
    Settings_dialog(QWidget *parent);
    void send_argument_string(void);
    QString get_gnuplot_binary(void) const;
    QString get_exiv2_binary(void) const;
    QString get_dcraw_binary(void) const;
    void check_gnuplot_binary(void);
    void check_exiv2_binary(void);
    void check_dcraw_binary(void);
    
    QSettings   settings;
    
    QLabel*     arguments_label;
    QLineEdit*  arguments_line;
    QLabel*     threshold_label;
    QLineEdit*  threshold_line;
    QLabel*     pixsize_label;
    QLineEdit*  pixsize_line;
    QPushButton* accept_button;
    QPushButton* cancel_button;
    
    QCheckBox*  cb_linear_gamma;
    QCheckBox*  cb_annotation;
    QCheckBox*  cb_profile;
    QCheckBox*  cb_grid;
    QCheckBox*  cb_focus;
    QCheckBox*  cb_lensprofile;
    QCheckBox*  cb_autocrop;
    QCheckBox*  cb_lpmm;
    QCheckBox*  cb_gnuplot_scaled;
    
    QRadioButton* rb_colour_none;
    QRadioButton* rb_colour_red;
    QRadioButton* rb_colour_green;
    QRadioButton* rb_colour_blue;

    QLabel*     gnuplot_label;
    QLineEdit*  gnuplot_line;
    QPushButton* gnuplot_button;

    QLabel*     exiv_label;
    QLineEdit*  exiv_line;
    QPushButton* exiv_button;

    QLabel*     dcraw_label;
    QLineEdit*  dcraw_line;
    QPushButton* dcraw_button;
    
    int gnuplot_img_width;
    
  signals:
    void argument_string(QString s);  
    
  public slots:
    void open();
    void save_and_close();
    void browse_for_gnuplot();
    void browse_for_exiv();
    void browse_for_dcraw();
    void set_gnuplot_img_width(int w);
};

#endif

