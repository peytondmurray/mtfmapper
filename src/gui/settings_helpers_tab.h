/*
Copyright 2020 Frans van den Bergh. All rights reserved.

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
#ifndef SETTINGS_HELPERS_TAB_H
#define SETTINGS_HELPERS_TAB_H

#include <QtWidgets>

class Settings_helpers_tab : public QWidget {
    Q_OBJECT

  public:
    Settings_helpers_tab(QWidget* parent = nullptr);

    QString get_gnuplot_binary(void) const;
    QString get_exiv2_binary(void) const;
    QString get_dcraw_binary(void) const;
    void check_gnuplot_binary(void);
    void check_exiv2_binary(void);
    void check_dcraw_binary(void);


    QLabel* gnuplot_label;
    QLineEdit* gnuplot_line;
    QPushButton* gnuplot_button;

    QLabel* exiv_label;
    QLineEdit* exiv_line;
    QPushButton* exiv_button;

    QLabel* dcraw_label;
    QLineEdit* dcraw_line;
    QPushButton* dcraw_button;

    static const QString setting_gnuplot;
    static const QString setting_exiv;
    static const QString setting_dcraw;
    static QString setting_gnuplot_default;
    static QString setting_exiv_default;
    static QString setting_dcraw_default;

  private:

  public slots:
    void browse_for_gnuplot();
    void browse_for_exiv();
    void browse_for_dcraw();

};

#endif

