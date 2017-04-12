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
#ifndef SFR_DIALOG_H
#define SFR_DIALOG_H

#include <QDialog>
#include <QtCharts>
using namespace QtCharts;

#include "sfr_entry.h"
#include "mouse_chart.h"

#include <vector>
using std::vector;

class Sfr_dialog : public QDialog {
  Q_OBJECT
  
  public:
    Sfr_dialog(QWidget *parent, const Sfr_entry& entry);
    void replace_entry(const Sfr_entry& entry);
    void add_entry(const Sfr_entry& entry);
    void notify_mouse_position(double value);
  
  protected:
    void paintEvent(QPaintEvent* event);
    void reject(void);

  private:
    void populate_series(const Sfr_entry& entry, QLineSeries* s);
  
    vector<Sfr_entry> entries;
    QChart* chart;
    Mouse_chart* chart_view;
    QValueAxis* x_axis;
    QValueAxis* y_axis;
    vector<QLineSeries*> series; 
    QGraphicsSimpleTextItem* mtf50_text;
    double cursor_domain_value;
    QLabel* cursor_label;
    
  public slots:
    
};

// least-squares fit (inverse of design matrix) of a cubic polynomial through 4 points [-1..2]/64.0
const double sfr_cubic_weights[4][4] = {
  { 0.00000000000000e+00,   1.00000000000000e+00,   0.00000000000000e+00,   0.00000000000000e+00},
  {-2.13333333333333e+01,  -3.20000000000000e+01,   6.40000000000000e+01,  -1.06666666666667e+01},
  { 2.04800000000000e+03,  -4.09600000000000e+03,   2.04800000000000e+03,   0.00000000000000e+00},
  {-4.36906666666667e+04,   1.31072000000000e+05,  -1.31072000000000e+05,   4.36906666666667e+04}
};

#endif

