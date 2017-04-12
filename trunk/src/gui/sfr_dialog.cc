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
#include "sfr_dialog.h"
#include "sfr_dialog.moc"

#include "common.h"
#include "config.h"


Sfr_dialog::Sfr_dialog(QWidget* parent ATTRIBUTE_UNUSED, const Sfr_entry& entry) {

    series.push_back(new QLineSeries());
    
    double maxval = 0;
    for (size_t i=0; i < 64; i++) {
        maxval = std::max(entry.sfr[i], maxval);
    }
    populate_series(entry, series[0]);
    
    chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(series[0]);
    
    x_axis = new QValueAxis();
    x_axis->setRange(0, 1.0);
    x_axis->setTickCount(11);
    x_axis->setTitleText("Frequency (c/p)");
    x_axis->setLabelFormat("%3.1f");
    chart->addAxis(x_axis, Qt::AlignBottom);
    
    y_axis = new QValueAxis();
    y_axis->setRange(0, maxval);
    y_axis->setTitleText("Contrast");
    chart->addAxis(y_axis, Qt::AlignLeft);
    
    series.back()->attachAxis(x_axis);
    series.back()->attachAxis(y_axis);
    
    chart_view = new Mouse_chart(chart, this);
    chart_view->setRenderHint(QPainter::Antialiasing);
    
    cursor_label = new QLabel(this);
    cursor_label->setText("Nothing near");
    cursor_label->setAlignment(Qt::AlignHCenter);
    cursor_label->setStyleSheet("QLabel { color : blue; }");
    
    QGridLayout* hlayout = new QGridLayout;
    hlayout->addWidget(chart_view, 0, 0);
    hlayout->addWidget(cursor_label, 1, 0);
    setLayout(hlayout);
    
    mtf50_text = new QGraphicsSimpleTextItem(chart);
    mtf50_text->setPen(QPen(series[0]->pen().color()));
   
    resize(600, 300);
    setMinimumHeight(300);
    setMinimumWidth(600);
    setWindowTitle("SFR / MTF curve");

    show();
}

void Sfr_dialog::reject(void) {
    chart->removeAllSeries();
    series.clear();
    QDialog::reject();
}

void Sfr_dialog::paintEvent(QPaintEvent* event) {
    

    QVector<QPointF> pts = series.front()->pointsVector();
    
    double prev_val  = pts[0].y();
    double mtf50 = 0;
    
    bool done = false;
    for (int i=0; i < pts.size() && !done; i++) {
        double mag = pts[i].y();
        if (prev_val > 0.5 && mag <= 0.5) {
            mtf50 = pts[i].x();
            done = true;
        }
        prev_val = mag;
    }
    
    double w = chart->size().width();
    
    QFontMetrics fm(QWidget::fontMetrics());
    char mtf50_str[200];
    sprintf(mtf50_str, "MTF50=%.3lf   ", mtf50);
    int tw = fm.width(mtf50_str);
    int th = fm.height();
    
    mtf50_text->setPos(w - tw, th); 
    if (mtf50 < 1.0) {
        mtf50_text->setText(mtf50_str);
    } else {
        mtf50_text->setText("");
    }
    
    double contrast = 0;
    for (int i=0; i < pts.size() && pts[i].x() < cursor_domain_value; i++) {
        contrast = pts[i].y();
    }
    sprintf(mtf50_str, "frequency: %5.3lf contrast: %5.3lf", cursor_domain_value, contrast);
    cursor_label->setText(mtf50_str);
    
    QDialog::paintEvent(event);
}

void Sfr_dialog::replace_entry(const Sfr_entry& entry) {
    if (series.size() == 0) {
        series.push_back(new QLineSeries());
    } else {
        chart->removeSeries(series.back());
        series.back() = new QLineSeries();
    }
    
    if (series.size() == 1) { // only recompute maxval for the first curve
        double maxval = 0;
        for (size_t i=0; i < 64; i++) {
            maxval = std::max(entry.sfr[i], maxval);
        }
        y_axis->setRange(0, maxval);
    }
    populate_series(entry, series.back());
    
    chart->addSeries(series.back());
    series.back()->attachAxis(x_axis);
    series.back()->attachAxis(y_axis);
    update();
    show();
    raise();
    activateWindow();
}

void Sfr_dialog::add_entry(const Sfr_entry& entry) {
    if (series.size() < 3) {
        series.push_back(new QLineSeries());
        
        populate_series(entry, series.back());
        
        chart->addSeries(series.back());
        series.back()->attachAxis(x_axis);
        series.back()->attachAxis(y_axis);
    }
    update();
    show();
    raise();
    activateWindow();
}

void Sfr_dialog::populate_series(const Sfr_entry& entry, QLineSeries* s) {
    for (size_t i=0; i < 64; i++) {
        double coef[4] = {0,0,0,0};
        for (int si=int(i)-1,ri=0; si <= int(i) + 2; si++,ri++) {
            int ei = si < 0 ? 0 : (si > 63 ? 63 : si); // assume MTF curve is constant outside of range
            
            for (int c=0; c < 4; c++) {
                coef[c] += entry.sfr[ei] * sfr_cubic_weights[c][ri];
            }
        }
        // now we can evaluate points at position x in [0,1) using the fitted polynomial
        for (int xi=0; xi < 20; xi++) {
            double dx = xi*(1.0/(20.0*64.0));
            double iy = coef[0] + coef[1]*dx + coef[2]*dx*dx + coef[3]*dx*dx*dx;
            s->append(i*(1.0/64.0) + dx, iy);
        }
    }
}

void Sfr_dialog::notify_mouse_position(double value) {
    cursor_domain_value = value;
    update();
}