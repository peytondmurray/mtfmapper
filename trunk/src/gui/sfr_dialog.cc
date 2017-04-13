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

static double multiple(double x) {
    double rval = 0;
    while (rval < (x+5e-4)) {
        rval += 0.2;
    }
    return rval;
}


Sfr_dialog::Sfr_dialog(QWidget* parent ATTRIBUTE_UNUSED, const Sfr_entry& entry) : cursor_domain_value(0) {

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
    
    double roundup_max = multiple(maxval);
    y_axis = new QValueAxis();
    y_axis->setRange(0, roundup_max);
    y_axis->setTickCount(roundup_max*5 + 1);
    y_axis->setTitleText("Contrast");
    chart->addAxis(y_axis, Qt::AlignLeft);
    
    series.back()->attachAxis(x_axis);
    series.back()->attachAxis(y_axis);
    
    chart_view = new Sfr_chartview(chart, this);
    chart_view->setRenderHint(QPainter::Antialiasing);
    
    label_layout = new QGridLayout();
    for (int i=0; i < 4; i++) {
        cursor_label.push_back(new QLabel(this));
        cursor_label.back()->setAlignment(i == 0 ? Qt::AlignRight : Qt::AlignLeft);
        cursor_label.back()->setStyleSheet("font-weight: bold;");
    }
    label_layout->addWidget(cursor_label[0], 0, 0);
    label_layout->addWidget(cursor_label[1], 0, 1);
    label_layout->addWidget(cursor_label[2], 1, 1);
    label_layout->addWidget(cursor_label[3], 2, 1);
    cursor_label[2]->setMinimumHeight(0);
    cursor_label[3]->setMinimumHeight(0);
    cursor_label[2]->hide();
    cursor_label[3]->hide();
    label_layout->rowMinimumHeight(0);
    label_layout->setVerticalSpacing(1);
    
    QGroupBox* gbox = new QGroupBox("");
    gbox->setLayout(label_layout);
    
    
    QVBoxLayout* hlayout = new QVBoxLayout;
    hlayout->addWidget(chart_view);
    hlayout->addWidget(gbox);
    setLayout(hlayout);
    
    for (int i=0; i < 3; i++) {
        mtf50_text.push_back(new QGraphicsSimpleTextItem(chart));
    }
    
    mtf50_rect = new QGraphicsRectItem(chart);
    mtf50_rect->setRect(0,0,0,0);
    QColor charcoal(54, 69, 79, 100);
    mtf50_rect->setBrush(QBrush(charcoal));
    mtf50_rect->setPen(charcoal);
    
    setAutoFillBackground(true); 
    QPalette palette;
    palette.setColor(backgroundRole(), Qt::white);
    setPalette(palette);
    
    chart->resize(650, 350);
    setMinimumHeight(400);
    setMinimumWidth(750);
    setWindowTitle("SFR / MTF curve");

    show();
}

void Sfr_dialog::clear(void) {
    chart->removeAllSeries();
    series.clear();
}

void Sfr_dialog::reject(void) {
    clear();
    QDialog::reject();
    for (auto l: cursor_label) {
        l->setText("");
        l->hide();
    }
    for (auto m: mtf50_text) {
        m->setText("");
    }
}

void Sfr_dialog::paintEvent(QPaintEvent* event) {
    QDialog::paintEvent(event);
    
    vector<double> contrast_list;
    vector<double> mtf50_list;
    if (series.size() > 0) {

        for (size_t si=0; si < series.size(); si++) {
        
            QVector<QPointF> pts = series[si]->pointsVector();

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
            mtf50_list.push_back(mtf50);
            
            double contrast = pts[0].y();
            for (int i = 0; i < pts.size() && pts[i].x() < cursor_domain_value; i++) {
                contrast = pts[i].y();
            }
            
            contrast_list.push_back(contrast);
        }
        

        QFontMetrics fm(QWidget::fontMetrics());
        int th = fm.height();
        char mtf50_str[200];
        int text_end = chart->size().width() - 2*fm.width("m");
        
        // add mtf50 tags in reverse
        for (int mi=mtf50_list.size()-1; mi >= 0; mi--) {
            if (mtf50_list[mi] < 1.0) {
                sprintf(mtf50_str, "MTF50=%.3lf", mtf50_list[mi]);
            } else {
                sprintf(mtf50_str, "MTF50=N/A");
            }
            
            int tw = fm.width(mtf50_str) + 2*fm.width("m");
            
            mtf50_text[mi]->setBrush(series[mi]->pen().color());
            mtf50_text[mi]->setPen(QPen(series[mi]->pen().color()));
            mtf50_text[mi]->setPos(text_end - tw, th);
            mtf50_text[mi]->setText(mtf50_str);
            text_end -= tw;
        }
    }
    
    QPointF tpos(cursor_domain_value, y_axis->max() - 0.001);
    QPointF bpos(cursor_domain_value, 0.001);
    tpos = chart->mapToPosition(tpos);
    bpos = chart->mapToPosition(bpos);
    tpos.rx() = std::max(tpos.x()-1, chart->mapToPosition(QPointF(0, 0)).x());
    bpos.rx() = std::min(bpos.x()+1, chart->mapToPosition(QPointF(1.0, 0)).x());
    
    mtf50_rect->setRect(tpos.x(), tpos.y(), bpos.x() - tpos.x(), bpos.y() - tpos.y());
    
    char mtf_str[200];
    sprintf(mtf_str, "frequency: %5.3lf ", cursor_domain_value); 
    cursor_label[0]->setText(mtf_str);
    cursor_label[0]->show();
    
    char contrast_str[200];
    for (size_t mi=0; mi < contrast_list.size(); mi++) {
        sprintf(contrast_str, "contrast: %5.3lf ", contrast_list[mi]);
        cursor_label[mi+1]->setText(contrast_str);
        
        QPalette palette = cursor_label[mi+1]->palette();
        palette.setColor(cursor_label[mi+1]->foregroundRole(), series[mi]->pen().color());
        cursor_label[mi+1]->setPalette(palette);
        cursor_label[mi+1]->show();
    }
    
    chart->update(); // this causes some lag, but it elliminates exposed cruft. a better solution would be nice
}

void Sfr_dialog::replace_entry(const Sfr_entry& entry) {
    if (series.size() == 0) {
        series.push_back(new QLineSeries());
    } else {
        chart->removeSeries(series.back());
        series.back() = new QLineSeries();
    }
    
    double maxval = 0;
    for (size_t i=0; i < 64; i++) {
        maxval = std::max(entry.sfr[i], maxval);
    }
    double roundup_max = multiple(maxval);
    y_axis->setRange(0, roundup_max);
    y_axis->setTickCount(roundup_max*5 + 1);
        
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
    } 
    replace_entry(entry);
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

void Sfr_dialog::notify_mouse_position(double value) { // we still need this to update the bottom label, but we could merge this in the chart?
    cursor_domain_value = std::min(1.0, std::max(0.0, value));
    
    QPointF top(cursor_domain_value, 1);
    QPointF bottom(cursor_domain_value, 0);
    
    top = chart->mapToPosition(top);
    bottom = chart->mapToPosition(bottom);
    
    update();
}