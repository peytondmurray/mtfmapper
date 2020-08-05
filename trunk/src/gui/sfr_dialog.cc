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
#include "include/logger.h"
#include <QtWidgets> 
#include "sfr_dialog.h"
#include "sfr_dialog.moc"

#include "common.h"
#include "config.h"
#include <locale.h>

Sfr_dialog::Sfr_dialog(QWidget* parent ATTRIBUTE_UNUSED, const Sfr_entry& entry) : cursor_domain_value(0), repainting(0) {
    
    chart = new QChart();
    chart->legend()->hide();
    
    chart->setAnimationOptions(QChart::NoAnimation);
    
    logger.info("sfr entry size: %ld\n", (*entry.info.sfr).size());

    x_axis = new QValueAxis();
    chart->addAxis(x_axis, Qt::AlignBottom);
    
    y_axis = new QValueAxis();
    chart->addAxis(y_axis, Qt::AlignLeft);
    
    entries.push_back(entry);
    update_lp_mm_mode();
    view.update(entries, series, *chart, *x_axis, *y_axis);
    
    chart_view = new Sfr_chartview(chart, this);
    chart_view->setRenderHint(QPainter::Antialiasing);
    
    
    label_layout = new QGridLayout();
    
    mtfmapper_logo = new QIcon;
    mtfmapper_logo->addFile(":/Icons/AppIcon256");
    
    QFontMetrics fm(QWidget::fontMetrics());
    int double_height = fm.height()*2 + 8;
    
    QString button_style(
        "QPushButton:flat{ background-color: rgba(0, 0, 0, 7%); border: 1px solid rgba(0,0,0,20%); border-radius: 2px; } "
        "QPushButton:flat::pressed{ background-color: rgba(0, 0, 0, 15%); border: 1px solid rgba(0,0,0,20%); border-radius: 2px; } "
    );
    
    save_img_button = new QPushButton("Save\nimage", this);
    save_img_button->setMinimumWidth(fm.width(save_img_button->text()));
    save_img_button->setMaximumWidth(fm.width(save_img_button->text()));
    save_img_button->setMinimumHeight(double_height);
    save_img_button->setMaximumHeight(double_height);
    save_img_button->setFlat(true);
    save_img_button->setStyleSheet(button_style);
    label_layout->addWidget(save_img_button, 0, 0, 3, 1);
    
    
    for (int i=0; i < 4; i++) {
        cursor_label.push_back(new QLabel(this));
        cursor_label.back()->setAlignment(i == 0 ? Qt::AlignRight : Qt::AlignLeft);
        cursor_label.back()->setStyleSheet("font-weight: bold;");
    }
    label_layout->addWidget(cursor_label[0], 0, 2);
    label_layout->addWidget(cursor_label[1], 0, 3);
    label_layout->addWidget(cursor_label[2], 1, 3);
    label_layout->addWidget(cursor_label[3], 2, 3);
    cursor_label[2]->setMinimumHeight(0);
    cursor_label[3]->setMinimumHeight(0);
    cursor_label[2]->hide();
    cursor_label[3]->hide();
    label_layout->rowMinimumHeight(0);
    label_layout->setVerticalSpacing(1);
    
    save_data_button = new QPushButton("Save\ndata", this);
    save_data_button->setMinimumWidth(fm.width(save_data_button->text()));
    save_data_button->setMaximumWidth(fm.width(save_data_button->text()));
    save_data_button->setMinimumHeight(double_height);
    save_data_button->setMaximumHeight(double_height);
    save_data_button->setFlat(true);
    save_data_button->setStyleSheet(button_style);
    label_layout->addWidget(save_data_button, 0, 1, 3, 1);
    
    // generate an alpha-blended version of the logo
    QPixmap logo_pixmap(mtfmapper_logo->pixmap(50));
    QImage logo_image(logo_pixmap.size(), QImage::Format_ARGB32_Premultiplied);
    logo_image.fill(Qt::transparent);
    QPainter p(&logo_image);
    p.setOpacity(0.5);
    p.drawPixmap(0, 0, logo_pixmap);
    p.end();
    QPixmap blended_pixmap = QPixmap::fromImage(logo_image);
    
    QLabel* logo = new QLabel;
    logo->setAlignment(Qt::AlignRight);
    logo->setPixmap(blended_pixmap);
    logo->setIndent(blended_pixmap.size().width());
    logo->setMinimumSize(QSize(blended_pixmap.size().width()*2, blended_pixmap.size().height()));
    logo->setMaximumSize(QSize(blended_pixmap.size().width()*2, blended_pixmap.size().height()));
    label_layout->addWidget(logo, 0, 5, 3, 1);
    
    box_view = new QComboBox;
    box_view->addItem("SFR");
    box_view->addItem("ESF");
    box_view->addItem("LSF");
    
    QLabel* view_label = new QLabel("Plot type:");
    QVBoxLayout* view_layout = new QVBoxLayout;
    
    view_layout->addWidget(view_label);
    view_layout->addWidget(box_view);
    QString view_tooltip("Hold down <ctrl> to display the ESF,\n or <alt> to display the LSF");
    box_view->setToolTip(view_tooltip);
    view_label->setToolTip(view_tooltip);
    
    label_layout->addLayout(view_layout, 0, 4, 3, 1);
    
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
    setWindowTitle(view.title().c_str());
    
    connect(save_img_button, SIGNAL(clicked()), this, SLOT(save_image()));
    connect(save_data_button, SIGNAL(clicked()), this, SLOT(save_data()));
    connect(box_view, SIGNAL(currentIndexChanged(int)), this, SLOT(plot_type_changed(int)));

    show();
}

void Sfr_dialog::clear(void) {
    chart->removeAllSeries();
    series.clear();
    for (auto l: cursor_label) {
        l->setText("");
        l->hide();
    }
    for (auto m: mtf50_text) {
        m->setText("");
    }
    // TODO: Not sure why we have to forcibly delete all the entries; maybe we can just flag them
    for (auto e: entries) {
        e.clear();
    }
    y_axis->setMax(0.1);
}

void Sfr_dialog::reject(void) {
    clear();
    emit sfr_dialog_closed();
    QDialog::reject();
}

void Sfr_dialog::paintEvent(QPaintEvent* event) {
    repainting.testAndSetAcquire(0, 1);
    QDialog::paintEvent(event);
    
    bool changed = update_lp_mm_mode();
    
    Qt::KeyboardModifiers current_modifier = QGuiApplication::queryKeyboardModifiers();
    if (current_modifier != last_modifier) {
        // update
        if (current_modifier & Qt::ControlModifier && !(current_modifier & Qt::AltModifier)) {
            view.push_view();
            view.set_view(Entry_view::VIEW_ESF);
        } else {
            if (current_modifier & Qt::AltModifier && !(current_modifier & Qt::ControlModifier)) {
                view.push_view();
                view.set_view(Entry_view::VIEW_LSF);
            } else {
                view.pop_view();
            }
        }
        box_view->setCurrentIndex((int)view.get_view());
        changed = true;
    }
    
    if (changed) {
        view.update(entries, series, *chart, *x_axis, *y_axis);
    }
    
    last_modifier = current_modifier;

    vector<double> contrast_list;
    if (series.size() > 0) {

        for (size_t si=0; si < series.size(); si++) {
        
            QVector<QPointF> pts = series[si]->pointsVector();
            double contrast = pts[0].y();
            for (int i = 0; i < pts.size() && pts[i].x() < cursor_domain_value; i++) {
                contrast = pts[i].y();
            }
            
            contrast_list.push_back(contrast);
        }
        

        QFontMetrics fm(QWidget::fontMetrics());
        int th = fm.height();
        int text_end = chart->size().width() - 2*fm.width("m");
        
        // add mtf50 tags in reverse
        for (int mi=series.size()-1; mi >= 0; mi--) {
            
            //string mtf50_str = view.mtf_xx(entries[mi].info.mtf_contrast, entries[mi].info.mtf50, entries[mi].info.pixel_pitch);
            string mtf50_str = view.mtf_xx(entries[mi]);
            int tw = fm.width(mtf50_str.c_str()) + 2*fm.width("m");
            
            mtf50_text[mi]->setBrush(series[mi]->pen().color());
            mtf50_text[mi]->setPen(QPen(series[mi]->pen().color()));
            mtf50_text[mi]->setPos(text_end - tw, th);
            mtf50_text[mi]->setText(mtf50_str.c_str());
            text_end -= tw;
        }
    }
    
    QPointF tpos(cursor_domain_value, y_axis->max() - 0.001);
    QPointF bpos(cursor_domain_value, y_axis->min() + 0.001);
    tpos = chart->mapToPosition(tpos);
    bpos = chart->mapToPosition(bpos);
    tpos.rx() = std::max(tpos.x()-1, chart->mapToPosition(QPointF(x_axis->min(), 0)).x());
    bpos.rx() = tpos.rx() + 2;
    
    mtf50_rect->setRect(tpos.x(), tpos.y(), bpos.x() - tpos.x(), bpos.y() - tpos.y());
    
    cursor_label[0]->setText(view.x_format(cursor_domain_value).c_str());
    cursor_label[0]->show();
    
    for (size_t mi=0; mi < contrast_list.size(); mi++) {
        cursor_label[mi+1]->setText(view.y_format(contrast_list[mi]).c_str());
        
        QPalette palette = cursor_label[mi+1]->palette();
        palette.setColor(cursor_label[mi+1]->foregroundRole(), series[mi]->pen().color());
        cursor_label[mi+1]->setPalette(palette);
        cursor_label[mi+1]->show();
    }
    
    chart->update(); // this causes some lag, but it elliminates exposed cruft. a better solution would be nice
    repainting.testAndSetRelease(1, 0);
}

void Sfr_dialog::replace_entry(const Sfr_entry& entry) {
    if (entries.size() < 3) {
        entries.push_back(entry);
    } else {
        entries.back() = entry;
    }
    view.update(entries, series, *chart, *x_axis, *y_axis);
    
    update();
    show();
    raise();
    activateWindow();
}

void Sfr_dialog::add_entry(const Sfr_entry& entry) {
    replace_entry(entry);
}

void Sfr_dialog::notify_mouse_position(double value, bool click) { // we still need this to update the bottom label, but we could merge this in the chart?
    lock_cursor ^= click;
    if (!lock_cursor) {
        cursor_domain_value = std::min(x_axis->max(), std::max(x_axis->min(), value));
        update();
    }
}

void Sfr_dialog::save_image(void) {
    QString savename = QFileDialog::getSaveFileName(
        this,
        tr("Save plot image"),
        QString(),
        QString("*.png")
    );
    
    // must append extension if none were specified
    if (!savename.contains('.')) {
        savename += ".png";
    }
    
    QImage img(size(), QImage::Format_RGB888);
    QPainter painter(&img);
    render(&painter);
    img.save(savename);
}

void Sfr_dialog::save_data(void) {
    setlocale(LC_ALL, "C");

    QString savename = QFileDialog::getSaveFileName(
        this,
        tr("Save CSV data"),
        QString(),
        QString("*.csv")
    );
    
    // must append extension if none were specified
    if (!savename.contains('.')) {
        savename += ".csv";
    }
    
    FILE* fout = fopen(savename.toLocal8Bit().data(), "wt");
    if (fout) {
        fprintf(fout, "frequency,");
        if (series.size() == 1) {
            fprintf(fout, "contrast\n");
            for (int i=0; i < series[0]->count()/20; i++) {
                fprintf(fout, "%.4lf,%.8lf\n", series[0]->at(i*20).x(), series[0]->at(i*20).y());
            }
        } else {
            size_t max_size = 0;
            for (size_t i = 0; i < series.size(); i++) {
                max_size = std::max(max_size, size_t(series[i]->count()));
            }
            max_size /= 20;
            int j;
            for (j=0; j < (int)series.size() - 1; j++) {
                fprintf(fout, "contrast_%d,", j);
            }
            fprintf(fout, "contrast_%d\n", j);
            
            for (size_t i=0; i < max_size; i++) {
                fprintf(fout, "%.4lf,", (int)i*20 < series[0]->count() ? series[0]->at(i*20).x() : 0.0);
                for (j=0; j < (int)series.size() - 1; j++) {
                    fprintf(fout, "%.8lf,", (int)i*20 < series[j]->count() ? series[j]->at(i*20).y() : 0.0);
                }
                fprintf(fout, "%.8lf\n", (int)i*20 < series[j]->count() ? series[j]->at(i*20).y() : 0.0);
            }
        }
        fclose(fout);
    } else {
        // display failure dialog
    }
}

bool Sfr_dialog::update_lp_mm_mode(void) {
    QSettings settings("mtfmapper", "mtfmapper");
    
    bool changed = false;
    changed |= view.set_lp_mm_mode((Qt::CheckState)settings.value("setting_lpmm").toInt() == Qt::Checked);
    changed |= view.set_default_pixel_pitch(settings.value("setting_pixelsize").toFloat());
    
    return changed;
}


void Sfr_dialog::plot_type_changed(int index) {
    
    switch(index) {
    case 0: view.set_view(Entry_view::VIEW_SFR); break;
    case 1: view.set_view(Entry_view::VIEW_ESF); break;
    case 2: view.set_view(Entry_view::VIEW_LSF); break;
    }
    
    view.update(entries, series, *chart, *x_axis, *y_axis);
    setWindowTitle(view.title().c_str());
}

void Sfr_dialog::keyPressEvent(QKeyEvent* /*event*/) {
    update();
}

void Sfr_dialog::keyReleaseEvent(QKeyEvent* /*event*/) {
    update();
}
