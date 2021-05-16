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
#include "include/logger.h"
#include "edge_select_dialog.h"

Edge_select_dialog::Edge_select_dialog(QWidget* parent)
 : QDialog(parent, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint), parent(parent) {

    cancel_button = new QPushButton("Cancel");
    cancel_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
    accept_button = new QPushButton("Accept");
    accept_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    
    apply_all_button = new QPushButton("Accept && batch queue");
    apply_all_button->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    help_button = new QPushButton("Help");
    help_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    load_button = new QPushButton("Load");
    load_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    
    save_button = new QPushButton("Save");
    save_button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    save_button->setEnabled(false);

    img_viewer = new GL_image_viewer(parent);
    img_viewer->set_clickable(true);
    img_viewer->set_resize_on_load(true);
    
    img_panel = new GL_image_panel_edges(img_viewer);
    img_panel->setMouseTracking(true);
    img_panel->set_max_scale_factor(4.0); // allow more zoom in this dialog
    
    img_viewer->setViewport(img_panel);   // TODO: could combine these
    img_viewer->set_GL_widget(img_panel);
    
    icon_image = std::shared_ptr<QImage>(new QImage(QSize(256, 256), QImage::Format_RGB888));
    icon_image->fill(Qt::red);
    img_panel->set_default_image(icon_image.get());

    help_dialog = new Manual_roi_help_dialog(this);
    help_dialog->setModal(true);
    
    img_filename = new QLabel("filename");
    img_type = new QLabel("type");
    img_progress = new QLabel("n remaining");
    edge_length = new QLabel("");
    text_edge_length = new QLabel("Edge length:");
    text_img_filename = new QLabel("Name:");
    text_img_type = new QLabel("Type:");
    text_img_progress = new QLabel("Queued:");
    
    histogram = new Histo_widget(this);
    
    QGridLayout* prop_layout = new QGridLayout;
    prop_layout->addWidget(text_img_filename, 0, 0);
    prop_layout->addWidget(img_filename, 0, 1);
    prop_layout->addWidget(text_img_type, 1, 0);
    prop_layout->addWidget(img_type, 1, 1);
    prop_layout->addWidget(text_img_progress, 2, 0);
    prop_layout->addWidget(img_progress, 2, 1);
    
    QGroupBox* prop_box = new QGroupBox("File properties");
    prop_box->setLayout(prop_layout);
    
    QGridLayout* histo_layout = new QGridLayout;
    histo_layout->addWidget(histogram, 0, 0, 1, 4);
    histo_layout->addWidget(text_edge_length, 1, 0);
    histo_layout->addWidget(edge_length, 1, 1, Qt::AlignLeft);
    
    QGroupBox* histo_box = new QGroupBox("Split edge histogram");
    histo_box->setLayout(histo_layout);
    
    QGridLayout* button_layout = new QGridLayout;
    button_layout->addWidget(help_button, 0, 1);
    button_layout->addWidget(apply_all_button, 1, 0, 1, 2, Qt::AlignRight);
    button_layout->addWidget(load_button, 2, 0);
    button_layout->addWidget(save_button, 2, 1);
    button_layout->addWidget(accept_button, 3, 0);
    button_layout->addWidget(cancel_button, 3, 1);
    
    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(histo_box);
    hlayout->addWidget(prop_box);
    hlayout->addStretch();
    hlayout->addLayout(button_layout);
    
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(img_viewer, 0, 0);
    vlayout->addLayout(hlayout, 1, 0);
    
    connect(cancel_button, SIGNAL(clicked()), this, SLOT( close() ));
    connect(accept_button, SIGNAL(clicked()), this, SLOT( export_roi() ));
    connect(help_button, SIGNAL(clicked()), this, SLOT( show_help() ));
    connect(load_button, SIGNAL(clicked()), this, SLOT( load_roi() ));
    connect(save_button, SIGNAL(clicked()), this, SLOT( save_roi() ));
    connect(img_panel, SIGNAL( update_edge_length(double) ), this, SLOT( update_edge_length(double) ));
    connect(img_panel, SIGNAL( update_histogram(histo_t, histo_t) ), this, SLOT( update_histogram(histo_t, histo_t) ));
    connect(img_panel, SIGNAL( enable_save_button() ), this, SLOT( enable_save_button() ));
    connect(img_panel, SIGNAL( disable_save_button() ), this, SLOT( disable_save_button() ));
    connect(apply_all_button, SIGNAL( clicked() ), this, SLOT( apply_all_button_clicked() ));
    
    // TODO: save button is disabled until at least one ROI is selected
    
    setModal(true);
    setLayout(vlayout);
    setWindowTitle("Select one or more edge ROIs");
}

bool Edge_select_dialog::load_image(QString img_name, 
    QStringList arguments, QString raw_filename) {
    
    bool load_success = img_viewer->load_image(img_name);
    
    if (!load_success) {
        return false;
    }
    
    // reset Bayer parameters (for histogram extraction later)
    bayer_channel = Bayer::bayer_t::NONE;
    cfa_pattern = Bayer::cfa_pattern_t::RGGB;
    cfa_mask = Bayer::to_cfa_mask(bayer_channel, cfa_pattern);
    
    QString type_string;
    if (get_panel()->image_channels() == 1) {
        if (arguments.contains("--bayer")) {
            type_string = "Bayer ";
            extract_bayer_info(arguments);
        } else {
            type_string = "Gray ";
        }
    } else {
        type_string = "RGB ";
    }
    if (get_panel()->image_channels() == 8) {
        type_string += "8-bit";
    } else {
        type_string += "16-bit";
    }
    img_type->setText(type_string);
    img_filename->setText(QFileInfo(raw_filename).fileName());
    
    img_panel->set_cfa_mask(cfa_mask);
    img_panel->clear_overlay();
    update();
    // force a repaint on the GL panel once the event loop runs
    QTimer::singleShot(0, [this]() { 
        img_panel->update();
        QCoreApplication::processEvents(); 
    });
    return true;
}

void Edge_select_dialog::open(void) {
    show();
}

void Edge_select_dialog::export_roi(void) {
    bool roi_success = img_panel->save_rois(roi_file);
    if (roi_success) {
        accept();
    } else {
        close();
    }
}

void Edge_select_dialog::load_roi(void) {
    QFileDialog dialog(this, "Select ROI file", QString(), QString());
    dialog.setFileMode(QFileDialog::ExistingFile);
    dialog.setNameFilter(tr("ROI files (*.roi)"));
    
    QStringList filenames;
    if (dialog.exec()) {
        filenames = dialog.selectedFiles();
        if (filenames.size() > 0 && filenames[0].length() > 0) {
            bool load_success = img_panel->load_rois(filenames[0]);
        }
    }
}

void Edge_select_dialog::save_roi(void) {
    QFileDialog dialog(this, "Save ROI file", QString(), QString("*.roi"));
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setNameFilter(tr("ROI files (*.roi)"));
    
    QStringList filenames;
    if (dialog.exec()) {
        filenames = dialog.selectedFiles();
        if (filenames.size() > 0 && filenames[0].length() > 0) {
            if (QFileInfo(filenames[0]).suffix().length() == 0) {
                filenames[0] += ".roi";
            }
            img_panel->save_rois(filenames[0]);
        }
    }
}

void Edge_select_dialog::show_help(void) {
    help_dialog->setModal(false);
    help_dialog->show();
}

void Edge_select_dialog::update_queue_size(int queue_size) {
    if (queue_size > 0) {
        img_progress->setText(QString("%1 files").arg(queue_size+1));
    } else {
        img_progress->setText(QString("%1 file").arg(queue_size+1));
    }
}

void Edge_select_dialog::update_edge_length(double length) {
    if (length > 0) {
        edge_length->setText(QString("%1").arg(length, 0, 'f', 1));
    } else {
        QTimer::singleShot(200, [this]() { 
            edge_length->setText(QString(""));
        });
        
    }
}

void Edge_select_dialog::update_histogram(histo_t dark, histo_t light) {
    histogram->set_histogram(dark, light);
}

void Edge_select_dialog::extract_bayer_info(const QStringList& arguments) {
    int bayer_idx = arguments.indexOf("--bayer");
    if (bayer_idx >= 0) {
        std::string bayer_channel_str = arguments[std::min(bayer_idx + 1, arguments.size())].toStdString();
        bayer_channel = Bayer::from_string(bayer_channel_str);
        cfa_pattern = Bayer::cfa_pattern_t::RGGB;
        int cfa_idx = arguments.indexOf("--cfa-pattern");
        if (cfa_idx >= 0) {
            std::string cfa_pattern_str = arguments[std::min(cfa_idx + 1, arguments.size())].toStdString();
            cfa_pattern = Bayer::from_cfa_string(cfa_pattern_str);
        }
        cfa_mask = Bayer::to_cfa_mask(bayer_channel, cfa_pattern);
    }
}

void Edge_select_dialog::disable_save_button(void) {
    save_button->setEnabled(false);
}

void Edge_select_dialog::enable_save_button(void) {
    save_button->setEnabled(true);
}

void Edge_select_dialog::apply_all_button_clicked(void) {
    static int manual_seq_number = 0;
    QString manual_roi_file = QString("%1/mtfmapper_batch_%2.roi").arg(QDir::tempPath()).arg(manual_seq_number++);
    bool save_success = img_panel->save_rois(manual_roi_file);
    if (save_success) {
        printf("emitting manual queue start\n");
        emit start_manual_queue_processing(manual_roi_file);
        accept();
    } else {
        logger.error("Could not save batch ROI filename %s\n", 
            manual_roi_file.toLocal8Bit().constData()
        );
    }
}


