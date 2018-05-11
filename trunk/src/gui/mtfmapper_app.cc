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
#include "mtfmapper_app.h"
#include "mtfmapper_app.moc"

#include "worker_thread.h"
#include "common.h"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <string>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::string; 
 
mtfmapper_app::mtfmapper_app(QWidget *parent ATTRIBUTE_UNUSED)
  : processor(this), sfr_dialog(nullptr)
{

    zoom_spinbox = new QSpinBox;
    zoom_spinbox->setRange(0, 1000);
    zoom_spinbox->setSingleStep(10);
    zoom_spinbox->setSuffix("%");
    zoom_spinbox->setSpecialValueText(tr("Automatic"));
    zoom_spinbox->setValue(100);
    zoom_spinbox->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

    img_comment_label = new QLabel("Comment: ");
    img_comment_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    img_comment_value = new QLabel("N/A");
    img_comment_value->setAlignment(Qt::AlignLeft);
    img_comment_value->setStyleSheet("QLabel { color : SteelBlue; }");
    
    af_ft_label = new QLabel("AF Tune Value: ");
    af_ft_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    af_ft_value = new QLabel("N/A");
    af_ft_value->setAlignment(Qt::AlignLeft);
    af_ft_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    af_ft_value->setStyleSheet("QLabel { color : SteelBlue; }");

    focal_length_label = new QLabel("Focal Length: ");
    focal_length_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focal_length_value = new QLabel("N/A");
    focal_length_value->setAlignment(Qt::AlignLeft);
    focal_length_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focal_length_value->setStyleSheet("QLabel { color : SteelBlue; }");

    focus_distance_label = new QLabel("Focus Distance: ");
    focus_distance_label->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focus_distance_value = new QLabel("N/A");
    focus_distance_value->setAlignment(Qt::AlignLeft);
    focus_distance_value->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );
    focus_distance_value->setStyleSheet("QLabel { color : SteelBlue; }");

    
    img_frame = new Img_frame(this);
    
    tb_img_annotated = new QCheckBox("Annotated image");
    tb_img_annotated->setChecked(true);
    tb_img_profile = new QCheckBox("Profile");
    tb_img_profile->setChecked(true);
    tb_img_gridimg = new QCheckBox("Grid");
    tb_img_gridimg->setChecked(true);
    tb_img_lensprofile = new QCheckBox("Lens profile");
    tb_img_lensprofile->setChecked(true);
    tb_img_orientation = new QCheckBox("Chart orientation");
    tb_img_orientation->setChecked(true);
    
    qgs = new QGraphicsScene;
    qgs->setSceneRect(0,0,400,400);
    qgv = new Imgviewer(qgs, this);
    qgv->setDragMode(QGraphicsView::ScrollHandDrag);
    qgv->setRenderHints(QPainter::SmoothPixmapTransform);
    qgpi = new QGraphicsPixmapItem;
    qgpi->setTransformationMode(Qt::SmoothTransformation);
    qgs->addItem(qgpi);
    qgv->resize(600,300);
    
    QStringList labels;
    labels.push_back(QString("Data set"));
    dataset_contents.setHorizontalHeaderLabels(labels);
    
    datasets = new QTreeView;
    datasets->resize(300,400);
    datasets->move(610,0);
    datasets->setModel(&dataset_contents);
    datasets->setRootIsDecorated(true);
    
    progress = new QProgressBar;

    clear_button = new QPushButton("Clear results");
    clear_button->setEnabled(false);

    save_button = new QPushButton("Save all results");
    save_button->setEnabled(false);

    save_subset_button = new QPushButton("Save result subset");
    save_subset_button->setEnabled(false);

    
    QGridLayout* tb_layout = new QGridLayout;
    tb_layout->addWidget(datasets, 0, 0);
    tb_layout->addWidget(clear_button, 1, 0);
    tb_layout->addWidget(save_button, 2, 0);
    tb_layout->addWidget(save_subset_button, 3, 0);
    QGroupBox* vbox2 = new QGroupBox(tr("Data set control"));
    vbox2->setLayout(tb_layout);
    vbox2->setMinimumWidth(200);
    
    QGroupBox* v3GroupBox = new QGroupBox(tr("Image properties"));
    QGridLayout* hlayout = new QGridLayout;
    hlayout->addWidget(img_comment_label, 0, 0);
    hlayout->addWidget(img_comment_value, 0, 1);
    hlayout->addWidget(af_ft_label, 0, 2);
    hlayout->addWidget(af_ft_value, 0, 3);
    hlayout->addWidget(focal_length_label, 1, 0);
    hlayout->addWidget(focal_length_value, 1, 1);
    hlayout->addWidget(focus_distance_label, 1, 2);
    hlayout->addWidget(focus_distance_value, 1, 3);
    v3GroupBox->setLayout(hlayout);

    
    QGroupBox* vGroupBox = new QGroupBox(tr("Current output"));
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(qgv, 0, 0);
    vlayout->addWidget(zoom_spinbox, 1, 0);
    vlayout->addWidget(vbox2, 0, 1);
    vlayout->addWidget(v3GroupBox);
    vGroupBox->setLayout(vlayout);
    
    splitter = new QSplitter(Qt::Horizontal);
    splitter->addWidget(vGroupBox);
    splitter->addWidget(vbox2);
    splitter->setStretchFactor(0, 1);
    splitter->setCollapsible(0, false);
    splitter->setCollapsible(1, false);

    
    abort_button = new QPushButton("Abort");
    abort_button->hide();
    
    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->addWidget(splitter, 0, 0, 1, 2);
    mainLayout->addWidget(progress, 1, 0);
    mainLayout->addWidget(abort_button, 1, 1);
    img_frame->setLayout(mainLayout);
    
    setCentralWidget(img_frame);
    
    worker_thread = new QThread;
    
    settings = new Settings_dialog(this);
    about    = new About_dialog(this);
    help     = new Help_dialog(this);
    
    create_actions();
    
    file_menu = new QMenu(tr("&File"), this);
    file_menu->addAction(open_act);
    file_menu->addAction(open_roi_act);
    file_menu->addAction(open_focus_act);
    file_menu->addSeparator();
    file_menu->addAction(exit_act);
    
    settings_menu = new QMenu(tr("&Settings"), this);
    settings_menu->addAction(prefs_act);

    help_menu = new QMenu(tr("&Help"), this);
    help_menu->addAction(help_act);
    help_menu->addAction(about_act);
    
    menuBar()->addMenu(file_menu);
    menuBar()->addMenu(settings_menu);
    menuBar()->addMenu(help_menu);
    
    connect(datasets, SIGNAL(clicked(const QModelIndex&)), this, SLOT(dataset_selected(const QModelIndex&)));
    connect(datasets->selectionModel(), SIGNAL(currentChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(dataset_selected_changed(const QModelIndex&, const QModelIndex&)));
    
    connect(&processor, SIGNAL(send_parent_item(QString, QString)), this, SLOT(parent_item(QString, QString)));
    connect(&processor, SIGNAL(send_child_item(QString, QString)), this, SLOT(child_item(QString, QString)));
    connect(&processor, SIGNAL(send_close_item()), this, SLOT(close_item()));
    connect(&processor, SIGNAL(send_delete_item(QString)), this, SLOT(item_for_deletion(QString)));
    connect(&processor, SIGNAL(send_exif_filename(QString, QString)), this, SLOT(populate_exif_info_from_file(QString, QString)));
    
    connect(zoom_spinbox, SIGNAL(valueChanged(int)), this, SLOT(zoom_changed(int)));
    connect(img_frame, SIGNAL(zoom_in()), this, SLOT(zoom_in()));
    connect(img_frame, SIGNAL(zoom_out()), this, SLOT(zoom_out()));
    connect(img_frame, SIGNAL(zoom_to_100()), this, SLOT(zoom_to_100()));
    
    connect(&processor, SIGNAL(send_progress_indicator(int)), progress, SLOT(setValue(int)));
    connect(settings, SIGNAL(argument_string(QString)), &processor, SLOT(receive_arg_string(QString)));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(hide_abort_button()));
    connect(abort_button, SIGNAL(clicked()), &processor, SLOT(receive_abort()));
    
    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_clear_button()));
    connect(clear_button, SIGNAL(clicked()), this, SLOT(clear_button_pressed()));
    connect(clear_button, SIGNAL(clicked()), this, SLOT(disable_save_button()));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_save_button()));
    connect(save_button, SIGNAL(clicked()), this, SLOT(save_button_pressed()));

    connect(save_subset_button, SIGNAL(clicked()), this, SLOT(save_subset_button_pressed()));

    connect(&processor, SIGNAL(send_all_done()), this, SLOT(enable_file_open()));

    mtfmapper_logo = new QIcon;
    mtfmapper_logo->addFile(":/Icons/AppIcon256");
    qgpi->setPixmap(mtfmapper_logo->pixmap(256));
    qgs->setSceneRect(QRectF(0, 0, 255, 255));


    setWindowTitle(tr("MTF Mapper"));
    resize(920,600);
    
    settings->send_argument_string();
    check_if_helpers_exist();
}

mtfmapper_app::~mtfmapper_app(void) {
    // clean up all the old folders
    clear_temp_files();
}

void mtfmapper_app::clear_temp_files(void) {
    QStringList dirnames;
    for (int i=0; i < tempfiles_to_delete.size(); i++) {
        QString fn(tempfiles_to_delete.at(i));
        QString dn(QFileInfo(fn).absolutePath());
        
        QFile().remove(fn);
        QDir().rmdir(dn);
            
    }
    clear_button->setEnabled(false);
}

void mtfmapper_app::create_actions(void) {
    open_act = new QAction(tr("&Open..."), this);
    open_act->setShortcut(tr("Ctrl+O"));
    connect(open_act, SIGNAL(triggered()), this, SLOT(open_auto()));
    
    open_roi_act = new QAction(tr("&Open single edge image(s)..."), this);
    open_roi_act->setShortcut(tr("Ctrl+R"));
    connect(open_roi_act, SIGNAL(triggered()), this, SLOT(open_roi()));
    
    open_focus_act = new QAction(tr("&Open Focus Position image(s)..."), this);
    open_focus_act->setShortcut(tr("Ctrl+F"));
    connect(open_focus_act, SIGNAL(triggered()), this, SLOT(open_focus()));
    
    exit_act = new QAction(tr("E&xit"), this);
    exit_act->setShortcut(tr("Ctrl+Q"));
    connect(exit_act, SIGNAL(triggered()), this, SLOT(close()));
    
    prefs_act = new QAction(tr("&Preferences"), this);
    prefs_act->setShortcut(tr("Ctrl-P"));
    connect(prefs_act, SIGNAL(triggered()), settings, SLOT( open() ));

    help_act = new QAction(tr("&Help"), this);
    help_act->setShortcut(tr("Ctrl-H"));
    connect(help_act, SIGNAL(triggered()), help, SLOT( open() ));

    about_act = new QAction(tr("&About"), this);
    connect(about_act, SIGNAL(triggered()), about, SLOT( open() ));
}

void mtfmapper_app::view_image(const QString& fname) {
    QImage image(fname);
    if (image.isNull()) {
        // Qt could not load the file. This might be because it is a 16-bit TIFF file.
        // Try loading it using opencv
        cv::Mat cvimg;
        bool opencv_succeeded = true;
        try {
            cvimg = cv::imread(fname.toStdString(), CV_LOAD_IMAGE_GRAYSCALE);
            if (cvimg.data) {
                // resize the image if it exceeds 10k pixels, since this seems to break QT 5.7 ?
                if (cvimg.cols > 9999 || cvimg.rows > 9999) {
                    cv::Mat smaller;
                    double sf = 1;
                    if (cvimg.cols >= cvimg.rows) {
                        sf = 8000.0 / double(cvimg.cols);
                    } else {
                        sf = 80000.0 / double(cvimg.rows);
                    }
                    cv::resize(cvimg, smaller, cv::Size(0, 0), sf, sf);
                    cvimg = smaller;
                    logger.debug("Input image %s resized to %d x %d\n", fname.toLocal8Bit().constData(), cvimg.cols, cvimg.rows);
                }
                image = QImage(cvimg.cols, cvimg.rows, QImage::Format_Grayscale8);
                memcpy(image.bits(), cvimg.data, cvimg.rows*cvimg.cols);
                logger.debug("Input image %s loaded through OpenCV\n", fname.toLocal8Bit().constData());
            } else {
                // too bad, OpenCV cannot handle it either ...
                opencv_succeeded = false;
            }
        }
        catch (const cv::Exception&) {
            // too bad, OpenCV cannot handle it either ...
            opencv_succeeded = false;
        }
        if (!opencv_succeeded) {
            QMessageBox::information(
                this, tr("Image Viewer"),
                tr("Cannot load %1.").arg(fname)
            );
            return;
        }
        
    }
    int rwidth, rheight;
    
    if (zoom_spinbox->value() == 0) {
        QSize vp_size = qgv->maximumViewportSize();
        rwidth  = int(vp_size.rwidth());
        rheight = int(vp_size.rwidth()*image.height()/image.width());
    } else {
        rwidth  = int(image.width() * (zoom_spinbox->value() / 100.0));
        rheight = int(image.height() * (zoom_spinbox->value() / 100.0));
    }
    qgpi->setPixmap(QPixmap::fromImage(image).scaled(QSize(rwidth,rheight), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    qgs->setSceneRect(QRectF(0,0,rwidth, rheight));
} 

void mtfmapper_app::open_auto() {
    open_action(false);
}

void mtfmapper_app::open_roi() {
    open_action(true, false);
}

void mtfmapper_app::open_focus() {
    open_action(false, true);
}
 
void mtfmapper_app::open_action(bool roi, bool focus) {

    QFileDialog* open_dialog = new QFileDialog(this, tr("Select input files"), QString::null, QString::null);
    open_dialog->setOption(QFileDialog::DontUseNativeDialog);
    
    if (!(roi || focus)) {
        QGroupBox* v4GroupBox = new QGroupBox(tr("Select desired MTF Mapper outputs to produce:"));
        QGridLayout* ft_gridbox = new QGridLayout();
        if (ft_gridbox) {
            ft_gridbox->addWidget(tb_img_annotated, 0, 0);
            ft_gridbox->addWidget(tb_img_profile, 0, 1);
            ft_gridbox->addWidget(tb_img_gridimg, 0, 2);
            ft_gridbox->addWidget(tb_img_lensprofile, 1, 0);
            ft_gridbox->addWidget(tb_img_orientation, 1, 1);
        }
        v4GroupBox->setLayout(ft_gridbox);

        QGridLayout* od_gridbox = qobject_cast<QGridLayout*>(open_dialog->layout());
        od_gridbox->addWidget(v4GroupBox);
    }
    
    open_dialog->setFileMode(QFileDialog::FileMode::ExistingFiles);
    
    // use the state from the settings menu as a starting point
    tb_img_annotated->setCheckState(settings->cb_annotation->checkState());
    tb_img_profile->setCheckState(settings->cb_profile->checkState());
    tb_img_gridimg->setCheckState(settings->cb_grid->checkState());
    tb_img_lensprofile->setCheckState(settings->cb_lensprofile->checkState());
    tb_img_orientation->setCheckState(settings->cb_orientation->checkState());

    if (open_dialog->exec()) {
        // write state back to settings menu
        settings->cb_annotation->setCheckState(tb_img_annotated->checkState());
        settings->cb_profile->setCheckState(tb_img_profile->checkState());
        settings->cb_grid->setCheckState(tb_img_gridimg->checkState());
        settings->cb_lensprofile->setCheckState(tb_img_lensprofile->checkState());
        settings->cb_orientation->setCheckState(tb_img_orientation->checkState());
        settings->set_gnuplot_img_width(int(qgv->size().height()*1.3));
        settings->send_argument_string();

        input_files = open_dialog->selectedFiles();
        if (input_files.size() > 0) {

            open_act->setEnabled(false);
            exit_act->setEnabled(false);

            progress->setRange(0, input_files.size() + 1);

            QStringList labels;
            labels.push_back(QString("Data set"));
            dataset_contents.setHorizontalHeaderLabels(labels);
            abort_button->show();
            processor.set_single_roi_mode(roi);
            processor.set_focus_mode(focus);
            processor.set_files(input_files);
            processor.set_gnuplot_binary(settings->get_gnuplot_binary());
            processor.set_dcraw_binary(settings->get_dcraw_binary());
            processor.set_exiv2_binary(settings->get_exiv2_binary());
            processor.start();
        }
    }
    
}
 
void mtfmapper_app::dataset_selected(const QModelIndex& index) {
    // Horrible hack to determine the index of the actual 
    // filename associated with this entry. There must be a 
    // better way ...
    int count_before = 0;
    int parent_row = 0;
    if (index.parent() != QModelIndex()) {
        parent_row = index.parent().row();
        count_before = index.row() + 1;
    } else {
        parent_row = index.row();
    }
    
    for (int row=parent_row-1; row >= 0; row--) {
        QStandardItem* current_dataset_item = dataset_contents.item(row);
        count_before += current_dataset_item->rowCount() + 1;
    }
    if (dataset_contents.itemFromIndex(index)->isEnabled()) {
        view_image(dataset_files.at(count_before));
        display_exif_properties(count_before);
        if (dataset_contents.itemFromIndex(index)->text().compare(QString("annotated")) == 0) {
            qgv->set_clickable(true);
            sfr_list.clear();
            
            QString sfr_source = QFileInfo(dataset_files.at(count_before)).dir().path() + QString("/edge_sfr_values.txt");
            // go and fetch the corresponding "edge_sfr" entries
            std::ifstream ifile(sfr_source.toLocal8Bit().data());
            string line;
            while (!ifile.fail()) {
                std::getline(ifile, line);
                std::istringstream is(line);
                vector<double> values;
                std::copy(std::istream_iterator<double>(is), std::istream_iterator<double>(), std::back_inserter(values));
                if (values.size() > 6 && values[5] > 0.0) {
                    vector<double> sfr;
                    std::copy(values.begin() + 5, values.end(), std::back_inserter(sfr));
                    sfr_list.push_back(Sfr_entry(values[1], values[2], sfr));
                }
            }
        } else {
            qgv->set_clickable(false);
            sfr_list.clear();
        }
    }
}

void mtfmapper_app::dataset_selected_changed(const QModelIndex& i1, const QModelIndex& i2 ATTRIBUTE_UNUSED) {
    dataset_selected(i1);
}
 
void mtfmapper_app::parent_item(QString s, QString f) {
    current_dataset_item = new QStandardItem(s);
    dataset_files.push_back(f);
}

void mtfmapper_app::child_item(QString s, QString f) {
    QStandardItem* child = new QStandardItem(s);
    child->setEditable(false);
    current_dataset_item->appendRow(child);
    dataset_files.push_back(f);
    exif_properties.push_back(exif_properties.back());
}

void mtfmapper_app::close_item(void) {
    dataset_contents.appendRow(current_dataset_item);
    datasets->setModel(&dataset_contents);
}

void mtfmapper_app::item_for_deletion(QString s) {
    tempfiles_to_delete.push_back(s);
}


void mtfmapper_app::zoom_changed(int i ATTRIBUTE_UNUSED) {
    if (datasets->selectionModel()->currentIndex().isValid()) {
        dataset_selected(datasets->selectionModel()->currentIndex());
    }
}

void mtfmapper_app::zoom_in(void) {
    zoom_spinbox->setValue(zoom_spinbox->value() + 10);
}

void mtfmapper_app::zoom_out(void) {
    zoom_spinbox->setValue(zoom_spinbox->value() - 10);
}

void mtfmapper_app::zoom_to_100(void) {
    zoom_spinbox->setValue(100);
}

void mtfmapper_app::hide_abort_button(void) {
    abort_button->hide();
}

void mtfmapper_app::enable_clear_button(void) {
    clear_button->setEnabled(true);
}

void mtfmapper_app::clear_button_pressed(void) {
    clear_temp_files();
    dataset_contents.clear(); 
    dataset_files.clear();
    exif_properties.clear();

    img_comment_value->setText("N/A");
    af_ft_value->setText("N/A");
    focal_length_value->setText("N/A");
    focus_distance_value->setText("N/A");

    qgpi->setPixmap(mtfmapper_logo->pixmap(256));
    qgs->setSceneRect(QRectF(0, 0, 255, 255));

    clear_button->setEnabled(false);
}

void mtfmapper_app::enable_save_button(void) {
    save_button->setEnabled(true);
    save_subset_button->setEnabled(true);
}

void mtfmapper_app::disable_save_button(void) {
    save_button->setEnabled(false);
    save_subset_button->setEnabled(false);
}


void mtfmapper_app::enable_file_open(void) {
    open_act->setEnabled(true);
    exit_act->setEnabled(true);
}

void mtfmapper_app::save_button_pressed(void) {
    save_action(false);
}

void mtfmapper_app::save_subset_button_pressed(void) {
    save_action(true);
}

void mtfmapper_app::save_action(bool subset) {
    // lock file->open to prevent messing with the dataset list
    bool open_was_enabled = open_act->isEnabled();
    open_act->setEnabled(false);
    exit_act->setEnabled(false);

    std::map<std::string, int> keepers;
    bool cancelled = false;
    if (subset) {
        QListWidget* subset_list = new QListWidget();
        for (int i = 0; i < dataset_contents.rowCount(); i++) {
            QStandardItem* current_dataset_item = dataset_contents.item(i);
            subset_list->addItem(current_dataset_item->text());
        }
        QDialog* subset_mb = new QDialog(this);
        subset_mb->setWindowTitle("Save result subset");
        QGridLayout* od_gridbox = new QGridLayout();
        od_gridbox->addWidget(new QLabel("Select result subsets to save:"), 0, 0);
        QSize vp_size = qgv->maximumViewportSize();
        int rowheight = subset_list->model()->rowCount() * subset_list->sizeHintForRow(0);
        rowheight = std::min(rowheight, vp_size.rheight());
        subset_list->setMaximumHeight(rowheight + subset_list->frameWidth()*2);
        subset_list->setMinimumHeight(rowheight + subset_list->frameWidth()*2);
        subset_list->setSelectionMode(QAbstractItemView::ExtendedSelection);
        od_gridbox->addWidget(subset_list, 1, 0, 1, 2);
        QPushButton* savebutton = new QPushButton("Save", subset_mb);
        QPushButton* cancelbutton = new QPushButton("Cancel", subset_mb);
        od_gridbox->addWidget(savebutton, 2, 0);
        od_gridbox->addWidget(cancelbutton, 2, 1);
        subset_mb->setLayout(od_gridbox);
        subset_mb->setModal(true);
        
        connect(cancelbutton, &QPushButton::clicked, [&cancelled, &subset_mb](){ cancelled = true;  subset_mb->close(); });
        connect(savebutton, &QPushButton::clicked, subset_mb, &QDialog::close);
        subset_mb->exec();

        if (!cancelled) {
            QList<QListWidgetItem *> selection = subset_list->selectedItems();
            for (int i = 0; i < selection.size(); i++) {
                keepers[selection[i]->text().toStdString()] = 1;
            }
        }
        delete subset_mb;
    }
    if (!cancelled) {
        QString save_path = QFileDialog::getExistingDirectory(
            this,
            tr("Choose directory to save results in"),
            QDir::homePath()
        );
        int overwrite_count = 0;
        int idx = 0;
        std::vector< std::pair<QString, QString> > copy_list;
        for (int i = 0; i < dataset_contents.rowCount(); i++) {
            QStandardItem* current_dataset_item = dataset_contents.item(i);
            QString prefix = current_dataset_item->text();
            idx++; // skip the actual image file
            for (int j = 0; j < current_dataset_item->rowCount(); j++) {
                QStandardItem* current_child = current_dataset_item->child(j);
                QString dest_fname = save_path + "/" + prefix;
                QString src_fname = dataset_files.at(idx++);
                dest_fname += "_" + current_child->text() + ".png";
                bool copy_allowed = false;
                if (subset) {
                    if (keepers.find(prefix.toStdString()) != keepers.end()) {
                        copy_allowed = true;
                    }
                }
                else {
                    copy_allowed = true;
                }
                if (copy_allowed) {
                    copy_list.push_back(std::make_pair(src_fname, dest_fname));
                    logger.info("cp [%s] [%s]\n", src_fname.toLocal8Bit().constData(), dest_fname.toLocal8Bit().constData());
                    if (QFile::exists(dest_fname)) {
                        overwrite_count++;
                    }
                }
            }
        }


        bool overwrite_ok = true;
        if (overwrite_count > 0) {
            QMessageBox::StandardButton response = QMessageBox::question(
                this,
                QString("Saving results"),
                overwrite_count == 1 ?
                tr("One of the output files already exists. Do you want to overwrite it?") :
                tr("Some output files already appear to exist. Do you want to overwrite %1 files?").arg(overwrite_count),
                QMessageBox::Yes | QMessageBox::No,
                QMessageBox::Yes
                );
            if (response == QMessageBox::No) {
                overwrite_ok = false;
            }
        }
        for (auto& cp : copy_list) {
            if (QFile::exists(cp.second)) {
                if (overwrite_ok) {
                    QFile::remove(cp.second);
                    QFile::copy(cp.first, cp.second);
                }
            }
            else {
                QFile::copy(cp.first, cp.second);
            }
        }
    }

    if (open_was_enabled) {
        open_act->setEnabled(true);
        exit_act->setEnabled(true);
    } // otherwise the worker thread will have to re-enable the open action
}

void mtfmapper_app::display_exif_properties(int index) {
    Exiv2_property* props = exif_properties.at(index);
    img_comment_value->setText(props->get_comment());
    af_ft_value->setText(props->get_af_tune());
    focus_distance_value->setText(props->get_focus_distance());

    focal_length_value->setText(props->get_focal_length() + " / " + props->get_aperture());
}

void mtfmapper_app::populate_exif_info_from_file(QString s, QString tempdir) {

    Exiv2_property* props = new Exiv2_property(settings->get_exiv2_binary(), s, tempdir + "/exifinfo.txt");
    exif_properties.push_back(props);

    // actually, we could delete it right away ...
    item_for_deletion(tempdir + QString("/exifinfo.txt"));
}

void mtfmapper_app::check_if_helpers_exist(void) {
    bool gnuplot_exists = QFile::exists(settings->get_gnuplot_binary());
    bool exiv_exists = QFile::exists(settings->get_exiv2_binary());
    bool dcraw_exists = QFile::exists(settings->get_dcraw_binary());

    if (!gnuplot_exists) {
        QMessageBox::warning(
            this, 
            QString("gnuplot helper"), 
            QString("gnuplot helper executable not found. Please configure this in the settings.")
        );
    }

    if (!exiv_exists) {
        QMessageBox::warning(
            this, 
            QString("Exiv2 helper"), 
            QString("Exiv2 helper executable not found. Please configure this in the settings.")
        );
    }

    if (!dcraw_exists) {
        QMessageBox::warning(
            this, 
            QString("dcraw helper"), 
            QString("dcraw helper executable not found. Please configure this in the settings.")
        );
    }
}

void mtfmapper_app::edge_selected(int px, int py, bool ctrl_down, bool shift_down) {
    px /= zoom_spinbox->value()/100.0;
    py /= zoom_spinbox->value()/100.0;
    if (sfr_list.size() > 0) {
        
        size_t close_idx = 0;
        double close_dist = 1e50;
        for (size_t i=0; i < sfr_list.size(); i++) {
            double d = sfr_list[i].distance(px, py);
            if (d < close_dist) {
                close_dist = d;
                close_idx = i;
            }
        }
        
        if (close_dist < 50) {
            // TODO: some form of visual confirmation of a click would be nice
            // we can probably just draw on the qgraphicsview of imgviewer?
    
            if (!sfr_dialog) {
                sfr_dialog = new Sfr_dialog(this, sfr_list[close_idx]);
            } else {
                if (shift_down) {
                    sfr_dialog->add_entry(sfr_list[close_idx]);
                } else {
                    delete sfr_dialog;
                    sfr_dialog = new Sfr_dialog(this, sfr_list[close_idx]);
                    //sfr_dialog->clear();
                    //sfr_dialog->replace_entry(sfr_list[close_idx]);
                }
            }
            
        }
    }
}

void mtfmapper_app::closeEvent(QCloseEvent* event) {
    if (sfr_dialog && sfr_dialog->isVisible()) {
        sfr_dialog->close();
    } 
    if (about && about->isVisible()) {
        about->close();
    }
    if (help && help->isVisible()) {
        help->close();
    }
    QMainWindow::closeEvent(event);
}

