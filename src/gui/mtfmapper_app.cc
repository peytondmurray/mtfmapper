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
#include <QtGui> 
#include "mtfmapper_app.h"
#include "mtfmapper_app.moc"

#include "worker_thread.h"

#include <string>
#include <iostream>

using std::cout;
using std::endl;
using std::string; 
 
mtfmapper_app::mtfmapper_app(QWidget *parent)
  : processor(this)
{
    
    
    //QLabel* zoomLabel = new QLabel(tr("Enter a zoom value between %1 and %2:").arg(0).arg(1000));
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
    
    tb_img_annotated = new QCheckBox("Annotated");
    tb_img_annotated->setChecked(true);
    tb_img_profile = new QCheckBox("Profile");
    tb_img_profile->setChecked(true);
    tb_img_gridimg = new QCheckBox("2D map");
    tb_img_gridimg->setChecked(true);
    tb_img_gridsurf = new QCheckBox("3D map");
    tb_img_gridsurf->setChecked(true);
    
    qgs = new QGraphicsScene;
    qgs->setSceneRect(0,0,400,400);
    qgv = new QGraphicsView(qgs);
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
    
    
    QGridLayout* tb_layout = new QGridLayout;
    tb_layout->addWidget(tb_img_annotated, 1, 0);
    tb_layout->addWidget(tb_img_profile, 1, 1);
    tb_layout->addWidget(tb_img_gridimg, 2, 0);
    tb_layout->addWidget(tb_img_gridsurf, 2, 1);
    tb_layout->addWidget(datasets, 3, 0, 1, 2);
    QGroupBox* vbox2 = new QGroupBox(tr("selection"));
    vbox2->setLayout(tb_layout);

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

    
    QGroupBox* vGroupBox = new QGroupBox(tr("output"));
    QGridLayout* vlayout = new QGridLayout;
    vlayout->addWidget(qgv, 0, 0);
    vlayout->addWidget(zoom_spinbox, 1, 0);
    vlayout->addWidget(vbox2, 0, 1);
    vlayout->addWidget(v3GroupBox);
    vGroupBox->setLayout(vlayout);
    
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(vGroupBox);
    mainLayout->addWidget(progress);
    img_frame->setLayout(mainLayout);
    
    setCentralWidget(img_frame);
    
    worker_thread = new QThread;
    
    settings = new Settings_dialog(this);
    
    create_actions();
    
    file_menu = new QMenu(tr("&File"), this);
    file_menu->addAction(open_act);
    file_menu->addSeparator();
    file_menu->addAction(exit_act);
    
    settings_menu = new QMenu(tr("&Settings"), this);
    settings_menu->addAction(prefs_act);
    
    menuBar()->addMenu(file_menu);
    menuBar()->addMenu(settings_menu);
    
    connect(datasets, SIGNAL(clicked(const QModelIndex&)), this, SLOT(dataset_selected(const QModelIndex&)));
    connect(datasets->selectionModel(), SIGNAL(currentChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(dataset_selected_changed(const QModelIndex&, const QModelIndex&)));
    
    connect(&processor, SIGNAL(send_parent_item(QString, QString)), this, SLOT(parent_item(QString, QString)));
    connect(&processor, SIGNAL(send_child_item(QString, QString)), this, SLOT(child_item(QString, QString)));
    connect(&processor, SIGNAL(send_close_item()), this, SLOT(close_item()));
    connect(&processor, SIGNAL(send_delete_item(QString)), this, SLOT(item_for_deletion(QString)));
    connect(&processor, SIGNAL(send_exif_filename(QString, QString)), this, SLOT(populate_exif_info_from_file(QString, QString)));
    
    connect(tb_img_annotated, SIGNAL(clicked()), this, SLOT(img_annotated_toggled()));
    connect(tb_img_profile, SIGNAL(clicked()), this, SLOT(img_profile_toggled()));
    connect(tb_img_gridimg, SIGNAL(clicked()), this, SLOT(img_gridimg_toggled()));
    connect(tb_img_gridsurf, SIGNAL(clicked()), this, SLOT(img_gridsurf_toggled()));
    
    connect(zoom_spinbox, SIGNAL(valueChanged(int)), this, SLOT(zoom_changed(int)));
    connect(img_frame, SIGNAL(zoom_in()), this, SLOT(zoom_in()));
    connect(img_frame, SIGNAL(zoom_out()), this, SLOT(zoom_out()));
    connect(img_frame, SIGNAL(zoom_to_100()), this, SLOT(zoom_to_100()));
    
    connect(&processor, SIGNAL(send_progress_indicator(int)), progress, SLOT(setValue(int)));
    connect(settings, SIGNAL(argument_string(QString)), &processor, SLOT(receive_arg_string(QString)));
    
    setWindowTitle(tr("Image Viewer"));
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
        
        /*bool fr = */QFile().remove(fn);
        /*bool dr = */QDir().rmdir(dn);
            
        //cout << "f:" << fn.toAscii().constData() << ":" << fr << endl;
        //cout << "d:" << dn.toAscii().constData() << ":" << dr << endl;
    }
}

void mtfmapper_app::create_actions(void) {
    open_act = new QAction(tr("&Open..."), this);
    open_act->setShortcut(tr("Ctrl+O"));
    connect(open_act, SIGNAL(triggered()), this, SLOT(open()));
    
    exit_act = new QAction(tr("E&xit"), this);
    exit_act->setShortcut(tr("Ctrl+Q"));
    connect(exit_act, SIGNAL(triggered()), this, SLOT(close()));
    
    prefs_act = new QAction(tr("&Preferences"), this);
    prefs_act->setShortcut(tr("Ctrl-P"));
    connect(prefs_act, SIGNAL(triggered()), settings, SLOT( open() ));
}

void mtfmapper_app::view_image(const QString& fname) {
    QImage image(fname);
    if (image.isNull()) {
        QMessageBox::information(
            this, tr("Image Viewer"),
            tr("Cannot load %1.").arg(fname)
        );
        return;
    }
    int rwidth  = int(image.width() * (zoom_spinbox->value() / 100.0));
    int rheight = int(image.height() * (zoom_spinbox->value() / 100.0));
    qgpi->setPixmap(QPixmap::fromImage(image).scaled(QSize(rwidth,rheight), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    qgs->setSceneRect(QRectF(0,0,rwidth, rheight));
} 
 
void mtfmapper_app::open()
{
    input_files = QFileDialog::getOpenFileNames(
        this,
        "Choose files to open",
        QString::null,
        QString::null);
        
    if (input_files.size() > 0) {
    
        clear_temp_files();
        dataset_contents.clear(); // this could be optional, but issues with overwriting the temp dirs may cause problems
        dataset_files.clear();
        exif_properties.clear();
        
        progress->setRange(0, input_files.size()+1);
        progress->setValue(1);
        
        QStringList labels;
        labels.push_back(QString("Data set"));
        dataset_contents.setHorizontalHeaderLabels(labels);
        processor.set_files(input_files);
        processor.set_gnuplot_binary(settings->get_gnuplot_binary());
        processor.start();
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
    
    view_image(dataset_files.at(count_before));
    display_exif_properties(count_before);
}

void mtfmapper_app::dataset_selected_changed(const QModelIndex& i1, const QModelIndex& i2) {
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

void mtfmapper_app::img_annotated_toggled(void) {
    for (int i=0; i < dataset_contents.rowCount(); i++) {
        QStandardItem* current_dataset_item = dataset_contents.item(i);
        for (int j=0; j < current_dataset_item->rowCount(); j++) {
            QStandardItem* current_child = current_dataset_item->child(j);
            if (current_child->text().compare("annotated") == 0) {
                current_child->setEnabled(tb_img_annotated->isChecked());
            }
        }
        
    }
}

void mtfmapper_app::img_profile_toggled(void) {
    for (int i=0; i < dataset_contents.rowCount(); i++) {
        QStandardItem* current_dataset_item = dataset_contents.item(i);
        for (int j=0; j < current_dataset_item->rowCount(); j++) {
            QStandardItem* current_child = current_dataset_item->child(j);
            if (current_child->text().compare("profile") == 0) {
                current_child->setEnabled(tb_img_profile->isChecked());
            }
        }
        
    }
}

void mtfmapper_app::img_gridimg_toggled(void) {
    for (int i=0; i < dataset_contents.rowCount(); i++) {
        QStandardItem* current_dataset_item = dataset_contents.item(i);
        for (int j=0; j < current_dataset_item->rowCount(); j++) {
            QStandardItem* current_child = current_dataset_item->child(j);
            if (current_child->text().compare("grid2d") == 0) {
                current_child->setEnabled(tb_img_gridimg->isChecked());
            }
        }
        
    }
}

void mtfmapper_app::img_gridsurf_toggled(void) {
    for (int i=0; i < dataset_contents.rowCount(); i++) {
        QStandardItem* current_dataset_item = dataset_contents.item(i);
        for (int j=0; j < current_dataset_item->rowCount(); j++) {
            QStandardItem* current_child = current_dataset_item->child(j);
            if (current_child->text().compare("grid3d") == 0) {
                current_child->setEnabled(tb_img_gridsurf->isChecked());
            }
        }
        
    }
}

void mtfmapper_app::zoom_changed(int i) {
    dataset_selected(datasets->selectionModel()->currentIndex());
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

void mtfmapper_app::display_exif_properties(int index) {
    Exiv2_property* props = exif_properties.at(index);
    img_comment_value->setText(props->get_comment());
    af_ft_value->setText(props->get_af_tune());
    focus_distance_value->setText(props->get_focus_distance());
    focal_length_value->setText(props->get_focal_length());
}

void mtfmapper_app::populate_exif_info_from_file(QString s, QString tempdir) {

    Exiv2_property* props = new Exiv2_property(s, tempdir + "/exifinfo.txt");
    exif_properties.push_back(props);

    // actually, we could delete it right away ...
    item_for_deletion(tempdir + QString("/exifinfo.txt"));
}

void mtfmapper_app::check_if_helpers_exist(void) {
    bool gnuplot_exists = QFile::exists(settings->get_gnuplot_binary());
    bool exiv_exists = QFile::exists(settings->get_exiv2_binary());

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
}
