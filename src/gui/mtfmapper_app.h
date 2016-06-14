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
#ifndef MTFMAPPER_APP_H
#define MTFMAPPER_APP_H
 
#include <QMainWindow> 
#include <QModelIndex>
#include <QStringList>
#include <QList>
#include "settings_dialog.h"
#include "worker_thread.h"
#include "img_frame.h"
#include "exiv2_property.h"
#include "about_dialog.h"
#include "help_dialog.h"
#include "imgviewer.h"

class QPushButton;
class QLabel;
class QGroupBox;
class QGraphicsScene;
class QGraphicsView;
class QStringList;
class QSpinBox;
class QGraphicsPixmapItem;
class QTreeView;
class QThread;
class QCheckBox;
class QProgressBar;
   
class mtfmapper_app : public QMainWindow
{
  Q_OBJECT
        
  public:
    mtfmapper_app(QWidget *parent = 0);
    virtual ~mtfmapper_app(void);

  private:
    void create_actions(void);
    void view_image(const QString& fname);
    void display_exif_properties(int index);
    void clear_temp_files(void);
    void check_if_helpers_exist(void);
    
    QMenu*          file_menu;
    QMenu*          settings_menu;
    QMenu*          help_menu;
    QAction*        open_act;
    QAction*        exit_act;
    QAction*        prefs_act;
    QAction*        about_act;
    QAction*        help_act;
    
    QTreeView*      datasets;
    
    QSpinBox*       zoom_spinbox;

    QLabel*         img_comment_label;
    QLabel*         img_comment_value;

    QLabel*         af_ft_label;
    QLabel*         af_ft_value;

    QLabel*         focal_length_label;
    QLabel*         focal_length_value;

    QLabel*         focus_distance_label;
    QLabel*         focus_distance_value;
    
    QGroupBox*      horizgroup;
    
    QGraphicsView*       qgv;
    QGraphicsScene*      qgs;
    QGraphicsPixmapItem* qgpi;
    
    Settings_dialog*    settings;
    About_dialog*       about;
    Help_dialog*        help;
    
    QStringList     input_files;
    QStringList     dataset_files;
    QStringList     tempfiles_to_delete;
    QList<Exiv2_property*>  exif_properties;
    
    QThread*        worker_thread;
    Worker_thread   processor;
    
    QStandardItemModel   dataset_contents;
    
    QStandardItem*  current_dataset_item;
    
    QCheckBox*      tb_img_annotated;
    QCheckBox*      tb_img_profile;
    QCheckBox*      tb_img_gridimg;
    QCheckBox*      tb_img_gridsurf;
    QCheckBox*      tb_img_focus;
    
    QProgressBar*   progress;
    QPushButton*    abort_button;
    
    Img_frame*      img_frame;
              
  public slots:
    void open();
    void dataset_selected(const QModelIndex&);
    void dataset_selected_changed(const QModelIndex&, const QModelIndex&);
    void parent_item(QString s, QString f);
    void child_item(QString s, QString f);
    void close_item(void);
    void item_for_deletion(QString s);
    void populate_exif_info_from_file(QString s, QString tempdir);
    void zoom_changed(int i);
    
    void img_annotated_toggled(void);
    void img_profile_toggled(void);
    void img_gridimg_toggled(void);
    void img_gridsurf_toggled(void);
    void img_focus_toggled(void);
  
    void zoom_in(void);
    void zoom_out(void);  
    void zoom_to_100(void);

    void hide_abort_button(void);
};
                               
                                
#endif
