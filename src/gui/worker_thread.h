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
#ifndef WORKER_THREAD_H
#define WORKER_THREAD_H

#include <QThread>
#include <QStringList>
#include <QStandardItemModel>
#include <QDir>
#include <QRunnable>

class mtfmapper_app;

class Worker_thread : public QThread 
{
  Q_OBJECT
  
  public:
    Worker_thread(QWidget *parent);
    void set_files(const QStringList& input_files);
    void run();
    const QStandardItemModel& get_output_files(void) {
        return output_files;
    }

    void set_gnuplot_binary(const QString& s) {
        gnuplot_binary = s;
    }

    void set_exiv2_binary(const QString& s) {
        exiv2_binary = s;
    }

    void set_dcraw_binary(const QString& s) {
        dcraw_binary = s;
    }
    
  signals:
    void send_parent_item(QString s, QString f);
    void send_child_item(QString s, QString f);
    void send_close_item(void); 
    void send_delete_item(QString s); 
    void send_exif_filename(QString s, QString tempdir);
    
    void send_progress_indicator(int p);
    void send_all_done(void);
    
  public slots:
    void receive_arg_string(QString s);
    void receive_abort(void);
    
  private:
    mtfmapper_app* parent;
    QStringList input_files;
    QStandardItemModel output_files;
    QString      arguments;
    QString      gnuplot_binary;
    QString      exiv2_binary;
    QString      dcraw_binary;

    bool abort;
};

#endif

