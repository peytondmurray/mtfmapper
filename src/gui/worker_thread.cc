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
#include "worker_thread.h"

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include <QFileInfo>
#include <QSharedPointer>
#include <QCoreApplication>
#include <QProcess>
#include "mtfmapper_app.h"

using std::cout;
using std::endl;
using std::string;

Worker_thread::Worker_thread(QWidget* parent) 
: parent(dynamic_cast<mtfmapper_app*>(parent)), tempdir_number(0), abort(false) {

}

void Worker_thread::set_files(const QStringList& files) {
    input_files = files;
}

void Worker_thread::run(void) {
    vector<std::pair<failure_t, QString>> failures;
    abort = false;
    output_files.clear();
    QString arguments = update_arguments(settings_arguments);
    for (int i=0; i < input_files.size() && !abort; i++) {
        emit send_progress_indicator(i+1);
        QString tempdir = tr("%1/mtfmappertemp_%2").arg(QDir::tempPath()).arg(tempdir_number++);
        QDir().mkdir(tempdir);

        QString input_file(input_files.at(i));

        QFileInfo fi(input_files.at(i));
        if ( fi.suffix().compare(QString("NEF"), Qt::CaseInsensitive) == 0 ||  // Nikon
             fi.suffix().compare(QString("ARW"), Qt::CaseInsensitive) == 0 ||  // Sony
             fi.suffix().compare(QString("PEF"), Qt::CaseInsensitive) == 0 ||  // Pentax
             fi.suffix().compare(QString("IIQ"), Qt::CaseInsensitive) == 0 ||  // Phase One
             fi.suffix().compare(QString("MOS"), Qt::CaseInsensitive) == 0 ||  // Leaf
             fi.suffix().compare(QString("ORF"), Qt::CaseInsensitive) == 0 ||  // Olympus
             fi.suffix().compare(QString("RW2"), Qt::CaseInsensitive) == 0 ||  // Panasonic
             fi.suffix().compare(QString("RAF"), Qt::CaseInsensitive) == 0 ||  // Fujifilm -> bayer mode will probably break horribly
             fi.suffix().compare(QString("DNG"), Qt::CaseInsensitive) == 0 ||  // Pentax/Ricoh, maybe others
             fi.suffix().compare(QString("CR2"), Qt::CaseInsensitive) == 0) { // Canon

            input_file = QString(tempdir + "/" + fi.completeBaseName() + QString(".tiff"));
            QProcess dcp(this);
            dcp.setProgram(dcraw_binary);
            dcp.setStandardOutputFile(input_file);
            if (arguments.contains(QString("--bayer"))) {
                dcp.setArguments(
                    QStringList() << "-4" << "-T" << "-D" << "-c" << input_files.at(i)
                );
            } else {
                dcp.setArguments(
                    QStringList() << "-w" <<  "-4" << "-T" << "-q" << "3" << "-c" << input_files.at(i)
                );
            }

            logger.debug("arguments to dcraw [%s]:\n", dcraw_binary.toLocal8Bit().constData());
            for (int kk = 0; kk < dcp.arguments().size(); kk++) {
                logger.debug("[%d]=%s\n", kk, dcp.arguments().at(kk).toLocal8Bit().constData());
            }
            dcp.start();
            dcp.waitForFinished(-1);
            int dc_rval = dcp.exitStatus() == QProcess::NormalExit && dcp.exitCode() == 0;
            if (!dc_rval) {
                logger.error("Error. dcraw call failed on input image %s [exit status=%d, exitcode=%d]\n", input_files.at(i).toLocal8Bit().constData(), dcp.exitStatus(), dcp.exitCode());
            }
            emit send_delete_item(input_file);
        }
        
        QStringList mma;
        
        // if a "--focal-ratio" setting is already present, then assume this
        // was a user-provided override, otherwise try to calculate it from the EXIF data
        if (!arguments.contains("--focal-ratio")) {
            Exiv2_property props(exiv2_binary, input_files.at(i), tempdir + "/exifinfo.txt");
            mma << "--focal-ratio" << props.get_focal_ratio();
        }

        mma << "--gnuplot-executable " + gnuplot_binary << input_file << tempdir << "--logfile " + tempdir + "/log.txt" 
            << arguments.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        
        Processing_command pc(
            QCoreApplication::applicationDirPath() + "/mtf_mapper",
            mma,
            input_file,
            tempdir,
            input_files.at(i)
        );
        
        if (force_manual_roi_mode) {
            send_processing_command(pc);
        } else {
            process_command(pc, failures);
        }
    }
    // turn off the transient modifiers
    set_single_roi_mode(false); 
    set_focus_mode(false);
    set_imatest_mode(false);
    set_manual_roi_mode(false);
    emit send_progress_indicator(input_files.size()+1);
    emit send_all_done();
    
    for (const auto& failure : failures) {
        emit mtfmapper_call_failed(failure.first, failure.second);
    }
}

void Worker_thread::receive_arg_string(QString s) {
    settings_arguments = s;
}

void Worker_thread::receive_abort() {
    abort = true;
}

QString Worker_thread::update_arguments(QString& s) {
    QString arguments(s);
    
    if (force_roi_mode) {
        if (!arguments.contains("--annotate") && !arguments.contains("-a ")) {
            arguments = arguments + QString(" -a");
        }
        if (!arguments.contains("--edges") && !arguments.contains("-q ")) {
            arguments = arguments + QString(" -q");
        }
        if (!arguments.contains("--single-roi")) {
            arguments = arguments + QString(" --single-roi");
        }
        
        // disable all the other mutually exclusive outputs
        arguments.replace("--lensprofile ", " ");
        arguments.replace("-s ", " ");
        arguments.replace("--surface ", " ");
        arguments.replace("-p ", " ");
        arguments.replace("--profile ", " ");
        arguments.replace("--focus ", " ");
        arguments.replace("--chart-orientation ", " ");
    }
    
    if (force_focus_mode) {
        if (!arguments.contains("--focus")) {
            arguments = arguments + QString(" --focus");
        }
        // disable all the other mutually exclusive outputs
        arguments.replace("--lensprofile ", " ");
        arguments.replace("-q ", " ");
        arguments.replace("--edges ", " ");
        arguments.replace("-a ", " ");
        arguments.replace("--annotate ", " ");
        arguments.replace("-s ", " ");
        arguments.replace("--surface ", " ");
        arguments.replace("-p ", " ");
        arguments.replace("--profile ", " ");
    }
    
    
    if (force_imatest_mode) {
        if (!arguments.contains("--imatest-chart")) {
            arguments = arguments + QString(" --imatest-chart");
        }
    }
    return arguments;
}

void Worker_thread::process_command(const Processing_command& command, vector<std::pair<failure_t, QString>>& failures) {
    
    QProcess process(this);
    process.setProgram(command.program);
    process.setArguments(command.arguments);
    const QString& tempdir = command.tmp_dirname;
    
    bool q_output_requested = process.arguments().contains("-q");
    bool e_output_requested = process.arguments().contains("-e") || process.arguments().contains("--esf");

    logger.debug("arguments to mtf mapper:\n");
    for (int kk = 0; kk < command.arguments.size(); kk++) {
        logger.debug("[%d]=%s\n", kk, command.arguments.at(kk).toLocal8Bit().constData());
    }
    process.start();
    process.waitForFinished(-1);
    int rval = process.exitStatus() == QProcess::NormalExit && process.exitCode() == 0;
    if (!rval) {
        failure_t failure = UNSPECIFIED;
        switch (process.exitCode()) {
        case 2: failure = IMAGE_OPEN_FAILURE; break;
        case 4: failure = NO_TARGETS_FOUND; break;
        case 5: failure = UNSUPPORTED_IMAGE_ENCODING; break;
        default: failure = UNSPECIFIED; break;
        }
        failures.push_back(std::pair<failure_t, QString>(failure, command.img_filename));
    } else {
        emit send_delete_item(tempdir + "/log.txt");

        // this call must come from within the worker thread, since we
        // may have to perform a raw conversion in the worker thread
        // which would cause the display image filename (root of each data set object)
        // to differ from the file containing the exif info
        
        emit send_exif_filename(command.exif_filename, tempdir);
        QString fname(QFileInfo(command.img_filename).completeBaseName());
        emit send_parent_item(fname, command.img_filename);
        
        if (q_output_requested) {
            emit send_delete_item(tempdir + QString("/edge_sfr_values.txt"));
            emit send_delete_item(tempdir + QString("/edge_mtf_values.txt"));
            emit send_delete_item(tempdir + QString("/edge_line_deviation.txt"));
            if (QFile().exists(tempdir + QString("/serialized_edges.bin"))) {
                emit send_delete_item(tempdir + QString("/serialized_edges.bin"));
            }
        }
        
        if (e_output_requested) {
            if (QFile().exists(tempdir + QString("/raw_esf_values.txt"))) {
                emit send_delete_item(tempdir + QString("/raw_esf_values.txt"));
            }
            if (QFile().exists(tempdir + QString("/raw_psf_values.txt"))) {
                emit send_delete_item(tempdir + QString("/raw_psf_values.txt"));
            }
        }
        
        QString an_file = QString("%1/annotated.png").arg(tempdir);
        if (QFile().exists(an_file)) {
            emit send_child_item(QString("annotated"), an_file);
            emit send_delete_item(an_file);
        }
        QString an_file_jpeg = QString("%1/annotated.jpg").arg(tempdir);
        if (QFile().exists(an_file_jpeg)) {
            emit send_child_item(QString("annotated"), an_file_jpeg);
            emit send_delete_item(an_file_jpeg);
        }
        QString pr_file = QString("%1/profile_image.png").arg(tempdir);
        if (QFile().exists(pr_file)) {
            emit send_child_item(QString("profile"), pr_file);
            emit send_delete_item(pr_file);
            emit send_delete_item(tempdir + QString("/profile.gnuplot"));
            emit send_delete_item(tempdir + QString("/profile.txt"));
            emit send_delete_item(tempdir + QString("/profile_peak.txt"));
        }
        QString gi_file = QString("%1/grid_image.png").arg(tempdir);
        if (QFile().exists(gi_file)) {
            emit send_child_item(QString("grid2d"), gi_file);
            emit send_delete_item(gi_file);
            emit send_delete_item(tempdir + QString("/grid.gnuplot"));
            emit send_delete_item(tempdir + QString("/grid.txt"));
        }
        QString gs_file = QString("%1/grid_surface.png").arg(tempdir);
        if (QFile().exists(gs_file)) {
            emit send_child_item(QString("grid3d"), gs_file);
            emit send_delete_item(gs_file);
        }
        QString fp_file = QString("%1/focus_peak.png").arg(tempdir);
        if (QFile().exists(fp_file)) {
            emit send_child_item(QString("focus"), fp_file);
            emit send_delete_item(fp_file);
            emit send_delete_item(tempdir + QString("/profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/profile_points.txt"));
            // we might as well try to delete all the bayer-channel specific outputs, whether they are generated, or not
            emit send_delete_item(tempdir + QString("/green_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/green_profile_points.txt"));
            emit send_delete_item(tempdir + QString("/blue_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/blue_profile_points.txt"));
            emit send_delete_item(tempdir + QString("/red_profile_curve.txt"));
            emit send_delete_item(tempdir + QString("/red_profile_points.txt"));
        }
        QString lp_file = QString("%1/lensprofile.png").arg(tempdir);
        if (QFile().exists(lp_file)) {
            emit send_child_item(QString("lensprofile"), lp_file);
            emit send_delete_item(lp_file);
            emit send_delete_item(tempdir + QString("/lensprofile.txt"));
            emit send_delete_item(tempdir + QString("/lensprofile.gnuplot"));
        }
        QString co_file = QString("%1/chart_orientation.png").arg(tempdir);
        if (QFile().exists(co_file)) {
            emit send_child_item(QString("chart orientation"), co_file);
            emit send_delete_item(co_file);
        }
        QString ca_file = QString("%1/ca_image.png").arg(tempdir);
        if (QFile().exists(ca_file)) {
            emit send_child_item(QString("chromatic aberration"), ca_file);
            emit send_delete_item(ca_file);
            emit send_delete_item(tempdir + QString("/ca_grid.gnuplot"));
            emit send_delete_item(tempdir + QString("/ca_grid.txt"));
        }
        if (QFile().exists(tempdir + QString("/chromatic_aberration.txt"))) {
            emit send_delete_item(tempdir + QString("/chromatic_aberration.txt"));
        }
        QString fids_file = QString("%1/fiducial_correspondence.txt").arg(tempdir);
        if (QFile().exists(fids_file)) {
            emit send_delete_item(fids_file);
        }
        
        emit send_close_item();
    }
}
