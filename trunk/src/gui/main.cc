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
#include <QApplication>

#include "include/logger.h"
Logger logger;

#include <QFileInfo>
#include <QSurfaceFormat>

#include "mtfmapper_app.h"

void message_output(QtMsgType type, const QMessageLogContext& /*context*/, const QString &msg) {
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        logger.debug("Debug: %s\n", localMsg.constData());
        break;
    case QtInfoMsg:
        logger.debug("Info: %s\n", localMsg.constData());
        break;
    case QtWarningMsg:
        logger.debug("Warning: %s\n", localMsg.constData());
        break;
    case QtCriticalMsg:
        logger.debug("Critical: %s\n", localMsg.constData());
        break;
    case QtFatalMsg:
        logger.error("Fatal: %s\n", localMsg.constData());
        abort();
    }
}
 
int main(int argc, char *argv[]) {

    string logname = (QDir::tempPath() + QDir::separator()).toStdString() + "mtfmapperlog.txt";
    logger.redirect(logname);
    logger.enable_level(Logger::DEBUG);
    qInstallMessageHandler(message_output);

    QApplication app(argc, argv);
    
    QSurfaceFormat fmt;
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    fmt.setMajorVersion(2);
    fmt.setMinorVersion(1);
    fmt.setSamples(1);
    QSurfaceFormat::setDefaultFormat(fmt);
    
    mtfmapper_app dialog;

    app.setAttribute(Qt::AA_UseHighDpiPixmaps);

    QIcon appIcon;
    appIcon.addFile(":/Icons/AppIcon64");
    appIcon.addFile(":/Icons/AppIcon128");
    appIcon.addFile(":/Icons/AppIcon256");
    app.setWindowIcon(appIcon);

    dialog.show();
    return app.exec();
}
