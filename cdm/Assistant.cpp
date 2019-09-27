/*******************************************************************
This file is part of cdm.

cdm is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
********************************************************************/

#include <QtCore/QByteArray>
#include <QtCore/QDir>
#include <QtCore/QLibraryInfo>
#include <QtCore/QProcess>
#if QT_VERSION < 0x050000
#include <QMessageBox>
#else
#include <QtWidgets>
#endif
#include "Assistant.h"


 Assistant::Assistant(QString _appDirPath)
     : proc(0)
 {
   appDirPath = _appDirPath;
 }

 Assistant::~Assistant()
 {
     if (proc && proc->state() == QProcess::Running) {
         proc->terminate();
         proc->waitForFinished(3000);
     }
     delete proc;
 }

 void Assistant::showDocumentation(const QString &page)
 {
     if (!startAssistant())
         return;

     QByteArray ba("SetSource ");
     //ba.append("qthelp://com.trolltech.examples.simpletextviewer/doc/");

     proc->write(page.toLocal8Bit() + '\n');
 }

 bool Assistant::startAssistant()
 {
     if (!proc)
         proc = new QProcess();

     if (proc->state() != QProcess::Running) {
         QString app = QLibraryInfo::location(QLibraryInfo::BinariesPath) + QDir::separator();
 #if !defined(Q_OS_MAC)
         app += QLatin1String("assistant");
 #else
         app += QLatin1String("Assistant.app/Contents/MacOS/Assistant");
 #endif
         QStringList args;
         args << QLatin1String("-collectionFile")
             << appDirPath + QLatin1String("/doc/cdmColl.qhc")
             << QLatin1String("-enableRemoteControl");
         proc->start(app, args);
         if (!proc->waitForStarted()) {
             QMessageBox::critical(0, QObject::tr("cdm"),
                 QObject::tr("Unable to launch Qt Assistant (%1)").arg(app));
             return false;
         }
     }
     return true;
 }
