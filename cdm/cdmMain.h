#ifndef CDMMAIN_H
#define CDMMAIN_H

#include <QApplication>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif
#include <QtSql>
#include <QtGui>

#include "cdmRunWidget.h"
#include "cdmOptions.h"
#include "cdmOptionsWidget.h"
#include "cdmOptionsWindow.h"
#include "cdmRun.h"
#include "QsLog.h"
#include "QsLogDest.h"
#include "Assistant.h"

class CdmMain : public QMainWindow
{
  
  Q_OBJECT
    
    public:
  
  CdmMain( QString appDirPath = 0, QMainWindow* parent = 0, Qt::WindowFlags fl = Qt::Window );
  virtual ~CdmMain();

 signals:
  void setDbPath(QString path);
  void updatewindow();

 public slots:
 
  void showAppWarning(QString);
  void openFDocWindow();
  void openEDocWindow();
  void openOptionsWindow();
  void openOptionsWidget();
  void saveOptionsWidget();
  void saveAsOptionsWidget();
  void openWidget(Options *options);
  void openFileWindow();
  void saveFileWindow();
  void closeWidgetTab(int index);
  void execute();
  void clearPlots();
  void colorPlots();
  void savePlots();
  void printPlots();
  void updateOptionsWindow();

 protected:
  void closeEvent(QCloseEvent *event);
  void keyPressEvent(QKeyEvent *e);

 private:
  QString            appDirPath;
  Assistant          *assistant;
  QTabWidget         *tab;
  Options            *options;
  OptionsWindow      *optionswindow;
};


#endif // CDMMAIN_H
