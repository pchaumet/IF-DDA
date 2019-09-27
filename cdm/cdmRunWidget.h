#include <iostream>
#include <string>
#include <fstream>
using namespace std;
#ifndef RUNWIDGET_H
#define RUNWIDGET_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#include <QtConcurrent>
#endif


#include "QsLog.h"
#include "cdmOptions.h"
#include "cdmRun.h"
#include "cdmlib.h"

#include <cdmPlot.h>

#define DELTA 1e-12

void cdmlibwrapper(Options *options, Run *run, QString *infoMessage, int *stopFlag);

class RunWidget : public QMainWindow
{
  Q_OBJECT
    
    public:
      RunWidget(Options *_options, Run *_run);
      ~RunWidget();

     void execute();

    public slots:
      void run_thread();
      void cancel_thread();
      void displayFinishedBox();
      void plotmicroscopy();
      void plotxnearfield();
      void plotynearfield();
      void plotznearfield();
      void plotxforce();
      void plotyforce();
      void plotzforce();
      void plotallxnearfield();
      void plotallynearfield();
      void plotallznearfield();
      void plotallxforce();
      void plotallyforce();
      void plotallzforce();
      void plotallx(QString field, QComboBox *xlist);
      void plotally(QString field, QComboBox *ylist);
      void plotallz(QString field, QComboBox *zlist);
      void plotx(QString field, double xvalue);
      void ploty(QString field, double yvalue);
      void plotz(QString field, double zvalue);
      void dipoles3Dplot();
      void poyntingplot();
      void cleanupPlots();
      void colorPlots();
      void savePlots();
      void printPlots();
      void closeWidgetTab(int index);
      void dockWidgetPlot_topLevelChanged(bool isFloating);
      void dockWidgetOutput_topLevelChanged(bool isFloating);
      void dockWidgetCentral_topLevelChanged(bool isFloating);

    private:
      Options              *options;
      Run                  *run;
      QFuture<void>        *future;
      QFutureWatcher<void> *watcher;
      QProgressDialog      *progressdlg;
      QProgressDialog      *canceldlg;
      QString              *infoMessage;
      int                  *stopFlag;
      QScrollArea          *centralArea;
      QScrollArea          *outputArea;
      QScrollArea          *plotArea;
      QWidget              *centralwidget;
      QDockWidget          *dockWidgetPlot;
      QDockWidget          *dockWidgetOutput;
      QDockWidget          *dockWidgetCentral;
      QFormLayout          *centrallayout;
      QWidget              *outputwidget;
      QFormLayout          *outputlayout;
      QTabWidget           *plotwidget;
      QComboBox            *nearfieldComboBox, *nearfieldtypeComboBox;
      QComboBox            *farfieldComboBox, *farfieldtypeComboBox;
      QComboBox            *forceComboBox, *forcetypeComboBox;
      QComboBox            *znearfieldComboBox, *ynearfieldComboBox, *xnearfieldComboBox;
      QComboBox            *zfarfieldComboBox, *yfarfieldComboBox, *xfarfieldComboBox;
      QComboBox            *zforceComboBox, *yforceComboBox, *xforceComboBox;
      PlotRaster           *plot, *plotpc, *plotmodulus, *plotphase;
      PlotVector           *plotv;

      double               xstepsize, num_col;
      int   nphi, ntheta, nxm, nym, nzm,  nxmp, nymp, nzmp, nxx, nyy, nzz;
      void displayResults();      
      void plotAxislist(QString field, double *X, double *Y, double *Z, 
			double *XWF, double *YWF, double *ZWF, 
			double *REFCOL, double *REFCOLWF, double refaxis,
			QString title1, QString title2, QString title3, QString title4);
      
};
#endif
