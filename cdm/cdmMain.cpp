#include "cdmMain.h"

#define DEBUG_LEVEL QsLogging::InfoLevel

CdmMain::CdmMain( QString _appDirPath, QMainWindow* parent, Qt::WindowFlags fl)
  : QMainWindow( parent, fl )
{
  appDirPath = _appDirPath;
  assistant = new Assistant(appDirPath);
  QDir qdir;
  
  options = new Options();
  options->initDb();

  optionswindow = new OptionsWindow(this,options,Qt::Window);
  optionswindow->setObjectName("OptionsWindow");
  connect(this,SIGNAL(updatewindow()),optionswindow,SLOT(updatewindow()));

  tab = new QTabWidget();
  tab->setObjectName("TabWidget");
  tab->setTabsClosable(true);
  connect(tab, SIGNAL(tabCloseRequested ( int )), this, SLOT(closeWidgetTab(int)));

  this->setCentralWidget(tab);
  this->setMinimumHeight(600);
  this->setMinimumWidth(800);
  
  // Creation of menuBar
  QMenuBar* menuBar = new QMenuBar(this);

  // Creating a menu "File"
  QMenu* menuFile = new QMenu("File");
  // Creating a menu "Calculation"
  QMenu* menuCalculation = new QMenu("Calculation");
  // Creating a menu "Plot"
  QMenu* menuPlot = new QMenu("Plot");
  // Creating a menu "Help"
  QMenu* menuHelp = new QMenu("Help");

  // And adding the menu "File" in the menu bar
  menuBar->addMenu(menuFile);
  // And adding the menu "Calculation" in the menu bar
  menuBar->addMenu(menuCalculation);
  // And adding the menu "Plot" in the menu bar
  menuBar->addMenu(menuPlot);
  // And adding the menu "Help" in the menu bar
  menuBar->addMenu(menuHelp);
  // Add a separator
  menuFile->addSeparator();
  // Add a separator
  menuCalculation->addSeparator();
  // Add a separator
  menuPlot->addSeparator();

 
  //  menuFile->addAction("Open File", this, SLOT(openFileWindow()) );
  //  menuFile->addAction("Save File", this, SLOT(saveFileWindow()) );
  menuFile->addAction("Exit", this, SLOT(close()) );
  menuHelp->addAction("Documentation-FR",this,SLOT(openFDocWindow()));
  menuHelp->addAction("Documentation-EN",this,SLOT(openEDocWindow()));
  menuCalculation->addAction("New", this, SLOT(openOptionsWidget()));
  menuCalculation->addAction("Load", this, SLOT(openOptionsWindow()));
  menuCalculation->addAction("Save", this, SLOT(saveOptionsWidget()));
  menuCalculation->addAction("Save as", this, SLOT(saveAsOptionsWidget()));
  menuCalculation->addAction("Start", this, SLOT(execute()));
  menuPlot->addAction("Clear", this, SLOT(clearPlots()));
  menuPlot->addAction("Color", this, SLOT(colorPlots()));
  menuPlot->addAction("Save", this, SLOT(savePlots()));
  menuPlot->addAction("Print", this, SLOT(printPlots()));
  // Then set the menu bar to the main window
  setMenuBar(menuBar);
}

void CdmMain::openFileWindow() {
}
void CdmMain::saveFileWindow() {
}
void CdmMain::openOptionsWidget() {
  //Options *options = new Options();
  //options->initDb();
  options->loadDb("");
  openWidget(options);
}
void CdmMain::saveOptionsWidget() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  OptionsWidget *optionswidget = NULL;
  if (widget) {
    optionswidget = widget->findChild<OptionsWidget *>("Options");
    optionswidget->updateOptions();
    optionswidget->saveName();
    updateOptionsWindow();
  }
}
void CdmMain::saveAsOptionsWidget() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  OptionsWidget *optionswidget = NULL;
  if (widget) {
    optionswidget = widget->findChild<OptionsWidget *>("Options");
    optionswidget->updateOptions();
    optionswidget->saveAsName();
    updateOptionsWindow();
  }
}
void CdmMain::updateOptionsWindow() {
 emit updatewindow();
}
void CdmMain::openOptionsWindow() {
   optionswindow->show();
   optionswindow->raise();
}
void CdmMain::execute() {
   QMainWindow *widget = (QMainWindow*)tab->currentWidget();
   RunWidget *runwidget = NULL;
   OptionsWidget *optionswidget = NULL;
   if (widget) {
    runwidget = widget->findChild<RunWidget *>("Run");
    optionswidget = widget->findChild<OptionsWidget *>("Options");
    optionswidget->updateOptions();
    runwidget->execute();
   }
}
void CdmMain::clearPlots() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  RunWidget *runwidget = NULL;
  if (widget) {
    runwidget = widget->findChild<RunWidget *>("Run");
    runwidget->cleanupPlots();
  }
}
void CdmMain::colorPlots() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  RunWidget *runwidget = NULL;
  if (widget) {
    runwidget = widget->findChild<RunWidget *>("Run");
    runwidget->colorPlots();
  }
}
void CdmMain::savePlots() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  RunWidget *runwidget = NULL;
  if (widget) {
    runwidget = widget->findChild<RunWidget *>("Run");
    runwidget->savePlots();
  }
}
void CdmMain::printPlots() {
  QMainWindow *widget = (QMainWindow*)tab->currentWidget();
  RunWidget *runwidget = NULL;
  if (widget) {
    runwidget = widget->findChild<RunWidget *>("Run");
    runwidget->printPlots();
  }
}
void CdmMain::openFDocWindow()
{
  //assistant->showDocumentation("index.html");
  QLOG_INFO() << "Open doc at " << "file:///"+QDir::currentPath()+"/../doc/userguide-FR.pdf";
  QDesktopServices::openUrl(QUrl("file:///"+ QDir::currentPath()+"/../doc/userguide-FR.pdf"));
}

void CdmMain::openEDocWindow()
{
  //assistant->showDocumentation("index.html");
  QLOG_INFO() << "Open doc at " << "file:///"+QDir::currentPath()+"/../doc/userguide-EN.pdf";
  QDesktopServices::openUrl(QUrl("file:///"+ QDir::currentPath()+"/../doc/userguide-EN.pdf"));
}
void CdmMain::openWidget(Options *options) {
    // check if tab is already there
    for (int i = 0; i < tab->count(); i++) { 
     if (tab->tabText(i) == options->getName())  {
        tab->setCurrentIndex(i);
        QMainWindow *widget = (QMainWindow*)tab->currentWidget();
        OptionsWidget *optionswidget = NULL;
        optionswidget = widget->findChild<OptionsWidget *>("Options");
        optionswidget->update();
        return;
     }
    }  
    QMainWindow *widget = new QMainWindow();
    widget->setDockNestingEnabled(true);
    QDockWidget *dockWidget = new QDockWidget(options->getName(), widget);
    dockWidget->setObjectName("DockWidget");
    dockWidget->setFeatures(QDockWidget::DockWidgetMovable |
                                      QDockWidget::DockWidgetFloatable);
    dockWidget->setAllowedAreas(Qt::LeftDockWidgetArea);
    QScrollArea *optionsArea = new QScrollArea(widget);
    OptionsWidget *optionswidget = new OptionsWidget(this,options);
    optionsArea->setWidget(optionswidget);
    optionsArea->setWidgetResizable(true);
    optionsArea->setMinimumWidth(400);
    dockWidget->setWidget(optionsArea);
    widget->addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    QScrollArea *runArea = new QScrollArea(widget);
    Run *run = new Run(options->getName());
    RunWidget *runwidget = new RunWidget(options, run);
    runArea->setWidget(runwidget);
    runArea->setWidgetResizable(true);
    widget->setCentralWidget(runArea);
    tab->setCurrentIndex(tab->addTab(widget,options->getName()));
   
}
void CdmMain::closeWidgetTab(int index)
{
  QMainWindow *widget = (QMainWindow*)tab->widget(index);
  RunWidget *runwidget = NULL;
  runwidget = widget->findChild<RunWidget *>("Run");
  if ( runwidget != NULL )
    delete runwidget;
  OptionsWidget *optionswidget = NULL;
  optionswidget = widget->findChild<OptionsWidget *>("Options");
  if ( optionswidget != NULL )
    delete optionswidget;
  tab->removeTab(index);
}
CdmMain::~CdmMain()
{
  delete tab;
  delete assistant;
  QLOG_INFO() << "cdm application ended";
}

void CdmMain::keyPressEvent(QKeyEvent *e) {
  if(e->type() == QKeyEvent::KeyPress) {
    if(e->matches(QKeySequence::Copy)) {
      delete this;
    } 
  }
}
  
void CdmMain::closeEvent(QCloseEvent* event)
{
  event->accept();
  delete this;
}

void 
CdmMain::showAppWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}

int main(int argc, char *argv[])
{
  QApplication app(argc, argv); 
  app.addLibraryPath("/usr/local/bin");
  // init the logging mechanism
  QsLogging::Logger& logger = QsLogging::Logger::instance();
  const QString sLogPath(QDir::currentPath() +  QDir::separator() +
                        "cdm_" + 
                         QDateTime::currentDateTime().toString("MMMdd,yy-hh:mm:ss") +
                         ".log");
  QsLogging::DestinationPtr fileDestination(QsLogging::DestinationFactory::MakeFileDestination(sLogPath) );
  QsLogging::DestinationPtr debugDestination(QsLogging::DestinationFactory::MakeDebugOutputDestination() );
  logger.addDestination(debugDestination.get());
  logger.addDestination(fileDestination.get());
  logger.setLoggingLevel(DEBUG_LEVEL);
  
  QLOG_INFO() << "cdm version " << CDMVERSION
	      << " started : " <<  app.applicationDirPath();
  QLOG_INFO() << "Built with Qt" << QT_VERSION_STR << "running on" << qVersion();
  QLOG_INFO() << "Qt User Data location : " 
              <<  QDir::currentPath();
  CdmMain* cdm = new CdmMain(app.applicationDirPath(),NULL,NULL);
  cdm->setWindowTitle("discrete dipole approximation  in free space");
  cdm->show();
  foreach (const QString &path, app.libraryPaths())
  QLOG_DEBUG() << path;
 
  return app.exec();
}
