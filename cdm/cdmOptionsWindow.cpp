#include "cdmOptionsWindow.h"
#include "cdmMain.h"

OptionsWindow::OptionsWindow( QMainWindow *_parent, Options *_options, Qt::WindowFlags fl)
  : QMainWindow( _parent, fl )
{
  parent = _parent;
  optionswidget = NULL;
  optionstable = NULL;
  options = _options;
  
  QWidget *centralWidget = new QWidget(this);
  this->setCentralWidget(centralWidget);
  vboxlayout = new QVBoxLayout();
  centralWidget->setLayout(vboxlayout);
  this->setWindowTitle("Configurations saved");
  this->setMinimumHeight(300);
  this->setMinimumWidth(400);
  
  optionstitle = new QLabel("Select a configuration:");
  vboxlayout->addWidget(optionstitle);
  optionsview = new QTableView();
  QPalette palette = optionsview->palette();
  palette.setBrush(QPalette::Base, Qt::transparent);
  optionsview->setPalette(palette);
  optionsview->setItemDelegate(new ColorDelegate());
  optionsview->setEditTriggers(QAbstractItemView::AllEditTriggers);
  //optionsview->setSelectionBehavior(QAbstractItemView::SelectRows);
  vboxlayout->addWidget(optionsview);

  saveButton = new QPushButton("Save",this);
  saveButton->setToolTip("Enter a new name or select a line to save the current calculation");
  saveButton->setFixedHeight(30);
  saveButton->setFixedWidth(150);
  QObject::connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  vboxlayout->addWidget(saveButton);
  removeButton = new QPushButton("Delete",this);
  removeButton->setToolTip("Select a line to delete");
  removeButton->setFixedHeight(30);
  removeButton->setFixedWidth(150);
  QObject::connect(removeButton, SIGNAL(clicked()), this, SLOT(remove()));
  vboxlayout->addWidget(removeButton);
  loadButton = new QPushButton("Load",this);
  loadButton->setToolTip("Select a line to load");
  loadButton->setFixedHeight(30);
  loadButton->setFixedWidth(150);
  QObject::connect(loadButton, SIGNAL(clicked()), this, SLOT(load()));
  vboxlayout->addWidget(loadButton);
  exportButton = new QPushButton("Export",this);
  exportButton->setToolTip("Select a line to export");
  exportButton->setFixedHeight(30);
  exportButton->setFixedWidth(150);
  QObject::connect(exportButton, SIGNAL(clicked()), this, SLOT(tofile()));
  vboxlayout->addWidget(exportButton);
  initOptionsTbl();
} 

OptionsWindow::~OptionsWindow()
{
  QLOG_DEBUG ( ) << "Deleting OptionsWindow";
  if (optionstable){ delete optionstable; optionstable = NULL;}
}

void OptionsWindow::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing OptionsWindow";
  this->hide();
}
void OptionsWindow::initOptionsTbl() {
  QTabWidget *tabwidget = NULL;
  QMainWindow *currentWindow = NULL;
  tabwidget = (QTabWidget*) parent->findChild<QTabWidget*>("TabWidget");
  if (tabwidget) {
   currentWindow = (QMainWindow*)tabwidget->currentWidget();
   optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
  }
  /*if (optionswidget)
    options = optionswidget->options;
  else {
   if (options) {delete options; options = NULL;}
   options = new Options();
   options->initDb();
  } */
  if (optionstable) {delete optionstable; optionstable = NULL;}
  dbpath = "options.db3";
  optionstable = new ManifestModel(this,QSqlDatabase::database(dbpath));
  optionstable->setTable("options_tbl");
  optionstable->setEditStrategy(QSqlTableModel::OnManualSubmit);
  optionsview->setModel(optionstable);
  optionstable->setFilter("name not like 'new'");
  optionstable->select();
  optionstable->insertRow(optionstable->rowCount());
  optionsview->resizeColumnsToContents();
  optionsview->resizeRowsToContents();
  optionsview->setCurrentIndex (optionsview->model()->index(optionstable->rowCount() - 1, 
		    0, QModelIndex()));
}
void 
OptionsWindow::save(){
  if (optionsview->selectionModel()->hasSelection() == false) {
     this->showWarning("Select a valid configuration to save");
     return;
  }
  QModelIndex selected = optionsview->selectionModel()->currentIndex();
  QSqlRecord record = optionstable->record(selected.row());
  QString name = record.value(0).toString();
  QString description = record.value(1).toString();
  if (name =="") {
   this->showWarning("Select a valid configuration to save");
   return;
  }
  QTabWidget *tabwidget = (QTabWidget*) parent->findChild<QTabWidget*>("TabWidget");
  QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
  optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
  if (optionswidget)
    optionswidget->updateOptions();
  else
    options->loadDb(name);
  options->saveDb(name,description);
  initOptionsTbl();
}
void 
OptionsWindow::load(){
  if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to load");
  return;
  }
  QModelIndexList selectedlist = optionsview->selectionModel()->selectedIndexes();
  for (int i = selectedlist.size() - 1 ; i >=0 ; i--) {
    QModelIndex selected = selectedlist.at(i);
    QSqlRecord record = optionstable->record(selected.row());
    QString name = record.value(0).toString();
    options->loadDb(name);
    CdmMain *mainwindow = (CdmMain*) parent; 
    mainwindow->openWidget(options);
 }
 this->hide();
}
void 
OptionsWindow::remove(){
  if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to delete");
    return;
  }
  QModelIndexList selectedlist = optionsview->selectionModel()->selectedIndexes();
  for (int i = selectedlist.size() - 1 ; i >=0 ; i--) {
    QModelIndex selected = selectedlist.at(i);
    QSqlRecord record = optionstable->record(selected.row());
    QString name = record.value(0).toString();
    QString description = record.value(1).toString();
    if (name =="") {
      this->showWarning("Select a valid configuration to delete");
      return;
    }
    options->removeDb(name);
  }
  initOptionsTbl();
}
void 
OptionsWindow::updatewindow(){
  initOptionsTbl();
}
void 
OptionsWindow::tofile(){
 if (optionsview->selectionModel()->hasSelection() == false) {
    this->showWarning("Select a valid configuration to export");
    return;
  }
  QModelIndex selected = optionsview->selectionModel()->currentIndex();
  QSqlRecord record = optionstable->record(selected.row());
  QString name = record.value(0).toString();
  if (name =="") {
   this->showWarning("Select a valid configuration to export");
   return;
  }
  options->loadDb(name);
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                     name + ".opt", tr("Options (*.opt)"));
  QFile optfile(fileName);
  if (!optfile.open(QIODevice::WriteOnly | QIODevice::Text))
     return;
  QTextStream opt(&optfile);
  // Fill ASCII options here
  opt << "Calculation options [0=rigorous, 1=renormalized Born, 2=Born, 3=Born order, 4=renormalized Rytov, 5=Rytov, 6=Beam propagation, 7=renormalized Beam propagation]:" << options->getNrig() << endl;
  opt << "Illumination properties options:" << endl;
  opt << " wavelength:" << options->getWavelength() << endl;
  opt << " P0: " << options->getP0() << endl;
  opt << " W0: " << options->getW0() << endl;
  opt << " Beam: " << options->getBeam() << endl;
  if (options->getBeam() == "Circular plane wave") {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization L (-1) R (1):" << options->getPolarizationRL() << endl;
  }
  else if (options->getBeam() == "Linear plane wave") {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization TM (1) TE (0):" << options->getPolarizationTM() << endl;
  }
  else if (options->getBeam() == "Circular Gaussian" || options->getBeam() == "Circular Gaussian (para)"|| options->getBeam() == "Circular Gaussian (FFT)" ) {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization L (-1) R (1):" << options->getPolarizationRL() << endl;
    opt << "  gaussian X:" << options->getXgaus() << endl;
    opt << "  gaussian Y:" << options->getYgaus() << endl;
    opt << "  gaussian Z:" << options->getZgaus() << endl;
  }
  else if (options->getBeam() == "Linear Gaussian" || options->getBeam() == "Linear Gaussian (para)" || options->getBeam() == "Linear Gaussian (FFT)" ) {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  incidence angle (phi with respect to x):" << options->getIncidenceangle_phi_x() << endl;
    opt << "  polarization TM (1) TE (0):" << options->getPolarizationTM() << endl;
    opt << "  gaussian X:" << options->getXgaus() << endl;
    opt << "  gaussian Y:" << options->getYgaus() << endl;
    opt << "  gaussian Z:" << options->getZgaus() << endl;
  }
  else if (options->getBeam() == "Multiple wave"  ) {
    for (int i = 0 ; i < options->getWaveMultiNumber(); i++) {
      opt << "  incidence angle (theta with respect to z):" << options->getThetam().at(i) << endl;
      opt << "  incidence angle (phi with respect to x):" << options->getPhim().at(i) << endl;
      opt << "  polarization TM (1) TE (0):" << options->getPpm().at(i) << endl;
      opt << "  Magnitude field real:" << (real(options->getE0m().at(i))) << endl;
      opt << "  Magnitude field imag:" << (imag(options->getE0m().at(i))) << endl;
    }
  }
  else if (options->getBeam() == "Antenna" ) {
    opt << "  incidence angle (theta with respect to z):" << options->getIncidenceangle_theta_z() << endl;
    opt << "  orientation angle (phi with respect to x)" << options->getIncidenceangle_phi_x() << endl;
    opt << "  position X:" << options->getXgaus() << endl;
    opt << "  position Y:" << options->getYgaus() << endl;
    opt << "  position Z:" << options->getZgaus() << endl;
  }

  opt << "Object properties options:" << endl;
  for (int i = 0 ; i < options->getObjectNumber(); i++) {
    opt << " object : " << i << " : " << options->getObject() << endl;
    if (options->getObject() == "sphere" || options->getObject() == "multiple spheres") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "inhomogeneous sphere") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  coherence length:" << options->getSpherecoherencelength() << endl;
      opt << "  standard deviation:" << options->getSpherestandardev() << endl;
    }
    else if (options->getObject() == "random spheres (length)") {
      opt << "  cube side X:" << options->getCubesidex() << endl;
      opt << "  cube side Y:" << options->getCubesidey() << endl;
      opt << "  cube side Z:" << options->getCubesidez() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  density:" << options->getDensity() << endl;
     }
    else if (options->getObject() == "random spheres (meshsize)") {
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  number of subunit X:" << options->getNxx() << endl;
      opt << "  number of subunit Y:" << options->getNyy() << endl;
      opt << "  number of subunit Z:" << options->getNzz() << endl;
      opt << "  meshsize:" << options->getMeshsize() << endl;      
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  density:" << options->getDensity() << endl;
     }
    else if (options->getObject() == "concentric spheres") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      if (i == 0) {
        opt << "  position X:" << options->getPositionx().at(i) << endl;
        opt << "  position Y:" << options->getPositiony().at(i) << endl;
        opt << "  position Z:" << options->getPositionz().at(i) << endl;
      }
    }
    else if (options->getObject() == "cube") {
      opt << "  cube side:" << options->getCubeside() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "inhomogeneous cuboid (length)") {
      opt << "  cube side X:" << options->getCubesidex() << endl;
      opt << "  cube side Y:" << options->getCubesidey() << endl;
      opt << "  cube side Z:" << options->getCubesidez() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  coherence length:" << options->getSpherecoherencelength() << endl;
      opt << "  standard deviation:" << options->getSpherestandardev() << endl;
    }
       else if (options->getObject() == "inhomogeneous cuboid (meshsize)") {
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  number of subunit X:" << options->getNxx() << endl;
      opt << "  number of subunit Y:" << options->getNyy() << endl;
      opt << "  number of subunit Z:" << options->getNzz() << endl;
      opt << "  meshsize:" << options->getMeshsize() << endl;
      opt << "  seed:" << options->getSphereseed() << endl;
      opt << "  coherence length:" << options->getSpherecoherencelength() << endl;
      opt << "  standard deviation:" << options->getSpherestandardev() << endl;
    }
      else if (options->getObject() == "cuboid (length)") {
      opt << "  cube side X:" << options->getCubesidex() << endl;
      opt << "  cube side Y:" << options->getCubesidey() << endl;
      opt << "  cube side Z:" << options->getCubesidez() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  theta:" << options->getThetaobj() << endl;
      opt << "  phi:" << options->getPhiobj() << endl;
      opt << "  psi:" << options->getPsiobj() << endl;
    }
      else if (options->getObject() == "cuboid (meshsize)") {
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  number of subunit X:" << options->getNxx() << endl;
      opt << "  number of subunit Y:" << options->getNyy() << endl;
      opt << "  number of subunit Z:" << options->getNzz() << endl;
      opt << "  meshsize:" << options->getMeshsize() << endl;      
    }
    else if (options->getObject() == "ellipsoid") {
      opt << "  half axe A:" << options->getDemiaxea() << endl;
      opt << "  half axe B:" << options->getDemiaxeb() << endl;
      opt << "  half axe C:" << options->getDemiaxec() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
    }
    else if (options->getObject() == "cylinder") {
      opt << "  radius:" << options->getSphereradius().at(i) << endl;
      opt << "  height:" << options->getHauteur() << endl;
      opt << "  position X:" << options->getPositionx().at(i) << endl;
      opt << "  position Y:" << options->getPositiony().at(i) << endl;
      opt << "  position Z:" << options->getPositionz().at(i) << endl;
      opt << "  theta:" << options->getThetaobj() << endl;
      opt << "  phi:" << options->getPhiobj() << endl;
    }
  }
    opt << " anisotropy:" << options->getAnisotropy() << endl;
  for (int i = 0 ; i < options->getObjectNumber(); i++) {
    if ( options->getAnisotropy() == "iso" ) {
      opt << "  material" << i << ":" << options->getMaterial()[i] << endl;
      opt << "  epsilon real:" << QString::number(real(options->getEpsilon().at(i))) << endl;
      opt << "  epsilon imag:" << QString::number(imag(options->getEpsilon().at(i))) << endl;
    }
    else if ( options->getAnisotropy() == "ani" ) {
      opt << "  epsilon11 real:" << QString::number(real(options->getEpsilon11().at(i))) << endl;
      opt << "  epsilon11 imag:" << QString::number(imag(options->getEpsilon11().at(i))) << endl;
      opt << "  epsilon12 real:" << QString::number(real(options->getEpsilon12().at(i))) << endl;
      opt << "  epsilon12 imag:" << QString::number(imag(options->getEpsilon12().at(i))) << endl;
      opt << "  epsilon13 real:" << QString::number(real(options->getEpsilon13().at(i))) << endl;
      opt << "  epsilon13 imag:" << QString::number(imag(options->getEpsilon13().at(i))) << endl;
      opt << "  epsilon21 real:" << QString::number(real(options->getEpsilon21().at(i))) << endl;
      opt << "  epsilon21 imag:" << QString::number(imag(options->getEpsilon21().at(i))) << endl;
      opt << "  epsilon22 real:" << QString::number(real(options->getEpsilon22().at(i))) << endl;
      opt << "  epsilon22 imag:" << QString::number(imag(options->getEpsilon22().at(i))) << endl;
      opt << "  epsilon23 real:" << QString::number(real(options->getEpsilon23().at(i))) << endl;
      opt << "  epsilon23 imag:" << QString::number(imag(options->getEpsilon23().at(i))) << endl;
      opt << "  epsilon31 real:" << QString::number(real(options->getEpsilon31().at(i))) << endl;
      opt << "  epsilon31 imag:" << QString::number(imag(options->getEpsilon31().at(i))) << endl;
      opt << "  epsilon32 real:" << QString::number(real(options->getEpsilon32().at(i))) << endl;
      opt << "  epsilon32 imag:" << QString::number(imag(options->getEpsilon32().at(i))) << endl;
      opt << "  epsilon33 real:" << QString::number(real(options->getEpsilon33().at(i))) << endl;
      opt << "  epsilon33 imag:" << QString::number(imag(options->getEpsilon33().at(i))) << endl;
    }
  }
  opt << "Study options:" << endl;
  opt << " Dipole/epsilon checked:" << options->getDipolepsilon() << endl;
  opt << " Farfield checked:" << options->getFarfield() << endl;
  if ( options->getFarfield() ) {
    opt << "  cross section checked:" << options->getCrosssection() << endl;
    opt << "  cross section + poynting checked:" << options->getCrosssectionpoynting() << endl;
    opt << "  quick computation:" << options->getQuickdiffract() << endl;
    opt << "  ntheta:" << options->getNtheta() << endl;
    opt << "  nphi:" << options->getNphi() << endl;
    opt << "  Emissivity:" << options->getNenergie() << endl;
  }
  opt << " Microscopy checked:" << options->getMicroscopy() << endl;
  if ( options->getMicroscopy() ) {
    opt << "  Microscope [0=Holographic, 1=Brightfield, 2=Darkfield & phase]:" << options->getNtypemic() << endl;
    opt << "  Side of observation [0 trans, 1 ref.]:" << options->getNside() << endl;
    opt << "  Quick computation:" << options->getMicroscopyFFT() << endl;
    opt << "  Numerical aperture lens:" << options->getNA() << endl;
    opt << "  Numerical aperture condenser:" << options->getNAinc() << endl;
    opt << "  Magnification:" << options->getGross() << endl;
    opt << "  Focal plane position:" << options->getZlens() << endl;
  }
  opt << " Force checked:" << options->getForce() << endl;
  if ( options->getForce() ) {
    opt << "  optical force checked:" << options->getOpticalforce() << endl;
    opt << "  optical force density checked:" << options->getOpticalforcedensity() << endl;
    opt << "  optical torque checked:" << options->getOpticaltorque() << endl;
    opt << "  optical torque density checked:" << options->getOpticaltorquedensity() << endl;
  }
  opt << " Nearfield checked:" << options->getNearfield() << endl;
  if ( options->getNearfield() ) {
    opt << "  localfield checked:" << options->getLocalfield() << endl;
    opt << "  macroscopic field checked:" << options->getMacroscopicfield() << endl;
    int nproche = options->getNproche();
    opt << "  range of study (-1,0,1):" << nproche << endl;
    opt << "  nxm:" << options->getNxm() << endl;
    opt << "  nym:" << options->getNym() << endl;
    opt << "  nzm:" << options->getNzm() << endl;
    opt << "  nxmp:" << options->getNxmp() << endl;
    opt << "  nymp:" << options->getNymp() << endl;
    opt << "  nzmp:" << options->getNzmp() << endl;
  }

  opt << "Numerical parameters options:" << endl;
  opt << " Tolerance:" << options->getTolerance() << endl;
  opt << " Methode:" << options->getMethodeit() << endl;
  opt << " Polarizability:" << options->getPolarizability() << endl;
  opt << " Quadrature:" << options->getQuad() << endl;
  opt << " FFT:" << options->getnfft2d() << endl;
  optfile.close();
  initOptionsTbl();
}
void 
OptionsWindow::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}

