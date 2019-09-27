#include "cdmOptionsWidget.h"

OptionsWidget::OptionsWidget(QMainWindow *_mainwindow, Options *_options)
{

  mainwindow = _mainwindow;
  options = _options;
  name = options->getName();
  description = options->getDescription();
  this->setObjectName("Options");
  connect(this, SIGNAL(updateOptionsWindow()),mainwindow,SLOT(updateOptionsWindow()));

  cfgWindow = new QDialog(this,Qt::Dialog);
  cfgWindow->setModal(true);
  QLabel *nameLabel = new QLabel("Name:");
  nameLineEdit = new QLineEdit;
  QLabel *descriptionLabel = new QLabel("Description:");
  descriptionLineEdit = new QLineEdit;
  QPushButton *finish = new QPushButton("Ok");
  connect(finish,SIGNAL(clicked()),this,SLOT(finish()));
  QPushButton *cancel = new QPushButton("Cancel");
  connect(cancel,SIGNAL(clicked()),this,SLOT(cancel()));
  QGridLayout *cfgWindowlayout = new QGridLayout;
  cfgWindowlayout->addWidget(nameLabel, 0, 0);
  cfgWindowlayout->addWidget(nameLineEdit, 0, 1);
  cfgWindowlayout->addWidget(descriptionLabel, 1, 0);
  cfgWindowlayout->addWidget(descriptionLineEdit, 1, 1);
  cfgWindowlayout->addWidget(finish, 2, 2);
  cfgWindowlayout->addWidget(cancel, 2, 3);
  cfgWindow->setLayout(cfgWindowlayout);
  cfgWindow->setWindowTitle("Save new configuration");

  layout = new QFormLayout(this);

  QHBoxLayout *calculationlayout = new QHBoxLayout();
  executeButton = new QPushButton("Start calculation");
  executeButton->setStyleSheet("background-color:red;");
 
  executeButton->setFixedWidth(150);
  connect(executeButton, SIGNAL(clicked()), this, SLOT(execute()));
  saveButton = new QPushButton("Save configuration");
  saveButton->setFixedWidth(150);
  connect(saveButton, SIGNAL(clicked()), this, SLOT(saveName()));
  calculationlayout->addWidget(executeButton);
  calculationlayout->addWidget(saveButton);

  QHBoxLayout *startnreadlayout = new QHBoxLayout();
  nrigLabel = new QLabel("Calculation options");
  nrig      = new QComboBox();
  nrig->addItems(options->nrigList);
  nrig->setCurrentIndex(options->getNrig());
  nrig->setFixedWidth(200);
  nreadLabel = new QLabel("Read local field from file");
  nread = new QCheckBox(this);
  connect(nread, SIGNAL(stateChanged(int)),this,
	SLOT(nreadCheckBoxStateChanged(int)));
  filerereadLabel = new QLabel("File name:");
  filereread = new QLineEdit(options->getFilereread());
  startnreadlayout->addWidget(filerereadLabel);
  startnreadlayout->addWidget(filereread);
  nmatlabLabel = new QLabel("Database file");
  nmatlab = new QComboBox();
  nmatlab->addItems(options->nmatlabList);
  nmatlab->setCurrentIndex(options->getNmatlab());
  nmatlab->setFixedWidth(200);
  fileh5Label= new QLabel("Name h5 file");
  fileh5 = new QLineEdit(options->getH5File());
  fileh5->setFixedWidth(120);
  advancedinterfaceLabel = new QLabel("Advanced interface");
  advancedinterface = new QCheckBox(this);
  connect(advancedinterface, SIGNAL(stateChanged(int)),this,
	SLOT(advancedinterfaceCheckBoxStateChanged(int)));

  
  wavelengthLabel = new QLabel("Wavelength (nm)");
  wavelength = new QLineEdit(QString::number(options->getWavelength()));
  wavelength->setFixedWidth(120);
  
  P0Label = new QLabel("Power (W)");
  P0 = new QLineEdit(QString::number(options->getP0()));
  P0->setFixedWidth(120);

  W0Label = new QLabel("Waist (nm)");
  W0 = new QLineEdit(QString::number(options->getW0()));
  W0->setFixedWidth(120);
  
  beamLabel = new QLabel("Beam");
  beam      = new QComboBox();
  beam->addItems(options->beamList);
  beam->setCurrentIndex(beam->findText(options->getBeam()));
  beam->setFixedWidth(200);
  beam->setCurrentIndex(-1);
  connect(beam , SIGNAL(currentIndexChanged(int)),this,
	SLOT(handleBeamSelectionChanged(int))); 
  beamButton = new QPushButton("Props");
  beamButton->setFixedWidth(70);
  connect(beamButton, SIGNAL(clicked()), this, SLOT(configureBeam()));

  wavemultinumberLabel = new QLabel("Number of plane waves");
  wavemultinumber      = new QSpinBox();
  wavemultinumber->setRange(1,MAX_WAVEMULTI_NUMBER);
  wavemultinumber->setFixedWidth(120);

  objectLabel = new QLabel("Object");
  object      = new QComboBox();
  object->addItems(options->objectList);
  object->setFixedWidth(200);
  object->setCurrentIndex(-1);
  connect(object , SIGNAL(currentIndexChanged(int)),this,
	SLOT(handleObjectSelectionChanged(int)));
  objectButton = new QPushButton("Props");
  objectButton->setFixedWidth(70);
  connect(objectButton, SIGNAL(clicked()), this, SLOT(configureObject()));

  objectnumberLabel = new QLabel("Number of objects");
  objectnumber      = new QSpinBox();
  objectnumber->setRange(1,MAX_OBJECT_NUMBER);
  objectnumber->setFixedWidth(120);

  anisotropyLabel = new QLabel("Anisotropy");
  anisotropy      = new QComboBox();
  anisotropy->addItems(options->anisotropyList);
  anisotropy->setCurrentIndex(anisotropy->findText(options->getAnisotropy()));
  anisotropy->setFixedWidth(120);

  epsilonButton = new QPushButton("Epsilon");
  epsilonButton->setFixedWidth(70);
  connect(epsilonButton, SIGNAL(clicked()), this, SLOT(configureEpsilon()));

  discretizationLabel = new QLabel("Discretization");
  discretization      = new QLineEdit(QString::number(options->getDiscretization()));
  discretization->setFixedWidth(120);

  toleranceLabel = new QLabel("Iterative Method Tolerance");
  tolerance      = new QLineEdit(QString::number(options->getTolerance()));
  tolerance->setFixedWidth(120);

  methodeitLabel = new QLabel("Iterative Method");
  methodeit      = new QComboBox();
  methodeit->addItems(options->methodeitList);
  methodeit->setCurrentIndex(methodeit->findText(options->getMethodeit()));
  methodeit->setFixedWidth(120);

  polarizabilityLabel = new QLabel("Polarizability");
  polarizability      = new QComboBox();
  polarizability->addItems(options->polarizabilityList);
  polarizability->setCurrentIndex(polarizability->findText(options->getPolarizability()));
  polarizability->setFixedWidth(120);

  quadLabel = new QLabel("Integration Green function");
  quad      = new QComboBox();
  quad->addItems(options->quadList);
  quad->setCurrentIndex(quad->findText(QString::number(options->getQuad())));
  quad->setFixedWidth(120);


  nfft2dLabel = new QLabel("FFT size for far field");
  nfft2d      = new QComboBox();
  nfft2d->addItems(options->nfft2dList);
  nfft2d->setCurrentIndex(nfft2d->findText(QString::number(options->getnfft2d())));
  nfft2d->setFixedWidth(120);
  
  QGridLayout *studyfarfieldlayout = new QGridLayout();
  emptycrosssectionLabel = new QLabel(" ");
  crosssectionLabel = new QLabel("Cross section");
  crosssection = new QCheckBox(this);
  connect(crosssection, SIGNAL(stateChanged(int)),this,
	SLOT(crosssectionCheckBoxStateChanged(int)));
  emptycrosssectionpoyntingLabel = new QLabel(" ");
  crosssectionpoyntingLabel = new QLabel("Poynting + asymmetry factor");
  crosssectionpoynting = new QCheckBox(this);
  connect(crosssectionpoynting, SIGNAL(stateChanged(int)),this,
	SLOT(crosssectionpoyntingCheckBoxStateChanged(int)));
  quickdiffractLabel = new QLabel("Quick computation");
  quickdiffract = new QCheckBox(this);
  ntheta = new QLineEdit(QString::number(options->getNtheta()));
  ntheta->setFixedWidth(40);
  nthetaLabel = new QLabel("Ntheta:");
  nphi = new QLineEdit(QString::number(options->getNphi()));
  nphi->setFixedWidth(40);
  nphiLabel = new QLabel("Nphi:");
  emptynenergieLabel = new QLabel(" ");
  nenergieLabel = new QLabel("Energy conservation");
  nenergie = new QCheckBox(this);
  connect(nenergie, SIGNAL(stateChanged(int)),this,
	SLOT(nenergieCheckBoxStateChanged(int)));

  studyfarfieldlayout->addWidget(emptycrosssectionLabel,0,0,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionLabel,0,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssection,0,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(emptycrosssectionpoyntingLabel,1,0,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionpoyntingLabel,1,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(crosssectionpoynting,1,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(quickdiffractLabel,1,3,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(quickdiffract,1,4,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nthetaLabel,2,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(ntheta,2,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nphiLabel,3,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nphi,3,2,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(emptynenergieLabel,4,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nenergieLabel,4,1,Qt::AlignLeft);
  studyfarfieldlayout->addWidget(nenergie,4,2,Qt::AlignLeft);
  

 QGridLayout *studymicroscopylayout = new QGridLayout();

 
  emptymicroscopyFFTLabel = new QLabel(" ");
  microscopyFFTLabel = new QLabel("Quick computation");
  microscopyFFT = new QCheckBox(this);
   connect(microscopyFFT, SIGNAL(stateChanged(int)),this,
	SLOT(microscopyFFTCheckBoxStateChanged(int)));

  na = new QLineEdit(QString::number(options->getNA()));
  na->setFixedWidth(60);
  naLabel = new QLabel("Numerical aperture [0,1] (objective lens):");

  nainc = new QLineEdit(QString::number(options->getNAinc()));
  nainc->setFixedWidth(60);
  naincLabel = new QLabel("Numerical aperture [0,1] (condenser lens):");

  gross = new QLineEdit(QString::number(options->getGross()));
  gross->setFixedWidth(60);
  grossLabel = new QLabel("Magnification:");

  zlens = new QLineEdit(QString::number(options->getZlens()));
  zlens->setFixedWidth(60);
  zlensLabel = new QLabel("Position of the focal plane (nm):");
 
  ntypemicLabel = new QLabel("Microscope");
  ntypemic      = new QComboBox();
  ntypemic->addItems(options->ntypemicList);
  ntypemic->setCurrentIndex(options->getNtypemic());
  ntypemic->setFixedWidth(150);
  connect(ntypemic, SIGNAL(currentIndexChanged(int)),this,
	SLOT(ntypemicStateChanged(int)));

  nsideLabel = new QLabel("Side of observation");
  nside      = new QComboBox();
  nside->addItems(options->nsideList);
  nside->setCurrentIndex(options->getNside());
  nside->setFixedWidth(150);
  connect(nside, SIGNAL(currentIndexChanged(int)),this,
	SLOT(nsideStateChanged(int)));

  studymicroscopylayout->addWidget(ntypemicLabel,0,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(ntypemic,0,2,Qt::AlignLeft);  
  studymicroscopylayout->addWidget(nsideLabel,1,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(nside,1,2,Qt::AlignLeft);  
  studymicroscopylayout->addWidget(emptymicroscopyFFTLabel,2,0,Qt::AlignLeft);
  studymicroscopylayout->addWidget(microscopyFFTLabel,2,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(microscopyFFT,2,2,Qt::AlignLeft);
  studymicroscopylayout->addWidget(naLabel,3,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(na,3,2,Qt::AlignLeft);
  studymicroscopylayout->addWidget(grossLabel,4,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(gross,4,2,Qt::AlignLeft);
  studymicroscopylayout->addWidget(zlensLabel,5,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(zlens,5,2,Qt::AlignLeft);
  studymicroscopylayout->addWidget(naincLabel,6,1,Qt::AlignLeft);
  studymicroscopylayout->addWidget(nainc,6,2,Qt::AlignLeft);

  
  QGridLayout *studyforcelayout = new QGridLayout();
  emptyopticalforceLabel = new QLabel(" ");
  opticalforceLabel = new QLabel("Optical force");
  opticalforce = new QCheckBox(this);
  connect(opticalforce, SIGNAL(stateChanged(int)),this,
	SLOT(opticalforceCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticalforceLabel,0,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforceLabel,0,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforce,0,2,Qt::AlignLeft);
  emptyopticalforcedensityLabel = new QLabel(" ");
  opticalforcedensityLabel = new QLabel("Optical force density");
  opticalforcedensity = new QCheckBox(this);
  connect(opticalforcedensity, SIGNAL(stateChanged(int)),this,
	SLOT(opticalforcedensityCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticalforcedensityLabel,1,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforcedensityLabel,1,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticalforcedensity,1,2,Qt::AlignLeft);
  emptyopticaltorqueLabel = new QLabel(" ");
  opticaltorqueLabel = new QLabel("Optical torque");
  opticaltorque = new QCheckBox(this);
  connect(opticaltorque, SIGNAL(stateChanged(int)),this,
	SLOT(opticaltorqueCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticaltorqueLabel,2,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorqueLabel,2,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorque,2,2,Qt::AlignLeft);
  emptyopticaltorquedensityLabel = new QLabel(" ");
  opticaltorquedensityLabel = new QLabel("Optical torque density");
  opticaltorquedensity = new QCheckBox(this);
  connect(opticaltorquedensity, SIGNAL(stateChanged(int)),this,
	SLOT(opticaltorquedensityCheckBoxStateChanged(int)));
  studyforcelayout->addWidget(emptyopticaltorquedensityLabel,3,0,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorquedensityLabel,3,1,Qt::AlignLeft);
  studyforcelayout->addWidget(opticaltorquedensity,3,2,Qt::AlignLeft);

  QGridLayout *studynearfieldlayout = new QGridLayout();
  emptylocalfieldLabel = new QLabel(" ");
  localfieldLabel = new QLabel("Local field");
  localfield = new QCheckBox();
  connect(localfield, SIGNAL(stateChanged(int)),this,
	SLOT(localfieldCheckBoxStateChanged(int)));
  studynearfieldlayout->addWidget(emptylocalfieldLabel,0,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(localfieldLabel,0,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(localfield,0,2,Qt::AlignLeft);
  emptymacroscopicfieldLabel = new QLabel(" ");
  macroscopicfieldLabel = new QLabel("Macroscopic field");
  macroscopicfield = new QCheckBox();
  connect(macroscopicfield, SIGNAL(stateChanged(int)),this,
	SLOT(macroscopicfieldCheckBoxStateChanged(int)));
  studynearfieldlayout->addWidget(emptymacroscopicfieldLabel,1,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(macroscopicfieldLabel,1,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(macroscopicfield,1,2,Qt::AlignLeft);
  emptyrangeofstudyLabel = new QLabel(" ");
  rangeofstudyLabel = new QLabel("Range of study");
  rangeofstudy = new QComboBox();
  rangeofstudy->addItems(options->rangeofstudyList);
  rangeofstudy->setCurrentIndex(options->getNproche());
  studynearfieldlayout->addWidget(emptyrangeofstudyLabel,2,0,Qt::AlignLeft);
  studynearfieldlayout->addWidget(rangeofstudyLabel,2,1,Qt::AlignLeft);
  studynearfieldlayout->addWidget(rangeofstudy,2,2,Qt::AlignLeft);

  nxmp = new QLineEdit(QString::number(options->getNxmp()));
  nxmp->setFixedWidth(40);
  nymp = new QLineEdit(QString::number(options->getNymp()));
  nymp->setFixedWidth(40);
  nzmp = new QLineEdit(QString::number(options->getNzmp()));
  nzmp->setFixedWidth(40);
  //  nxx = new QLineEdit(QString::number(options->getNxx()));
  //    nxx->setFixedWidth(40);
  //    nyy = new QLineEdit(QString::number(options->getNyy()));
  //    nyy->setFixedWidth(40);
  //    nzz = new QLineEdit(QString::number(options->getNzz()));
  //    nzz->setFixedWidth(40);
  //    QLOG_INFO() << " UPDATE NXX = " << options->getNxx();

  QBoxLayout   *beamlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  beamlayout->addWidget(beam);
  beamlayout->addWidget(beamButton);

  QBoxLayout   *objectlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  objectlayout->addWidget(object);
  objectlayout->addWidget(objectButton);

  QBoxLayout   *epsilonlayout = new QBoxLayout(QBoxLayout::LeftToRight);
  epsilonlayout->addWidget(anisotropy);
  epsilonlayout->addWidget(epsilonButton);

  QFrame *hsep000 = new QFrame(this);
  hsep000->setFrameShape(QFrame::HLine);
  hsep000->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep000);
  layout->addRow(calculationlayout);
  layout->addRow(nrigLabel,nrig);
  layout->addRow(nreadLabel,nread);
  layout->addRow(startnreadlayout);
  layout->addRow(nmatlabLabel,nmatlab);
  layout->addRow(fileh5Label,fileh5);
  layout->addRow(advancedinterfaceLabel->text(),advancedinterface);

  QFrame *hsep0 = new QFrame(this);
  QFrame *hsep00 = new QFrame(this);
  hsep0->setFrameShape(QFrame::HLine);
  hsep0->setFrameShadow(QFrame::Sunken);
  hsep00->setFrameShape(QFrame::HLine);
  hsep00->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep0);
  
QLabel* title1 = new QLabel("Illumination properties");
title1->setFont(QFont( "Helvetica", 14, QFont::Bold, TRUE ) );
layout->addRow(title1);

  layout->addRow(hsep00);
  layout->addRow(wavelengthLabel,wavelength);
  layout->addRow(P0Label,P0);
  layout->addRow(W0Label,W0);
  layout->addRow(beamLabel,beamlayout);
  layout->addRow(wavemultinumberLabel,wavemultinumber);
  QFrame *hsep01 = new QFrame(this);
  QFrame *hsep001 = new QFrame(this);
  hsep01->setFrameShape(QFrame::HLine);
  hsep01->setFrameShadow(QFrame::Sunken);
  hsep001->setFrameShape(QFrame::HLine);
  hsep001->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep01);

QLabel* title2 = new QLabel("Object properties");
title2->setFont(QFont( "Helvetica", 14, QFont::Bold, TRUE ) );
layout->addRow(title2);

  layout->addRow(hsep001);
  layout->addRow(objectLabel,objectlayout);
  layout->addRow(objectnumberLabel,objectnumber);
  layout->addRow(anisotropyLabel,epsilonlayout);
  layout->addRow(discretizationLabel,discretization);
  QFrame *hsep02 = new QFrame(this);
  QFrame *hsep002 = new QFrame(this);
  hsep02->setFrameShape(QFrame::HLine);
  hsep02->setFrameShadow(QFrame::Sunken);
  hsep002->setFrameShape(QFrame::HLine);
  hsep002->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep02);

QLabel* title3 = new QLabel("Study");
title3->setFont(QFont( "Helvetica", 14, QFont::Bold, TRUE ) );
layout->addRow(title3);
 
  layout->addRow(hsep002);
  dipolepsilon = new QCheckBox(this);
  connect(dipolepsilon, SIGNAL(stateChanged(int)),this,
	SLOT(dipolepsilonCheckBoxStateChanged(int)));
  QLabel* subtitle1 = new QLabel("Only dipoles");
subtitle1->setFont(QFont( "Helvetica", 12, QFont::Bold, TRUE ) );
  layout->addRow(subtitle1,dipolepsilon);
  farfield = new QCheckBox(this);
  connect(farfield , SIGNAL(stateChanged(int)),this,
	SLOT(farfieldCheckBoxStateChanged(int)));
QLabel* subtitle2 = new QLabel("Far field");
subtitle2->setFont(QFont( "Helvetica", 12, QFont::Bold, TRUE ) );
  layout->addRow(subtitle2,farfield);
  layout->addRow(studyfarfieldlayout);

  microscopy = new QCheckBox(this);
  connect(microscopy, SIGNAL(stateChanged(int)),this,
	SLOT(microscopyCheckBoxStateChanged(int)));

  
QLabel* subtitle5 = new QLabel("Microscopy");
subtitle5->setFont(QFont( "Helvetica", 12, QFont::Bold, TRUE ) );
  layout->addRow(subtitle5,microscopy);
  layout->addRow(studymicroscopylayout);

  
  force = new QCheckBox(this);
  connect(force , SIGNAL(stateChanged(int)),this,
	SLOT(forceCheckBoxStateChanged(int)));
QLabel* subtitle3 = new QLabel("Optical Force");
subtitle3->setFont(QFont( "Helvetica", 12, QFont::Bold, TRUE ) );
  layout->addRow(subtitle3,force);
  layout->addRow(studyforcelayout);

  nearfield = new QCheckBox(this);
  connect(nearfield , SIGNAL(stateChanged(int)),this,
	SLOT(nearfieldCheckBoxStateChanged(int)));
  QLabel* subtitle4 = new QLabel("Near field");
subtitle4->setFont(QFont( "Helvetica", 12, QFont::Bold, TRUE ) );
  layout->addRow(subtitle4,nearfield);
  layout->addRow(studynearfieldlayout);
  QFrame *hsep1 = new QFrame(this);
  QFrame *hsep2 = new QFrame(this);
  hsep1->setFrameShape(QFrame::HLine);
  hsep1->setFrameShadow(QFrame::Sunken);
  hsep2->setFrameShape(QFrame::HLine);
  hsep2->setFrameShadow(QFrame::Sunken);
  layout->addRow(hsep1);

  
QLabel* title4 = new QLabel("Numerical parameters");
title4->setFont(QFont( "Helvetica", 14, QFont::Bold, TRUE ) );
layout->addRow(title4);
 
  layout->addRow(hsep2);
  layout->addRow(toleranceLabel,tolerance);
  layout->addRow(methodeitLabel,methodeit);
  layout->addRow("Additional sideband x",nxmp);
  layout->addRow("Additional sideband y",nymp);
  layout->addRow("Additional sideband z",nzmp);
  layout->addRow(polarizabilityLabel,polarizability);
  layout->addRow(quadLabel,quad);
  layout->addRow(nfft2dLabel,nfft2d);
  this->setLayout(layout);
  update();
  beamconfigdlg = NULL;
  objectconfigdlg = NULL;
  epsilonconfigdlg = NULL;

  // init now object combo index (triggers the handle)
  QLOG_DEBUG ( ) << "Number of objects : " << options->getObjectNumber();
  objectnumber->setValue(options->getObjectNumber());
  object->setCurrentIndex(object->findText(options->getObject()));
  connect(objectnumber , SIGNAL(valueChanged(int)),this,
	SLOT(handleObjectNumberSelectionChanged(int)));

  // init now wavemulti combo index (triggers the handle)
  QLOG_DEBUG ( ) << "Number of plane waves : " << options->getWaveMultiNumber();
  wavemultinumber->setValue(options->getWaveMultiNumber());
  beam->setCurrentIndex(beam->findText(options->getBeam()));
  connect(wavemultinumber , SIGNAL(valueChanged(int)),this,
	SLOT(handleWaveMultiNumberSelectionChanged(int)));
  
}
OptionsWidget::~OptionsWidget()
{
  QLOG_DEBUG ( ) << "Deleting OptionsWidget"; 
}
void 
OptionsWidget::execute() {
  RunWidget *runwidget = NULL;
  OptionsWidget *optionswidget = NULL;
  if (mainwindow) {
    QTabWidget *tabwidget = (QTabWidget*) mainwindow->findChild<QTabWidget*>("TabWidget");
    QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
    runwidget = currentWindow->findChild<RunWidget *>("Run");
    optionswidget = currentWindow->findChild<OptionsWidget *>("Options");
    optionswidget->updateOptions();
    runwidget->execute();
   }
}
void
OptionsWidget::saveName() {
  QLOG_INFO() << "PRESENT SAVE NAME " << name;
  this->updateOptions();
  if (name == "new") {
    nameLineEdit->clear();
    descriptionLineEdit->clear();
    cfgWindow->exec();
    QLOG_INFO() << "NEW SAVE NAME " << name;
  }
  else {
   nameLineEdit->setText(name);
   this->finish();
  }
}
void
OptionsWidget::saveAsName() {
  QLOG_INFO() << "PRESENT SAVE NAME " << name;
  this->updateOptions();
  nameLineEdit->clear();
  descriptionLineEdit->clear();
  cfgWindow->exec();
  QLOG_INFO() << "NEW SAVE NAME " << name;
}
void 
OptionsWidget::finish() {
 name = nameLineEdit->text();
 QLOG_INFO() << "OptionsWidget::finish() " << name;
 if ( name == "" ) {
   QMessageBox::warning(this, "Warning:", "Enter a valid name");
   return;
 }
 options->setName(name);
 options->setDescription(descriptionLineEdit->text());
 description = options->getDescription();
 options->saveDb(name,description);
 cfgWindow->hide();
 if (mainwindow) {
   QTabWidget *tabwidget = (QTabWidget*) mainwindow->findChild<QTabWidget*>("TabWidget");
   QMainWindow *currentWindow = (QMainWindow*)tabwidget->currentWidget();
   QDockWidget *dockWidget = currentWindow->findChild<QDockWidget*>("DockWidget");
   RunWidget *runwidget = currentWindow->findChild<RunWidget *>("Run");
   QDockWidget *dockWidgetCentral = runwidget->findChild<QDockWidget*>("DockWidgetCentral");
   QDockWidget *dockWidgetOutput = runwidget->findChild<QDockWidget*>("DockWidgetOutput");
   QDockWidget *dockWidgetPlot = runwidget->findChild<QDockWidget*>("DockWidgetPlot");
   dockWidget->setWindowTitle(name);
   dockWidgetCentral->setWindowTitle(name);
   dockWidgetOutput->setWindowTitle(name);
   dockWidgetPlot->setWindowTitle(name);
   tabwidget->setTabText(tabwidget->currentIndex(),name);
 }
  emit updateOptionsWindow();
}
void 
OptionsWidget::cancel() {
 cfgWindow->hide();
}
void 
OptionsWidget::handleObjectSelectionChanged(int index){

 QLOG_DEBUG() << "OptionsWidget::handleObjectSelectionChanged";
 if ( object->currentText() != "multiple spheres" && 
      object->currentText() != "concentric spheres")  {
   this->setObjectNumber(1);
   objectnumber->setEnabled( false );
 }
 else {
   objectnumber->setValue(options->getObjectNumber());
   objectnumber->setEnabled( true );
 }
}
void 
OptionsWidget::handleObjectNumberSelectionChanged(int index) {
//if ( object->currentText() == "nspheres" ||
//     object->currentText() == "concentricsphere")
   options->setObjectNumber(this->getObjectNumber());
}
void 
OptionsWidget::handleBeamSelectionChanged(int index){

 QLOG_DEBUG() << "OptionsWidget::handleBeamSelectionChanged";
 if ( beam->currentText() != "Multiplane wave" )  {
   this->setWaveMultiNumber(1);
   wavemultinumber->setEnabled( false );
 }
 else {
   wavemultinumber->setValue(options->getWaveMultiNumber());
   wavemultinumber->setEnabled( true );
 }
}
void 
OptionsWidget::handleWaveMultiNumberSelectionChanged(int index) {
//if ( object->currentText() == "nspheres" ||
//     object->currentText() == "concentricsphere")
   options->setWaveMultiNumber(this->getWaveMultiNumber());
}
void 
OptionsWidget::nreadCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   filerereadLabel->show();
   filereread->show();
  }
  else if (state == Qt::Unchecked) {
   filerereadLabel->hide();
   filereread->hide(); 
  }
}
void
OptionsWidget::advancedinterfaceCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
    methodeitLabel->show();
    methodeit->show();
    toleranceLabel->show();
    tolerance->show();
    polarizabilityLabel->show();
    polarizability->show();
    quadLabel->show();
    quad->show();
    nfft2dLabel->show();
    nfft2d->show();
    nrigLabel->show();
    nrig->show();
    nreadLabel->show();
    nread->show();
    nmatlabLabel->show();
    nmatlab->show();
    fileh5Label->show();
    fileh5->show();
  }
   else if (state == Qt::Unchecked) {
     methodeitLabel->hide();
     methodeit->hide();
     toleranceLabel->hide();
     tolerance->hide();
     polarizabilityLabel->hide();
     polarizability->hide();
     quadLabel->hide();
     quad->hide();
     nfft2dLabel->hide();
     nfft2d->hide();
     nrigLabel->hide();
     nrig->hide();
     nreadLabel->hide();
     nread->hide();
     nmatlabLabel->hide();
     nmatlab->hide();
     fileh5Label->hide();
     fileh5->hide();
   }
}
void 
OptionsWidget::dipolepsilonCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   farfield->setChecked(false);
   microscopy->setChecked(false);
   nearfield->setChecked(false);
   force->setChecked(false);
  }
}
void 
OptionsWidget::farfieldCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptycrosssectionLabel->show();
   crosssectionLabel->show();
   crosssection->show();
   emptycrosssectionpoyntingLabel->show();
   crosssectionpoyntingLabel->show();
   crosssectionpoynting->show();
   emptynenergieLabel->show();
   nenergieLabel->show();
   nenergie->show(); 
   dipolepsilon->setChecked(false);
  }
  else if (state == Qt::Unchecked) {
   emptycrosssectionLabel->hide();
   crosssection->setChecked(false);
   crosssectionLabel->hide();
   crosssection->hide();  
   emptycrosssectionpoyntingLabel->hide();
   crosssectionpoynting->setChecked(false);
   crosssectionpoyntingLabel->hide();
   crosssectionpoynting->hide();
   crosssectionpoyntingCheckBoxStateChanged(Qt::Unchecked);
   emptynenergieLabel->hide(); 
   nenergie->setChecked(false);
   nenergieLabel->hide();
   nenergie->hide();
   nenergieCheckBoxStateChanged(Qt::Unchecked);
  }
}
void 
OptionsWidget::crosssectionCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::nenergieCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::crosssectionpoyntingCheckBoxStateChanged(int state) {
QLOG_DEBUG() << "OptionsWidget::crosssectionpoyntingCheckBoxStateChanged" << state;
if (state == Qt::Checked) {
   quickdiffractLabel->show();
   quickdiffract->show();
   nthetaLabel->show();
   ntheta->show();
   nphiLabel->show();
   nphi->show();
  }
  else if (state == Qt::Unchecked) {
   quickdiffractLabel->hide();
   quickdiffract->setChecked(false);
   quickdiffract->hide();
   nthetaLabel->hide();
   ntheta->hide();
   nphiLabel->hide();
   nphi->hide();
  }
}
void 
OptionsWidget::microscopyCheckBoxStateChanged(int state) {
if (state == Qt::Checked) {
   microscopyFFT->setChecked(true);
   microscopyFFTLabel->show();
   microscopyFFT->show();
   naLabel->show();
   na->show();
   grossLabel->show();
   gross->show();
   zlensLabel->show();
   zlens->show();
   ntypemicLabel->show();
   ntypemic->show();
   nsideLabel->show();
   nside->show();
   naincLabel->show();
   nainc->show();
   dipolepsilon->setChecked(false);
 }
  else if (state == Qt::Unchecked) {
   microscopyFFTLabel->hide();
   microscopyFFT->setChecked(false);
   microscopyFFT->hide();
   naLabel->hide();
   na->hide();   
   grossLabel->hide();
   gross->hide();
   zlensLabel->hide();
   zlens->hide();
   naincLabel->hide();
   nainc->hide();
   ntypemicLabel->hide();
   ntypemic->hide();
   nsideLabel->hide();
   nside->hide();
  }
}
void 
OptionsWidget::ntypemicStateChanged(int state) {
QLOG_INFO() << ntypemic->currentText() ;
if ( ntypemic->currentText() == "Holographic" )  {
  beam->setEnabled( true );
  beamButton->setEnabled( true );
 }
 else {
   beam->setEnabled( false );
   beamButton->setEnabled( false );
}  
}
void
OptionsWidget::nsideStateChanged(int state) {
}
void 
OptionsWidget::microscopyFFTCheckBoxStateChanged(int state) {
}
void 
OptionsWidget::forceCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptyopticalforceLabel->show();
   opticalforceLabel->show();
   opticalforce->show();
   emptyopticalforcedensityLabel->show();
   opticalforcedensityLabel->show();
   opticalforcedensity->show();
   emptyopticaltorqueLabel->show();
   opticaltorqueLabel->show();
   opticaltorque->show();
   emptyopticaltorquedensityLabel->hide();
   opticaltorquedensityLabel->show();
   opticaltorquedensity->show();
   dipolepsilon->setChecked(false);
  }
  else if (state == Qt::Unchecked) {
   emptyopticalforceLabel->hide();
   opticalforce->setChecked(false);
   opticalforceLabel->hide();
   opticalforce->hide();
   emptyopticalforcedensityLabel->hide();
   opticalforcedensity->setChecked(false);
   opticalforcedensityLabel->hide();
   opticalforcedensity->hide();
   emptyopticaltorqueLabel->hide();
   opticaltorque->setChecked(false);
   opticaltorqueLabel->hide();
   opticaltorque->hide();
   emptyopticaltorquedensityLabel->hide();
   opticaltorquedensity->setChecked(false);
   opticaltorquedensityLabel->hide();
   opticaltorquedensity->hide();
  }
}
void 
OptionsWidget::opticalforceCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    opticalforcedensity->setChecked(false);
    opticaltorque->setChecked(false);
    if ( opticaltorque->isChecked() == false &&
        opticaltorquedensity->isChecked() == false )
     force->setChecked(false);
   }
}
void 
OptionsWidget::opticalforcedensityCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    if (opticalforce->isChecked() == false &&
        opticaltorque->isChecked() == false &&
        opticaltorquedensity->isChecked() == false )
     force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticalforce->setChecked(true);
   
}
void 
OptionsWidget::opticaltorqueCheckBoxStateChanged(int state) {
    if (state == Qt::Unchecked) {
      opticaltorquedensity->setChecked(false);
      if ( opticalforcedensity->isChecked() == false &&
        opticalforce->isChecked() == false )
     force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticalforce->setChecked(true);
}
void 
OptionsWidget::opticaltorquedensityCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked) {
    if (opticalforcedensity->isChecked() == false &&
        opticaltorque->isChecked() == false &&
        opticalforce->isChecked() == false )
      force->setChecked(false);
   }
   else if (state == Qt::Checked)
     opticaltorque->setChecked(true);
}
void 
OptionsWidget::nearfieldCheckBoxStateChanged(int state) {
  if (state == Qt::Checked) {
   emptyrangeofstudyLabel->show();
   rangeofstudyLabel->show();
   rangeofstudy->show();
   emptylocalfieldLabel->show();
   localfieldLabel->show();
   localfield->show();
   emptymacroscopicfieldLabel->show();
   macroscopicfieldLabel->show();
   macroscopicfield->show();
   dipolepsilon->setChecked(false);
   nxmp->setText(QString::number(0));
   nymp->setText(QString::number(0));
   nzmp->setText(QString::number(0));
   nxmp->show();
   nymp->show();
   nzmp->show();
  }
  else if (state == Qt::Unchecked) {
   emptylocalfieldLabel->hide();
   localfield->setChecked(false);
   localfield->hide();  
   localfieldLabel->hide();
   emptymacroscopicfieldLabel->hide();
   macroscopicfield->setChecked(false);
   macroscopicfield->hide();
   macroscopicfieldLabel->hide();
   emptyrangeofstudyLabel->hide();
   rangeofstudy->hide();  
   rangeofstudyLabel->hide();
   nxmp->setText(QString::number(0));
   nymp->setText(QString::number(0));
   nzmp->setText(QString::number(0));
   nxmp->hide();
   nymp->hide();
   nzmp->hide();
  }
}
void 
OptionsWidget::localfieldCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked)
    if (macroscopicfield->isChecked() == false )
     nearfield->setChecked(false);
}
void 
OptionsWidget::macroscopicfieldCheckBoxStateChanged(int state) {
   if (state == Qt::Unchecked)
    if (localfield->isChecked() == false )
     nearfield->setChecked(false);
}
QString 
OptionsWidget::getFilereread(){
  return filereread->text();
}
QString 
OptionsWidget::getH5File(){
  return fileh5->text();
}
double 
OptionsWidget::getWavelength(){
  return wavelength->text().toDouble();
}
double 
OptionsWidget::getP0(){
  return P0->text().toDouble();
}
double 
OptionsWidget::getW0(){
  return W0->text().toDouble();
}
QString 
OptionsWidget::getBeam(){
  return beam->currentText();
}
QString 
OptionsWidget::getObject(){
  return object->currentText();
}
int 
OptionsWidget::getObjectNumber(){
  return objectnumber->value();
}
int 
OptionsWidget::getWaveMultiNumber(){
  return wavemultinumber->value();
}
QString 
OptionsWidget::getAnisotropy(){
  return anisotropy->currentText();
}

int 
OptionsWidget::getDiscretization(){
  return discretization->text().toInt();
}
double 
OptionsWidget::getTolerance(){
  return tolerance->text().toDouble();
}
QString 
OptionsWidget::getMethodeit(){
  return methodeit->currentText();
}
QString 
OptionsWidget::getPolarizability(){
  return polarizability->currentText();
}
int 
OptionsWidget::getQuad(){
  return quad->currentText().toInt();
}
int 
OptionsWidget::getnfft2d(){
  return nfft2d->currentText().toInt();
}
void 
OptionsWidget::setWavelength(double _wavelength){
  wavelength->setText(QString::number(_wavelength));
}
void 
OptionsWidget::setP0(double _P0){
  P0->setText(QString::number(_P0));
}
void 
OptionsWidget::setW0(double _W0){
  W0->setText(QString::number(_W0));
}
void 
OptionsWidget::setBeam(QString _beam){
  beam->setCurrentIndex(beam->findText(_beam));
}
void 
OptionsWidget::setObject(QString _object){
  object->setCurrentIndex(object->findText(_object));
}
void 
OptionsWidget::setObjectNumber(int _objectnumber){
  objectnumber->setValue(_objectnumber);
}
void 
OptionsWidget::setWaveMultiNumber(int _wavemultinumber){
  wavemultinumber->setValue(_wavemultinumber);
}
void 
OptionsWidget::setAnisotropy(QString _anisotropy){
anisotropy->setCurrentIndex(anisotropy->findText(_anisotropy));
}
void 
OptionsWidget::setDiscretization(int _discretization){
discretization->setText(QString::number(_discretization));
}
void 
OptionsWidget::setTolerance(double _tolerance){
tolerance->setText(QString::number(_tolerance));
}
void 
OptionsWidget::setMethodeit(QString _methodeit){
methodeit->setCurrentIndex(methodeit->findText(_methodeit));
}
void 
OptionsWidget::setPolarizability(QString _polarizability){
polarizability->setCurrentIndex(polarizability->findText(_polarizability));
}
void 
OptionsWidget::setQuad(int _quad){
quad->setCurrentIndex(quad->findText(QString::number(_quad)));
}
void 
OptionsWidget::setnfft2d(int _nfft2d){
nfft2d->setCurrentIndex(nfft2d->findText(QString::number(_nfft2d)));
}
void OptionsWidget::configureBeam() {
   this->updateOptions();
   if (beamconfigdlg) delete beamconfigdlg;
   beamconfigdlg = new BeamConfigDialog(this,Qt::Widget,options);
   beamconfigdlg->exec();
}
void OptionsWidget::configureObject(){
   this->updateOptions();
   if (objectconfigdlg) delete objectconfigdlg;
   objectconfigdlg = new ObjectConfigDialog(this,Qt::Widget,options);
   objectconfigdlg->exec();
}
void OptionsWidget::configureEpsilon() {
   this->updateOptions();
   if (epsilonconfigdlg) delete epsilonconfigdlg;
   epsilonconfigdlg = new EpsilonConfigDialog(this,Qt::Widget,options);
   epsilonconfigdlg->exec();
}
void 
OptionsWidget::updateFarfield() {
  crosssection->setChecked(options->getCrosssection());
  crosssectionpoynting->setChecked(options->getCrosssectionpoynting());
  quickdiffract->setChecked(options->getQuickdiffract());
  nenergie->setChecked(options->getNenergie());
  ntheta->setText(QString::number(options->getNtheta()));
  nphi->setText(QString::number(options->getNphi()));
  if (crosssection->isChecked() || crosssectionpoynting->isChecked() ||
      nenergie->isChecked())
   farfield->setChecked(true);
  farfieldCheckBoxStateChanged(farfield->isChecked());
}
void 
OptionsWidget::updateMicroscopy() {
  microscopy->setChecked(options->getMicroscopy());
  microscopyFFT->setChecked(options->getMicroscopyFFT());
  na->setText(QString::number(options->getNA()));  
  gross->setText(QString::number(options->getGross()));
  zlens->setText(QString::number(options->getZlens()));
  nainc->setText(QString::number(options->getNAinc()));
  if ( microscopy->isChecked() )
  microscopy->setChecked(true);
  microscopyCheckBoxStateChanged(microscopy->isChecked());
  }
void 
OptionsWidget::updateNearfield() {
  localfield->setChecked(options->getLocalfield());
  macroscopicfield->setChecked(options->getMacroscopicfield());
  rangeofstudy->setCurrentIndex(options->getNproche());
  if (localfield->isChecked() || macroscopicfield->isChecked())
   nearfield->setChecked(true);
  nearfieldCheckBoxStateChanged(nearfield->isChecked()); 
}
void 
OptionsWidget::updateForce() {
  opticalforce->setChecked(options->getOpticalforce());
  opticalforcedensity->setChecked(options->getOpticalforcedensity());
  opticaltorque->setChecked(options->getOpticaltorque());
  opticaltorquedensity->setChecked(options->getOpticaltorquedensity());
  if (opticalforce->isChecked() || opticalforcedensity->isChecked() ||
      opticaltorque->isChecked() || opticaltorquedensity->isChecked())
   force->setChecked(true);
  forceCheckBoxStateChanged(force->isChecked());
}
void 
OptionsWidget::updateAdvancedinterface() {
  //  methodeit->setChecked(options->getMethodeit());
}
void 
OptionsWidget::update() {
  filereread->setText(options->getFilereread());
  fileh5->setText(options->getH5File());
  this->setWavelength(options->getWavelength());
  this->setP0(options->getP0());
  this->setW0(options->getW0());
  this->setBeam(options->getBeam());
  this->setObject(options->getObject());
  this->setObjectNumber(options->getObjectNumber());
  this->setWaveMultiNumber(options->getWaveMultiNumber());
  this->setAnisotropy(options->getAnisotropy());
  this->setDiscretization(options->getDiscretization());
  this->setTolerance(options->getTolerance());
  this->setMethodeit(options->getMethodeit());
  this->setPolarizability(options->getPolarizability());
  this->setQuad(options->getQuad());
  this->setnfft2d(options->getnfft2d());

  nread->setChecked(options->getNread());
  nreadCheckBoxStateChanged(nread->isChecked());
  //  nmatlab->setChecked(options->getNmatlab());
  advancedinterface->setChecked(options->getAdvancedinterface());
  advancedinterfaceCheckBoxStateChanged(advancedinterface->isChecked());  
  dipolepsilon->setChecked(options->getDipolepsilon());
  updateFarfield();
  updateAdvancedinterface();
  updateForce();
  updateNearfield();
  updateMicroscopy();
  nxmp->setText(QString::number(options->getNxmp()));
  nymp->setText(QString::number(options->getNymp()));
  nzmp->setText(QString::number(options->getNzmp()));
  //    nxx->setText(QString::number(options->getNxx()));
  //    nyy->setText(QString::number(options->getNyy()));
  //    nzz->setText(QString::number(options->getNzz()));
  //    QLOG_INFO() << " UPDATE NXX = " << options->getNxx();
}
void 
OptionsWidget::updateOptions() {

  options->setNread(nread->isChecked());
  options->setNmatlab(nmatlab->currentIndex());
   options->setH5File(this->getH5File());
  //  options->setNmatlab(nmatlab->isChecked());
  options->setAdvancedinterface(advancedinterface->isChecked());  
  options->setWavelength(this->getWavelength());
  options->setP0(this->getP0());
  options->setW0(this->getW0());
  options->setBeam(this->getBeam());
  options->setObject(this->getObject());
  options->setObjectNumber(this->getObjectNumber());
  options->setWaveMultiNumber(this->getWaveMultiNumber());
  options->setAnisotropy(this->getAnisotropy());
  options->setDiscretization(this->getDiscretization());
  options->setTolerance(this->getTolerance());
  options->setMethodeit(this->getMethodeit());
  options->setPolarizability(this->getPolarizability());
  options->setQuad(this->getQuad());
  options->setDipolepsilon(dipolepsilon->isChecked());
  options->setFarfield(farfield->isChecked());
  options->setNearfield(nearfield->isChecked());
  options->setForce(force->isChecked());
  options->setLocalfield(localfield->isChecked());
  options->setMacroscopicfield(macroscopicfield->isChecked());
  options->setCrosssection(crosssection->isChecked());
  options->setCrosssectionpoynting(crosssectionpoynting->isChecked());
  options->setQuickdiffract(quickdiffract->isChecked());
  options->setNrig(nrig->currentIndex());
  options->setNenergie(nenergie->isChecked());
  options->setMicroscopy(microscopy->isChecked());
  options->setMicroscopyFFT(microscopyFFT->isChecked());
  options->setNtypemic(ntypemic->currentIndex());  
  options->setNside(nside->currentIndex());  
  options->setNenergie(nenergie->isChecked());
  options->setOpticalforce(opticalforce->isChecked());
  options->setOpticalforcedensity(opticalforcedensity->isChecked());
  options->setOpticaltorque(opticaltorque->isChecked());
  options->setOpticaltorquedensity(opticaltorquedensity->isChecked());
  options->setNproche(rangeofstudy->currentIndex());

  if ( options->getNread() == true )
   options->setFilereread(this->getFilereread());
  if ( options->getNearfield() == true ) {
    options->setNxm(this->getDiscretization()+2*nxmp->text().toInt());
    options->setNym(this->getDiscretization()+2*nymp->text().toInt());
    options->setNzm(this->getDiscretization()+2*nzmp->text().toInt());
    options->setNxmp(nxmp->text().toInt());
    options->setNymp(nymp->text().toInt());
    options->setNzmp(nzmp->text().toInt());
  //      options->setNxx(nxx->text().toInt());
  //      options->setNyy(nyy->text().toInt());
   //     options->setNzz(nzz->text().toInt());
  //      QLOG_INFO() << " UPDATE NXX = " << options->getNxx();
  }
  else {
    options->setNxm(this->getDiscretization());
    options->setNym(this->getDiscretization());
    options->setNzm(this->getDiscretization());
    options->setNxmp(nxmp->text().toInt());
    options->setNymp(nymp->text().toInt());
    options->setNzmp(nzmp->text().toInt());
  //      options->setNxx(nxx->text().toInt());
  //      options->setNyy(nyy->text().toInt());
   //     options->setNzz(nzz->text().toInt());
   //     QLOG_INFO() << " UPDATE NXX = " << options->getNxx();
  }
  options->setNtheta(ntheta->text().toInt());
  options->setNphi(nphi->text().toInt());
  options->setNA(na->text().toDouble());
  options->setGross(gross->text().toDouble());
  options->setZlens(zlens->text().toDouble());  
  options->setNAinc(nainc->text().toDouble());
  options->setnfft2d(this->getnfft2d());  
  if ( options->getLocalfield() == false && options->getMacroscopicfield() == false)
    options->setNearfield(false);
  if ( options->getOpticalforce() == false && options->getOpticalforcedensity() == false &&
        options->getOpticaltorque() == false && options->getOpticaltorquedensity() == false )
    options->setForce(false);
  if ( options->getNenergie() == false && 
       options->getCrosssection() == false &&
       options->getCrosssectionpoynting() == false)
    options->setFarfield(false);
  if ( options->getMicroscopyFFT() == false ) 
  QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getGross());
    QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getZlens());
  QLOG_DEBUG() << "OptionsWidget::updateOptions> " << QString::number(options->getNrig());
  QLOG_DEBUG() << "OptionsWidget::updateOptions> crosssection " << crosssection->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> crosssectionpoynting " << crosssectionpoynting->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> quickdiffract " << quickdiffract->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> localfield " << localfield->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> macroscopicfield " << macroscopicfield->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticalforce " << opticalforce->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticalforcedensity " << opticalforcedensity->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticaltorque " << opticaltorque->isChecked();
  QLOG_DEBUG() << "OptionsWidget::updateOptions> opticaltorquedensity " << opticaltorquedensity->isChecked();
}

