#include "cdmRunWidget.h"

RunWidget::RunWidget(Options *_options, Run *_run)
{
   options = _options;
   run = _run;
   ntheta = options->getNtheta();
   nphi = options->getNphi();
   //   nxx = options->getNxx();
   //     nyy = options->getNyy();
    //    nzz = options->getNzz();
   nxm = options->getNxm();
   nym = options->getNym();
   nzm = options->getNzm();   
   nxmp = options->getNxmp();
   nymp = options->getNymp();
   nzmp = options->getNzmp();
   stopFlag = NULL;
   infoMessage = NULL;
   canceldlg = NULL;
   progressdlg = NULL;
   this->setObjectName("Run");
   this->setDockNestingEnabled(true);
   centralwidget = new QWidget(this);
   dockWidgetCentral = 	new QDockWidget(options->getName(), this);
   dockWidgetCentral->setObjectName("DockWidgetCentral");
   dockWidgetCentral->setFeatures(QDockWidget::DockWidgetMovable |
                                 QDockWidget::DockWidgetFloatable);
   connect(dockWidgetCentral, SIGNAL(topLevelChanged(bool)), this, 
				SLOT(dockWidgetCentral_topLevelChanged(bool)));
   dockWidgetCentral->setFeatures(QDockWidget::DockWidgetMovable |
                                 QDockWidget::DockWidgetFloatable);
   centralArea = new QScrollArea(this);
   centralArea->setWidget(centralwidget);
   centralArea->setWidgetResizable(true);
   dockWidgetCentral->setWidget(centralArea);
   this->addDockWidget(Qt::TopDockWidgetArea, dockWidgetCentral);
   outputwidget = new QWidget(this);
   dockWidgetOutput = new QDockWidget(options->getName(), this);
   dockWidgetOutput->setObjectName("DockWidgetOutput");
   connect(dockWidgetOutput, SIGNAL(topLevelChanged(bool)), this, 
				SLOT(dockWidgetOutput_topLevelChanged(bool)));
   dockWidgetOutput->setFeatures(QDockWidget::DockWidgetMovable |
                                 QDockWidget::DockWidgetFloatable);
   outputArea = new QScrollArea(this);
   outputArea->setWidget(outputwidget);
   outputArea->setWidgetResizable(true);
   dockWidgetOutput->setWidget(outputArea);
   this->addDockWidget(Qt::BottomDockWidgetArea, dockWidgetOutput);
   plotwidget = new QTabWidget(this);
   plotwidget->setMovable(true);
   plotwidget->setTabsClosable(true);
   connect(plotwidget, SIGNAL(tabCloseRequested ( int )), this, SLOT(closeWidgetTab(int)));
   dockWidgetPlot = new QDockWidget(options->getName(),this);
   dockWidgetPlot->setObjectName("DockWidgetPlot");
   connect(dockWidgetPlot, SIGNAL(topLevelChanged(bool)), this, 
                SLOT(dockWidgetPlot_topLevelChanged(bool)));
   dockWidgetPlot->setFeatures(QDockWidget::DockWidgetMovable |
                               QDockWidget::DockWidgetFloatable);
   plotArea = new QScrollArea(this);
   plotArea->setWidget(plotwidget);
   plotArea->setWidgetResizable(true);
   dockWidgetPlot->setWidget(plotArea);
   this->addDockWidget(Qt::BottomDockWidgetArea, dockWidgetPlot);
}
RunWidget::~RunWidget()
{
  delete run;
  QLOG_DEBUG ( ) << "Deleting RunWidget"; 
}
void 
RunWidget::closeWidgetTab(int index)
{
  QWidget * tmpwidget = plotwidget->widget(index);
  Plot *plot = NULL;
  plot = tmpwidget->findChild<Plot *>("plot3d");
  if ( plot != NULL )
     delete plot;
  PlotRaster *raster = NULL;
  raster = tmpwidget->findChild<PlotRaster *>("plotraster");
  if ( raster != NULL )
     delete raster;
  PlotVector *plotvect = NULL;
  plotvect = tmpwidget->findChild<PlotVector *>("plotvector");
  if ( plotvect != NULL )
     delete plotvect;
  plotwidget->removeTab(index);
  delete tmpwidget;
}
void RunWidget::
dockWidgetCentral_topLevelChanged(bool isFloating)  {
 if (isFloating) {
   dockWidgetCentral->setWindowFlags(Qt::Window);
   dockWidgetCentral->setWindowFlags(dockWidgetCentral->windowFlags() 
				& ~Qt::WindowCloseButtonHint);
   dockWidgetCentral->show();
 }
}
void RunWidget::
dockWidgetOutput_topLevelChanged(bool isFloating)  {
 if (isFloating) {
   
   dockWidgetOutput->setWindowFlags(Qt::Window);
   dockWidgetOutput->setWindowFlags(dockWidgetOutput->windowFlags() 
				& ~Qt::WindowCloseButtonHint);
   dockWidgetOutput->show();
 }
}
void RunWidget::
dockWidgetPlot_topLevelChanged(bool isFloating)  {
 if (isFloating) {
   dockWidgetPlot->setWindowFlags(Qt::Window);
   dockWidgetPlot->setWindowFlags(dockWidgetPlot->windowFlags()
                & ~Qt::WindowCloseButtonHint);
   dockWidgetPlot->show();
 }
}
void 
RunWidget::execute() {
   
   // Clean up widgets
   delete centralwidget->layout();
   qDeleteAll(centralwidget->children());
   centrallayout = new QFormLayout();
   centralwidget->setLayout(centrallayout);

   delete outputwidget->layout();
   qDeleteAll(outputwidget->children());
   outputlayout = new QFormLayout();
   outputwidget->setLayout(outputlayout);
   this->cleanupPlots();

   if (infoMessage) {delete infoMessage; infoMessage = NULL;}
   infoMessage = new QString("");
   if (stopFlag) {free(stopFlag); stopFlag = NULL;}
   stopFlag = (int*) malloc(sizeof(int));
   future = new QFuture<void>;
   watcher = new QFutureWatcher<void>;
   progressdlg = new QProgressDialog("Calculation progress", "Cancel", 0, 0, this,Qt::Widget);
   progressdlg->setWindowModality(Qt::WindowModal);
   connect(progressdlg, SIGNAL(canceled()),this, SLOT(cancel_thread()));
   connect(watcher, SIGNAL(finished()),this, SLOT(displayFinishedBox()));
   progressdlg->show();
   run_thread();
}
void 
RunWidget::cancel_thread()
{
   QLOG_DEBUG () << "Cancelling thread...";
   canceldlg = new QProgressDialog("Cancelling calculation...", "Cancel",0, 0, this,Qt::Widget);
   canceldlg->setCancelButton(0);
   canceldlg->setWindowModality(Qt::WindowModal);
   canceldlg->show();
   
   // Cancel thread
   *stopFlag = -1;
}
void 
RunWidget::run_thread()
{
   // Launch fortran lib calculation
   *future = QtConcurrent::run(cdmlibwrapper,options, run, infoMessage, stopFlag);
    watcher->setFuture(*future);
}
void cdmlibwrapper(Options *options, Run *run, QString *infoMessage, int *stopFlag) 
{
    QLOG_DEBUG () << " Launching calculations with following parameters:";
    // cdm.in file
    double wavelength;
    wavelength = options->getWavelength();
    QLOG_DEBUG () << "Wavelength:" << QString::number(wavelength,'g',5);
    char beam[64] = {' '};
    if (options->getBeam() == "Circular plane wave") strcpy (beam,"pwavecircular");
    if (options->getBeam() == "Linear plane wave") strcpy (beam,"pwavelinear");
    if (options->getBeam() == "Multiplane wave") strcpy (beam,"wavelinearmulti");
    if (options->getBeam() == "Circular Gaussian") strcpy (beam,"gwavecircular");
    if (options->getBeam() == "Linear Gaussian") strcpy (beam,"gwavelinear");
    if (options->getBeam() == "Circular Gaussian (FFT)") strcpy (beam,"gfftwavecircular");
    if (options->getBeam() == "Linear Gaussian (FFT)") strcpy (beam,"gfftwavelinear");
    if (options->getBeam() == "Circular Gaussian (para)") strcpy (beam,"gparawavecircular");
    if (options->getBeam() == "Linear Gaussian (para)") strcpy (beam,"gparawavelinear");
    if (options->getBeam() == "Antenna") strcpy (beam,"antenna");
    if (options->getBeam() == "Arbitrary wave (file)") strcpy (beam,"arbitrary");

    QLOG_DEBUG () << "Beam:" << options->getBeam();
    char namefileinc[64];
    for(int i = 0; i < 64; i++)
     namefileinc[i] = ' ';
    strncpy(namefileinc,(char*)options->getBeamFile().toStdString().c_str(),options->getBeamFile().size());
    QLOG_DEBUG () << "BeamFile:" << options->getBeamFile();
    char object[64] = {' '};
    if (options->getObject() == "sphere") strcpy (object,"sphere");
    if (options->getObject() == "inhomogeneous sphere") strcpy (object,"inhomosphere");
    if (options->getObject() == "random spheres (length)") strcpy (object,"randomsphere1");
    if (options->getObject() == "random spheres (meshsize)") strcpy (object,"randomsphere2");
    if (options->getObject() == "cube") strcpy (object,"cube");
    if (options->getObject() == "cuboid (length)") strcpy (object,"cuboid1");
    if (options->getObject() == "cuboid (meshsize)") strcpy (object,"cuboid2");
    if (options->getObject() == "inhomogeneous cuboid (length)") strcpy (object,"inhomocuboid1");
    if (options->getObject() == "inhomogeneous cuboid (meshsize)") strcpy (object,"inhomocuboid2");
    if (options->getObject() == "ellipsoid") strcpy (object,"ellipsoid");
    if (options->getObject() == "multiple spheres") strcpy (object,"nspheres");
    if (options->getObject() == "cylinder") strcpy (object,"cylinder");
    if (options->getObject() == "concentric spheres") strcpy (object,"concentricsphere");
    if (options->getObject() == "arbitrary") strcpy (object,"arbitrary");

    QLOG_DEBUG () << "Object:" << options->getObject();
    char namefileobj[64];
    for(int i = 0; i < 64; i++)
     namefileobj[i] = ' ';
    strncpy(namefileobj,(char*)options->getObjectFile().toStdString().c_str(),options->getObjectFile().size());
    QLOG_DEBUG () << "ObjectFile:" << options->getObjectFile();
    char anisotropy[3];
    for(int i = 0; i < 3; i++)
     anisotropy[i] = ' ';
    strncpy(anisotropy,(char*)options->getAnisotropy().toStdString().c_str(),options->getAnisotropy().size());
    QLOG_DEBUG () << "Anisotropy:" << options->getAnisotropy();
    char material[MAX_OBJECT_NUMBER][64];
    for( int i = 0; i < options->getObjectNumber(); i++ ) {
      for( int j = 0; j < 64; j++ ) {
       material[i][j] = ' ';
       strncpy(material[i],(char*)options->getMaterial().at(i).toStdString().c_str(),
	       options->getMaterial().at(i).size());
       QLOG_DEBUG () << "Material (" << i << ":" << options->getMaterial().at(i);
      }
    }
    char filereread[64];
    for(int i = 0; i < 64; i++)
     filereread[i] = ' ';
    strncpy(filereread,(char*)options->getFilereread().toStdString().c_str(),
			options->getFilereread().size());
    int discretization;
    discretization = options->getDiscretization();
    QLOG_DEBUG () << "Discretization:" << discretization;
    double tolerance;
    tolerance = options->getTolerance();
    QLOG_DEBUG () << "Tolerance:" << tolerance;
    char methodeit[12];
    for(int i = 0; i < 12; i++)
     methodeit[i] = ' ';
    strncpy(methodeit,(char*)options->getMethodeit().toStdString().c_str(),
			options->getMethodeit().size());
    QLOG_DEBUG () << "Methode:" << options->getMethodeit();
    char polarizability[2];
    for(int i = 0; i < 2; i++)
     polarizability[i] = ' ';
    strncpy(polarizability,(char*)options->getPolarizability().toStdString().c_str(),
			options->getPolarizability().size());
    QLOG_DEBUG () << "Polarizability:" << options->getPolarizability();
    int quad;
    quad = options->getQuad();
    // cdm.out file
    int localfieldCheck;
    localfieldCheck = options->getLocalfield();
    int macroscopicfieldCheck;
    macroscopicfieldCheck = options->getMacroscopicfield();
    int crosssectionCheck;
    crosssectionCheck = options->getCrosssection();
    int crosssectionpoyntingCheck;
    crosssectionpoyntingCheck = options->getCrosssectionpoynting();
    int quickdiffractCheck;
    quickdiffractCheck = options->getQuickdiffract();
    int nrigCheck;
    nrigCheck = options->getNrig();
    int opticalforceCheck;
    opticalforceCheck = options->getOpticalforce();
    int opticalforcedensityCheck;
    opticalforcedensityCheck = options->getOpticalforcedensity();
    int opticaltorqueCheck;
    opticaltorqueCheck = options->getOpticaltorque();
    int opticaltorquedensityCheck;
    opticaltorquedensityCheck = options->getOpticaltorquedensity();
    int nprocheCheck;
    options->setNproche(options->getNproche());
    if ( options->getObjectNumber() > 1 ) {
       if ( options->getNproche() == 0 )
          options->setNproche(1);
    }
    if ( options->getNearfield() == 0) {
      options->setNproche(0);
    }
    nprocheCheck = options->getNproche();    
    int microscopyCheck;
    microscopyCheck = options->getMicroscopy();
    int microscopyFFTCheck;
    microscopyFFTCheck = options->getMicroscopyFFT();
    int nenergieCheck;
    nenergieCheck = options->getNenergie();
    int dipolepsilon;
    dipolepsilon = options->getDipolepsilon();
    int nobjet;
    int farfieldCheck, nearfieldCheck, forceCheck;
    farfieldCheck = options->getFarfield();
    forceCheck = options->getForce();
    nearfieldCheck = options->getNearfield();
    if ( dipolepsilon == true )
      nobjet = true;
    else
      nobjet = false;
    int nxm,nym,nzm,nxmp,nymp,nzmp,nxx,nyy,nzz,ntheta,nphi;
    //     nxx = options->getNxx();
    //     nyy = options->getNyy();
   //      nzz = options->getNzz();
    nxm = options->getNxm();
    nym = options->getNym();
    nzm = options->getNzm();
    nxmp = options->getNxmp();
    nymp = options->getNymp();
    nzmp = options->getNzmp();
    ntheta = options->getNtheta();
    nphi = options->getNphi();
    if (options->getObject() == "cuboid (meshsize)" || options->getObject() == "random spheres (meshsize)" || options->getObject() == "inhomogeneous cuboid (meshsize)") {
      if (options->getNproche() !=2) {
	nxm = options->getNxx();
	nym = options->getNyy();
	nzm = options->getNzz();
	options->setNxm(options->getNxx());
	options->setNym(options->getNyy());
	options->setNzm(options->getNzz());
	nxmp = 0;
	nymp = 0;
	nzmp = 0;
      }
      else {
	nxm = options->getNxx()+2*options->getNxmp();
	nym = options->getNyy()+2*options->getNymp();
	nzm = options->getNzz()+2*options->getNzmp();
	options->setNxm(options->getNxx()+2*options->getNxmp());
	options->setNym(options->getNyy()+2*options->getNymp());
	options->setNzm(options->getNzz()+2*options->getNzmp());		
      }
    }

 
    if (options->getObject() == "arbitrary") {
      int nlength = options->getObjectFile().size();     
      string objfile(namefileobj);
      string objfilel = objfile.substr(0,nlength);
      cout<< objfilel  <<endl;
      ifstream monFlux(objfilel.c_str());
     
      if(monFlux){
    	monFlux  >> nxm >> nym >> nzm;
    	QLOG_DEBUG() << "Valeur max :" << nxm << nym << nzm;

	if (options->getNproche() !=2) {
	  options->setNxm(nxm);
	  options->setNym(nym);
	  options->setNzm(nzm);
	  nxmp = 0;
	  nymp = 0;
	  nzmp = 0;
	}
	else {
	  nxm =  nxm+2*options->getNxmp();
	  nym =  nym+2*options->getNymp();
	  nzm =  nzm+2*options->getNzmp();
	  options->setNxm(nxm);
	  options->setNym(nym);
	  options->setNzm(nzm);
	}
      }
      else{
    	cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
      }
    }
    
    double meshsize;
    meshsize = options->getMeshsize();
    QLOG_DEBUG() << "NMAX:" << nxm*nym*nzm << "NTHETA:" << ntheta << "NPHI:" << nphi;
    int nfft2d;
    nfft2d = options->getnfft2d();
    double numaper;
    numaper = options->getNA();
    double numaperinc;
    numaperinc = options->getNAinc();    
    double gross;   
    gross = options->getGross();
    double zlens;   
    zlens = options->getZlens();
    int nside;
    nside = options->getNside();
    int ntypemic;
    ntypemic = options->getNtypemic();

    QLOG_DEBUG () << "Beam:" << options->getNtypemic() << ntypemic;
    QLOG_DEBUG () << "Beam:" << options->getNside() << nside;

    int nreadCheck;
    nreadCheck = options->getNread();
    int nmatlabCheck;
    nmatlabCheck = options->getNmatlab();
    char fileh5[100];
    for(int i = 0; i < 100; i++)
      fileh5[i] = ' ';
    strncpy(fileh5,(char*)options->getH5File().toStdString().c_str(),
	    options->getH5File().size());
    
    int advancedinterfaceCheck;
    advancedinterfaceCheck = options->getAdvancedinterface();    
   // sphere.in / cube.in files ( includes multiple spheres)
   double density;
   density  = options->getDensity();
   QLOG_DEBUG () << "Density:" << QString::number(density,'g',5);
   double side;
   side  = options->getCubeside();
   QLOG_DEBUG () << "Cube side:" << QString::number(side,'g',5);
   double sidex;
   sidex  = options->getCubesidex();
   QLOG_DEBUG () << "Cuboid side x:" << QString::number(sidex,'g',5);
   double sidey;
   sidey  = options->getCubesidey();
   QLOG_DEBUG () << "Cuboid side y:" << QString::number(sidey,'g',5);
   double sidez;
   sidez  = options->getCubesidez();
   QLOG_DEBUG () << "Cuboid side z:" << QString::number(sidez,'g',5);
   double hauteur;
   hauteur = options->getHauteur();
   QLOG_DEBUG () << "Hauteur :" << QString::number(hauteur,'g',5);
   int objectnumber;
   objectnumber = options->getObjectNumber();
   int seed;
   seed = options->getSphereseed();
   QLOG_DEBUG () << "Seed :" << QString::number(seed,'g',5);
   double lc;
   lc = options->getSpherecoherencelength();
   QLOG_DEBUG () << "Coherence length :" << QString::number(lc,'g',5);
   double hc;
   hc = options->getSpherestandardev();
   QLOG_DEBUG () << "Standard deviation :" << QString::number(hc,'g',5);
   double rayon[MAX_OBJECT_NUMBER];
   double xg[MAX_OBJECT_NUMBER];
   double yg[MAX_OBJECT_NUMBER];
   double zg[MAX_OBJECT_NUMBER];
   dcmplx eps[MAX_OBJECT_NUMBER];
   dcmplx epsani[MAX_OBJECT_NUMBER][3][3];
  
   for (int i = 0 ; i < objectnumber; i++) {
     rayon[i] = options->getSphereradius().at(i);
     QLOG_DEBUG () << "Sphere radius[" << i << "]=" << QString::number(rayon[i],'g',5);
     xg[i] = options->getPositionx().at(i);
     QLOG_DEBUG () << "Position X[" << i << "]=" << QString::number(xg[i],'g',5);
     yg[i] = options->getPositiony().at(i);
     QLOG_DEBUG () << "Position Y[" << i << "]=" << QString::number(yg[i],'g',5);
     zg[i] = options->getPositionz().at(i);
     QLOG_DEBUG () << "Position Z[" << i << "]=" << QString::number(zg[i],'g',5);
     eps[i] = options->getEpsilon().at(i);
     QLOG_DEBUG () << "Epsilon [:" << i << "] = " 
		 << QString::number(real(eps[i]),'g',5) 
                 << " " <<  QString::number(imag(eps[i]),'g',5);
     epsani[i][0][0] = options->getEpsilon11().at(i);
     QLOG_DEBUG () << "Epsilon[0][0][" << i << "] = " 
		 << QString::number(real(epsani[i][0][0]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][0][0]),'g',5);
     epsani[i][1][0] = options->getEpsilon12().at(i);
     QLOG_DEBUG () << "Epsilon[1][0][" << i << "] = " 
		 << QString::number(real(epsani[i][1][0]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][1][0]),'g',5);
     epsani[i][2][0] = options->getEpsilon13().at(i);
     QLOG_DEBUG () << "Epsilon[2][0][" << i << "] = " 
		 << QString::number(real(epsani[i][2][0]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][2][0]),'g',5);
     epsani[i][0][1] = options->getEpsilon21().at(i);
     QLOG_DEBUG () << "Epsilon[0][1][" << i << "] = " 
		 << QString::number(real(epsani[i][0][1]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][0][1]),'g',5);
     epsani[i][1][1] = options->getEpsilon22().at(i);
     QLOG_DEBUG () << "Epsilon[1][1][" << i << "] = " 
		 << QString::number(real(epsani[i][1][1]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][1][1]),'g',5);
     epsani[i][2][1] = options->getEpsilon23().at(i);
     QLOG_DEBUG () << "Epsilon[2][1][" << i << "] = " 
		 << QString::number(real(epsani[i][2][1]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][2][1]),'g',5);
     epsani[i][0][2] = options->getEpsilon31().at(i);
     QLOG_DEBUG () << "Epsilon[0][2][" << i << "] = " 
		 << QString::number(real(epsani[i][0][2]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][0][2]),'g',5);
     epsani[i][1][2] = options->getEpsilon32().at(i);
     QLOG_DEBUG () << "Epsilon[1][2][" << i << "] = " 
		 << QString::number(real(epsani[i][1][2]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][1][2]),'g',5);
     epsani[i][2][2] = options->getEpsilon33().at(i);
     QLOG_DEBUG () << "Epsilon[2][2][" << i << "] = " 
		 << QString::number(real(epsani[i][2][2]),'g',5) 
                 << " " <<  QString::number(imag(epsani[i][2][2]),'g',5);
   }
   // ellipsoid
   double demiaxea;
   demiaxea = options->getDemiaxea();
   QLOG_DEBUG () << "Half axe A:" << QString::number(demiaxea,'g',5);
   double demiaxeb;
   demiaxeb = options->getDemiaxeb();
   QLOG_DEBUG () << "Half axe B:" << QString::number(demiaxeb,'g',5);
   double demiaxec; 
   demiaxec = options->getDemiaxec();
   QLOG_DEBUG () << "Half axe C:" << QString::number(demiaxec,'g',5);
   int nfftadin;
   int i;
   if (advancedinterfaceCheck == 0 && dipolepsilon == false) {
     if (quickdiffractCheck == 1 || microscopyFFTCheck == 1 || nenergieCheck == 1) {
       QLOG_DEBUG () << "fft:" << nfft2d << wavelength << rayon[0] << nxm;
       
     
       if (options->getObject() == "cuboid (meshsize)" || options->getObject() == "random spheres (meshsize)" || options->getObject() == "inhomogeneous cuboid (meshsize)") {
	 nfftadin = 10*wavelength/meshsize;
       }
       if (options->getObject() == "sphere" || options->getObject() == "inhomogeneous sphere") {
	 nfftadin = 10*wavelength/rayon[0]*nxm;
       }
       if (options->getObject() == "cube") {
	 nfftadin = 10*wavelength/side*nxm;
       }
       if (options->getObject() == "random spheres (length)" || options->getObject() == "cuboid (length)" || options->getObject() == "inhomogeneous cuboid (length)") {
	 side = max(sidex,sidey);
	 side = max(side,sidez);
	 nfftadin = 10*wavelength/side*nxm;
       }
       if (options->getObject() == "ellipsoid") {
	 side = max(demiaxea,demiaxeb);
	 side = max(side,demiaxec);
	 nfftadin = 5*wavelength/side*nxm;
       }
       if (options->getObject() == "cylinder") {
	 side = max(hauteur,rayon[0]);
	 nfftadin = 10*wavelength/side*nxm;
       }
       if (options->getObject() == "concentric spheres") {
	 nfftadin = 10*wavelength/rayon[objectnumber-1]*nxm;
       }
       if (options->getObject() == "multiple spheres") {
	 double xmax,xmin,ymax,ymin,zmax,zmin;
	 xmax=xg[0];
	 xmin=xg[0];
	 ymax=yg[0];
	 ymin=yg[0];
	 zmax=zg[0];
	 zmin=zg[0];
	 for (int i = 0 ; i < objectnumber; i++) {
	   xmax=max(xmax,xg[i]+rayon[i]);
	   xmin=min(xmin,xg[i]-rayon[i]);
	   ymax=max(ymax,yg[i]+rayon[i]);
	   ymin=max(ymin,yg[i]-rayon[i]);
	   zmax=max(zmax,zg[i]+rayon[i]);
	   zmin=min(zmin,zg[i]-rayon[i]);
	 }
	 side=max(xmax-xmin,ymax-ymin);
	 side=max(side,zmax-zmin);
	 nfftadin = 10*wavelength/(side)*nxm;	  
       }
       i=1;
       while( nfftadin > pow(2,i) ) {
	 i++;
	 QLOG_DEBUG () << "fft:" << i,nxm;
       }
       nfft2d=pow(2,i);
       options->setnfft2d(nfft2d);
       nmatlabCheck=1;
       options->setNmatlab(nmatlabCheck);
     }
   }
   QLOG_DEBUG () << "fftv:" << nfft2d;
   if (nfft2d > 4096) {
     *infoMessage = QString("Meshsize too small for microscopy");
     return;
   }
   if (options->getBeam() == "Arbitrary wave (file)"){
     *infoMessage = QString("Arbitrary wave only with advanced interface");
     return;
   }
   if (options->getBeam() == "arbitrary"){
     *infoMessage = QString("Arbitrary object only with advanced interface");
     return;
   }
   double thetaobj;
   thetaobj = options->getThetaobj();
   QLOG_DEBUG () << "Theta object:" << QString::number(thetaobj,'g',5);
   double phiobj;
   phiobj = options->getPhiobj();
   QLOG_DEBUG () << "Phi object:" << QString::number(phiobj,'g',5);
   double psiobj;
   psiobj = options->getPsiobj();
   QLOG_DEBUG () << "Psi object:" << QString::number(psiobj,'g',5);

   // planewavecircular.in files
   double theta;
   theta = options->getIncidenceangle_theta_z();
   QLOG_DEBUG () << "Theta:" << QString::number(theta,'g',5);
   double phi;
   phi = options->getIncidenceangle_phi_x();
   QLOG_DEBUG () << "Phi:" << QString::number(phi,'g',5);
   double pp;
   pp = options->getPolarizationTM();
   QLOG_DEBUG () << "Polarization TM:" << QString::number(pp,'g',5);
   double ss;
   if (options->getBeam() == "Linear plane wave")
     ss = options->getPolarizationTE();
   else if (options->getBeam() == "Circular plane wave")
     ss = options->getPolarizationRL();
   QLOG_DEBUG () << "Polarization TE/RL:" << QString::number(ss,'g',5);
   double P0;
   P0 = options->getP0();
   QLOG_DEBUG () << "P0:" << QString::number(P0,'g',5);
   double W0;
   W0 = options->getW0();
   QLOG_DEBUG () << "W0:" << QString::number(W0,'g',5);
   double xgaus;
   xgaus = options->getXgaus();
   QLOG_DEBUG () << "Xgaus:" << QString::number(xgaus,'g',5);
   double ygaus;
   ygaus = options->getYgaus();
   QLOG_DEBUG () << "Ygaus:" << QString::number(ygaus,'g',5);
   double zgaus;
   zgaus = options->getZgaus();
   QLOG_DEBUG () << "Zgaus:" << QString::number(zgaus,'g',5);

   // wave multi
   double thetam[MAX_WAVEMULTI_NUMBER];
   double phim[MAX_WAVEMULTI_NUMBER];
   double ppm[MAX_WAVEMULTI_NUMBER];
   double ssm[MAX_WAVEMULTI_NUMBER];
   dcmplx E0m[MAX_WAVEMULTI_NUMBER];
   int wavemultinumber;
   wavemultinumber = options->getWaveMultiNumber();
   for (int i = 0 ; i < wavemultinumber; i++) {
     thetam[i] = options->getThetam().at(i);
     phim[i] = options->getPhim().at(i);
     ppm[i] = options->getPpm().at(i);
     ssm[i] = options->getSsm().at(i);
     E0m[i] = options->getE0m().at(i);
   }
   // Initialize info string
   char infostr[64];
   for(int i = 0; i < 64; i++) 
     infostr[i] = ' ';
   // return results
   int objectsubunits;
   int meshsubunits;
   double lambda10n;
   double k0;
   int nmaxpp;
   double toleranceobtained;
   int numberofax1, numberofax2;
   double reflectivity;
   double absorptivity;
   double transmittivity;
   double extinctioncrosssection;
   double absorbingcrosssection;
   double crosssection;
   double crosssectionpoynting;
   double assymetricparam;
   double irra;
   dcmplx E0;
   double opticalforce[3], opticalforcedensity;
   double opticaltorque[3], opticaltorquedensity;
   double *incidentfield,*localfield,*macroscopicfield;
   double *forcex,*forcey,*forcez;
   double *forcexmulti,*forceymulti,*forcezmulti;
   double *torquex,*torquey,*torquez;
   double *torquexmulti,*torqueymulti,*torquezmulti;
   double *xc,*yc,*zc;
   double *kxy, *xy;
   double *xcwf,*ycwf,*zcwf;
   double *thetafield,*phifield,*poyntingfield;
   dcmplx *incidentfieldx, *incidentfieldy, *incidentfieldz;
   dcmplx *localfieldx, *localfieldy, *localfieldz;
   dcmplx *macroscopicfieldx, *macroscopicfieldy, *macroscopicfieldz;
   dcmplx *polarisafield, *epsilonfield;
   dcmplx *eimagex, *eimagey, *eimagez;
   dcmplx *efourierx, *efouriery, *efourierz;
   dcmplx *efourierincx, *efourierincy, *efourierincz;
   dcmplx *eimageincx, *eimageincy, *eimageincz;
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
   dcmplx *FF, *FF0, *FFloc, *xr, *xi;
   dcmplx *wrk;
   dcmplx *FFTTENSORxx, *FFTTENSORxy, *FFTTENSORxz;
   dcmplx *FFTTENSORyy, *FFTTENSORyz, *FFTTENSORzz;
   dcmplx *vectx, *vecty, *vectz;
   dcmplx *Ediffkzpos, *Ediffkzneg;
   int *Tabdip, *Tabmulti;

   // Clean up memory
   run->cleanVectorsMemory();
   // Allocate memory
   int available_mem = run->checkAvailableMemorySize() * 0.9;
   int needed_mem = run->allocateVectorsMemory(nxm*nym*nzm, ntheta, nphi, nfft2d, objectnumber);
   if ( needed_mem == -1 ) {
      *infoMessage = QString("Memory allocation failed !");
      QLOG_FATAL() << "Not enough memory ! Available memory=" << available_mem << "MB";
      return;
    }
   else if ( needed_mem >= available_mem ) {
     *infoMessage = QString("Memory allocation failed !");
     QLOG_FATAL() << "Not enough memory ! Needed  (" << needed_mem
		  << "MB) exceeds available (" << available_mem << "MB)";
     return;
    }
   else
     QLOG_INFO() << "Memory used=" << needed_mem << "MB (available memory="
                 << available_mem << "MB)";
    
   incidentfield = run->getIncidentField();
   localfield = run->getLocalField();
   macroscopicfield = run->getMacroscopicField();
   xc = run->getXc();
   yc = run->getYc();
   zc = run->getZc();
   xcwf = run->getXcWF();
   ycwf = run->getYcWF();
   zcwf = run->getZcWF();
   forcex = run->getForceX();
   forcey = run->getForceY();
   forcez = run->getForceZ();
   forcexmulti = run->getForceXMulti();
   forceymulti = run->getForceYMulti();
   forcezmulti = run->getForceZMulti();
   torquex = run->getTorqueX();
   torquey = run->getTorqueY();
   torquez = run->getTorqueZ();
   torquexmulti = run->getTorqueXMulti();
   torqueymulti = run->getTorqueYMulti();
   torquezmulti = run->getTorqueZMulti();
   phifield = run->getPhiField();
   thetafield = run->getThetaField();
   poyntingfield = run->getPoyntingField();
   incidentfieldx = run->getIncidentFieldX();
   localfieldx = run->getLocalFieldX();
   macroscopicfieldx = run->getMacroscopicFieldX();
   incidentfieldy = run->getIncidentFieldY();
   localfieldy = run->getLocalFieldY();
   macroscopicfieldy = run->getMacroscopicFieldY();
   incidentfieldz = run->getIncidentFieldZ();
   localfieldz = run->getLocalFieldZ();
   macroscopicfieldz = run->getMacroscopicFieldZ();
   polarisafield = run->getPolarisaField();
   epsilonfield = run->getEpsilonField();
   xy = run->getXY();
   kxy = run->getKXY();
   eimagex = run->getEimageX();
   eimagey = run->getEimageY();
   eimagez = run->getEimageZ();
   efourierx = run->getEfourierX();
   efouriery = run->getEfourierY();
   efourierz = run->getEfourierZ();
   efourierincx = run->getEfourierincX();
   efourierincy = run->getEfourierincY();
   efourierincz = run->getEfourierincZ();
   eimageincx = run->getEimageincX();
   eimageincy = run->getEimageincY();
   eimageincz = run->getEimageincZ();
   
   FF = run->getFF();
   FF0 = run->getFF0();
   FFloc = run->getFFloc();
   xr = run->getxr();
   xi = run->getxi();
   wrk = run->getwrk();
   FFTTENSORxx = run->getFFTTENSORxx();
   FFTTENSORxy = run->getFFTTENSORxy();
   FFTTENSORxz = run->getFFTTENSORxz();
   FFTTENSORyy = run->getFFTTENSORyy();
   FFTTENSORyz = run->getFFTTENSORyz();
   FFTTENSORzz = run->getFFTTENSORzz();
   vectx = run->getvectx();
   vecty = run->getvecty();
   vectz = run->getvectz();
   Ediffkzpos = run->getEdiffkzpos();
   Ediffkzneg = run->getEdiffkzneg();
   Tabdip = run->getTabdip();
   Tabmulti = run->getTabmulti();
   
   cdmlib_(&wavelength, beam, object, anisotropy, material,
	   &discretization, &tolerance, methodeit, polarizability, &quad, &nreadCheck, filereread, &nmatlabCheck, fileh5,
	   &localfieldCheck, &macroscopicfieldCheck, &crosssectionCheck,
	   &crosssectionpoyntingCheck, &quickdiffractCheck, &nrigCheck, 
           &opticalforceCheck, &opticalforcedensityCheck,
	   &opticaltorqueCheck, &opticaltorquedensityCheck, &nprocheCheck, 
	   &microscopyCheck, &microscopyFFTCheck, &nenergieCheck, &nobjet,
           // cube,sphere (includes multiple)
	   & density, &side, &sidex, &sidey, &sidez, &hauteur,
           &objectnumber, rayon, xg, yg, zg, eps, (dcmplx*)epsani, &lc, &hc, &seed,
	   // ellipsoid
           &demiaxea, &demiaxeb, &demiaxec, 
	   &thetaobj, &phiobj, &psiobj, namefileobj,
           // planewave
	   &theta, &phi, &pp, &ss, &P0, &W0, &xgaus, &ygaus, &zgaus, namefileinc,
           // wave multi
           thetam, phim, ppm, ssm, E0m, &wavemultinumber,
	   infostr, stopFlag,
	   &objectsubunits, &meshsubunits, &meshsize, 
           &lambda10n, &k0, &toleranceobtained, &numberofax1, &numberofax2,
           &absorptivity, &reflectivity, &transmittivity,
           &extinctioncrosssection, &absorbingcrosssection, 
           &crosssection, &crosssectionpoynting, &assymetricparam,
           &irra, &E0,
           opticalforce, &opticalforcedensity,
           opticaltorque, &opticaltorquedensity,
           &nxm, &nym, &nzm, &nxmp, &nymp, &nzmp, &nmaxpp,
           incidentfield, localfield, macroscopicfield,
	   xc, yc, zc, xcwf, ycwf, zcwf, 
           &ntheta, &nphi, thetafield, phifield, poyntingfield, 
	   forcex, forcey, forcez, forcexmulti, forceymulti, forcezmulti,
           torquex, torquey, torquez, torquexmulti, torqueymulti, torquezmulti,
           (dcmplx*)incidentfieldx, (dcmplx*)incidentfieldy, (dcmplx*)incidentfieldz,
	   (dcmplx*)localfieldx, (dcmplx*)localfieldy, (dcmplx*)localfieldz,
           (dcmplx*)macroscopicfieldx, (dcmplx*)macroscopicfieldy, (dcmplx*)macroscopicfieldz,
	   (dcmplx*)polarisafield, (dcmplx*)epsilonfield, &nfft2d, (dcmplx*)eimagex, (dcmplx*)eimagey,
           (dcmplx*)eimagez, (dcmplx*)eimageincx, (dcmplx*)eimageincy, (dcmplx*)eimageincz, 
	   (dcmplx*)efourierx, (dcmplx*)efouriery, (dcmplx*)efourierz,
	   (dcmplx*)efourierincx, (dcmplx*)efourierincy, (dcmplx*)efourierincz,
	   kxy, xy, &numaper, &numaperinc, &gross, &zlens, &ntypemic, &nside,
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
//     taille double complex (3*nxm*nym*nzm)
           (dcmplx*)FF, (dcmplx*)FF0, (dcmplx*)FFloc, (dcmplx*)xr, (dcmplx*)xi,
//     taille double complex (3*nxm*nym*nzm,12)
           (dcmplx*)wrk,
//     taille double complex (8*nxm*nym*nzm)
           (dcmplx*)FFTTENSORxx, (dcmplx*)FFTTENSORxy, (dcmplx*)FFTTENSORxz, 
           (dcmplx*)FFTTENSORyy, (dcmplx*)FFTTENSORyz,
           (dcmplx*)FFTTENSORzz, (dcmplx*)vectx, (dcmplx*)vecty, (dcmplx*)vectz,
//     taille double complex (nfft2d,nfft2d,3)
           (dcmplx*)Ediffkzpos,(dcmplx*)Ediffkzneg,
//     taille entier (nxm*nym*nzm)
           (int*)Tabdip, (int*)Tabmulti
	   );
   /*dcmplx *_epstest = run->getEpsilonField();
       for (int k = 0 ; k < nxm*nym*nzm*3*3; k++)
         QLOG_DEBUG() << " RESULT fabs(Epsilonfield [" << k << "])=" 
		  << fabs(_epstest[k]) ;*/
   *infoMessage = QString(infostr);
   run->setObjectSubunits(objectsubunits);
   run->setMeshSubunits(meshsubunits);
   run->setMeshSize(meshsize);
   run->setLambda10n(lambda10n);
   run->setK0(k0);
   run->setNmaxpp(nmaxpp);
   run->setToleranceObtained(toleranceobtained);
   run->setNumberofAx1(numberofax1);
   run->setNumberofAx2(numberofax2);
   run->setReflectivity(reflectivity);
   run->setAbsorptivity(absorptivity);
   run->setTransmittivity(transmittivity);
   run->setExtinctionCrossection(extinctioncrosssection);
   run->setAbsorbingCrossection(absorbingcrosssection);
   run->setScatteringCrossection(crosssection);
   run->setScatteringCrossectionWithIntegration(crosssectionpoynting);
   run->setScatteringAssymetricParam(assymetricparam);
   run->setIrra(irra);
   run->setE0(E0);
   run->setOpticalForcex(opticalforce[0]);
   run->setOpticalForcey(opticalforce[1]);
   run->setOpticalForcez(opticalforce[2]);
   run->setOpticalForceModulus(opticalforcedensity);
   run->setOpticalTorquex(opticaltorque[0]);
   run->setOpticalTorquey(opticaltorque[1]);
   run->setOpticalTorquez(opticaltorque[2]);
   run->setOpticalTorqueModulus(opticaltorquedensity);
//   QLOG_INFO() << "cdmlibwrapper> Thread finished";
}
void 
RunWidget::displayFinishedBox()
{
   if (canceldlg != NULL) {delete canceldlg; canceldlg = NULL;}
   if (progressdlg != NULL) {delete progressdlg; progressdlg = NULL;}
   QMessageBox::about(this,"Result",(*infoMessage).split("!")[0]);
   if (*stopFlag == 0)
     displayResults();
   
}
void 
RunWidget::displayResults()
{
   QLOG_INFO () << "!!! THE RESULTS !!!"; 
   int nmax = options->getNxm()*options->getNym()*options->getNzm();
   int nmaxs = (options->getNxm() - 1)*(options->getNym() - 1)*(options->getNzm() - 1);
   // Scalar results
   QFrame *hsep0 = new QFrame(this);
   QFrame *hsep00 = new QFrame(this);
   hsep0->setFrameShape(QFrame::HLine);
   hsep0->setFrameShadow(QFrame::Sunken);
   hsep00->setFrameShape(QFrame::HLine);
   hsep00->setFrameShadow(QFrame::Sunken);
   QLabel *objsub = new QLabel(QString::number(run->getObjectSubunits()),this);
   outputlayout->addRow("Object subunits:", objsub);
   QLabel *meshsub = new QLabel(QString::number(run->getMeshSubunits()),this);
   outputlayout->addRow("Mesh subunits:", meshsub);
   QLabel *meshsize = new QLabel(QString::number(run->getMeshSize()),this);
   QLabel *meshsizeLabel = new QLabel("<html><body>Mesh size [m]:</body></html>",this);
   outputlayout->addRow(meshsizeLabel, meshsize);
   QLabel *lambda10n = new QLabel(QString::number(run->getLambda10n()),this);
   QLabel *lambdaLabel = new QLabel("<html><body><sup>&lambda;</sup>&frasl;<sub>(10|n|)</sub> [m]:</body></html>",this);
   outputlayout->addRow(lambdaLabel, lambda10n);
   if (options->getDipolepsilon() == false ) {
     QLabel *k0 = new QLabel(QString::number(run->getK0()),this);
     QLabel *kappa0Label = 
	new QLabel("<html><body>k<sub>0</sub> [m<sup>-1</sup>]:</body></html>",this);
     outputlayout->addRow(kappa0Label, k0);
     QLabel *irra = new QLabel(QString::number(run->getIrra()),this);
     outputlayout->addRow("Irradiance [W/m<sup>2</sup>]:", irra);
     QLabel *E0 = new QLabel(QString::number(abs(run->getE0())),this);
     outputlayout->addRow("Field modulus [V/m]:", E0);
     QLabel *toleranceobtained = 
	new QLabel(QString::number(run->getToleranceObtained()),this);
     outputlayout->addRow("Tolerance obtained:", toleranceobtained);
     QLabel *numberofax = new QLabel(QString::number(run->getNumberofAx1()) + 
			"(" + QString::number(run->getNumberofAx2()) + ")",this);
     outputlayout->addRow("Number of products Ax (iterations):", numberofax);
   }
   if (options->getNenergie()) {
      QLabel *absorptivity = new 
      QLabel(QString::number(1 - run->getAbsorptivity()),this);
      QLabel *absorptivityLabel = 
	new QLabel("<html><body>Absorptivity:</body></html>",this);
      outputlayout->addRow(absorptivityLabel, absorptivity);
      QLabel *reflectivity = new 
      QLabel(QString::number(run->getReflectivity()),this);
      QLabel *reflectivityLabel = 
	new QLabel("<html><body>Reflectivity:</body></html>",this);    
      outputlayout->addRow(reflectivityLabel, reflectivity);
      QLabel *transmittivity = new 
      QLabel(QString::number(run->getTransmittivity()),this);
      QLabel *transmittivityLabel = 
	new QLabel("<html><body>Transmittivity:</body></html>",this);
      outputlayout->addRow(transmittivityLabel, transmittivity);
   }
   if (options->getCrosssection()) {
      QLabel *extinctioncrosssection = new 
      QLabel(QString::number(run->getExtinctionCrossection()),this);
      QLabel *extinctionLabel = 
	new QLabel("<html><body>Extinction cross section [m<sup>2</sup>]:</body></html>",this);
      outputlayout->addRow(extinctionLabel, extinctioncrosssection);
      QLabel *absorbingcrosssection = new 
      QLabel(QString::number(run->getAbsorbingCrossection()),this);
      QLabel *absorbingLabel = 
	new QLabel("<html><body>Absorbing cross section [m<sup>2</sup>]:</body></html>",this);    
      outputlayout->addRow(absorbingLabel, absorbingcrosssection);
      QLabel *crosssection = new 
      QLabel(QString::number(run->getScatteringCrossection()),this);
      QLabel *scatteringLabel = 
	new QLabel("<html><body>Scattering cross section [m<sup>2</sup>]:</body></html>",this);
      outputlayout->addRow(scatteringLabel, crosssection);
   }
   if (options->getCrosssectionpoynting()) {
   QLabel *crosssectionpoynting = new 
      QLabel(QString::number(run->getScatteringCrossectionWithIntegration()),this);
      QLabel *scatteringpoyntingLabel = 
	new QLabel("<html><body>Scattering cross section  with integration [m<sup>2</sup>]:</body></html>",this);
      outputlayout->addRow(scatteringpoyntingLabel, crosssectionpoynting);
      QLabel *assymetricparam = new 
	QLabel(QString::number(run->getScatteringAssymetricParam()),this);
      outputlayout->addRow("Scattering assymetric parameter:", assymetricparam);
   }
   if (options->getOpticalforce()) {
      QLabel *opticalforcex = new 
	QLabel(QString::number(run->getOpticalForcex()),this);
      outputlayout->addRow("Optical force x [N]:", opticalforcex);
      QLabel *opticalforcey = new 
	QLabel(QString::number(run->getOpticalForcey()),this);
      outputlayout->addRow("Optical force y [N]:", opticalforcey);
      QLabel *opticalforcez = new 
	QLabel(QString::number(run->getOpticalForcez()),this);
      outputlayout->addRow("Optical force z [N]:", opticalforcez);
      if (options->getObjectNumber() > 1 ) {
       for ( int i = 0 ; i < options->getObjectNumber(); i++ ) {
         QLabel *opticalforcex = new 
	   QLabel(QString::number(run->getForceXMulti()[i]),this);
         outputlayout->addRow("Optical force x (object "+QString::number(i+1)+") [N]:",
	 opticalforcex);
         QLabel *opticalforcey = new 
	   QLabel(QString::number(run->getForceYMulti()[i]),this);
         outputlayout->addRow("Optical force y (object "+QString::number(i+1)+") [N]:",
	 opticalforcey);
          QLabel *opticalforcez = new 
	   QLabel(QString::number(run->getForceZMulti()[i]),this);
         outputlayout->addRow("Optical force z (object "+QString::number(i+1)+") [N]:",
	 opticalforcez);
       }
      }
      QLabel *opticalforcedensity = new 
	QLabel(QString::number(run->getOpticalForceModulus()),this);
      outputlayout->addRow("Optical force modulus [N]:", opticalforcedensity);
   }
   if (options->getOpticaltorque()) {
       QLabel *opticaltorquex = new 
	QLabel(QString::number(run->getOpticalTorquex()),this);
      outputlayout->addRow("Optical torque x [N/m]:", opticaltorquex);
      QLabel *opticaltorquey = new 
	QLabel(QString::number(run->getOpticalTorquey()),this);
      outputlayout->addRow("Optical torque y [N/m]:", opticaltorquey);
      QLabel *opticaltorquez = new 
	QLabel(QString::number(run->getOpticalTorquez()),this);
      outputlayout->addRow("Optical torque z [N/m]:", opticaltorquez);
      if (options->getObjectNumber() > 1 ) {
        for ( int i = 0 ; i < options->getObjectNumber(); i++ ) {
         QLabel *opticaltorquex = new 
	   QLabel(QString::number(run->getTorqueXMulti()[i]),this);
         outputlayout->addRow("Optical torque x (object "+QString::number(i+1)+") [N/m]:",
	 opticaltorquex);
         QLabel *opticaltorquey = new 
	   QLabel(QString::number(run->getTorqueYMulti()[i]),this);
         outputlayout->addRow("Optical torque y (object "+QString::number(i+1)+") [N/m]:",
	 opticaltorquey);
          QLabel *opticaltorquez = new 
	   QLabel(QString::number(run->getTorqueZMulti()[i]),this);
         outputlayout->addRow("Optical torque z (object "+QString::number(i+1)+") [N/m]:",
	 opticaltorquez);
        }
      }
      QLabel *opticaltorquedensity = new 
	QLabel(QString::number(run->getOpticalTorqueModulus()),this);
      outputlayout->addRow("Optical torque modulus [N/m]:", opticaltorquedensity);
   }
   ////////////////////////////////////////////////////////////////////////////////////
   // Dipole epsilon buttons
   QPushButton *dipoles3DButton = new QPushButton("Plot epsilon/dipoles",this);
   dipoles3DButton->setFixedWidth(200);
   connect(dipoles3DButton, SIGNAL(clicked()), this, SLOT(dipoles3Dplot()));
   ////////////////////////////////////////////////////////////////////////////////////
   // Near Field buttons
   // Plot All Buttons
   QPushButton *plotallxnearfieldButton = new QPushButton("Plot all X",this);
   connect(plotallxnearfieldButton, SIGNAL(clicked()), this, SLOT(plotallxnearfield()));
   QPushButton *plotallynearfieldButton = new QPushButton("Plot all Y",this);
   connect(plotallynearfieldButton, SIGNAL(clicked()), this, SLOT(plotallynearfield()));
   QPushButton *plotallznearfieldButton = new QPushButton("Plot all Z",this);
   connect(plotallznearfieldButton, SIGNAL(clicked()), this, SLOT(plotallznearfield()));
   // Plot X,Y,Z Buttons
   QPushButton *plotxnearfieldButton = new QPushButton("Plot X",this);
   connect(plotxnearfieldButton, SIGNAL(clicked()), this, SLOT(plotxnearfield()));
   QPushButton *plotynearfieldButton = new QPushButton("Plot Y",this);
   connect(plotynearfieldButton, SIGNAL(clicked()), this, SLOT(plotynearfield()));
   QPushButton *plotznearfieldButton = new QPushButton("Plot Z",this);
   connect(plotznearfieldButton, SIGNAL(clicked()), this, SLOT(plotznearfield()));
   ////////////////////////////////////////////////////////////////////////////////////
   // Far Field buttons (microscopy, poynting)
   QPushButton *poyntingButton = new QPushButton("Plot Poynting",this);
   poyntingButton->setFixedWidth(200);
   connect(poyntingButton, SIGNAL(clicked()), this, SLOT(poyntingplot()));
   QPushButton *plotmicroscopyButton = new QPushButton("Plot microscopy",this);
   connect(plotmicroscopyButton, SIGNAL(clicked()), this, SLOT(plotmicroscopy()));
   ////////////////////////////////////////////////////////////////////////////////////
   // Force buttons
   // Plot All Buttons
   QPushButton *plotallxforceButton = new QPushButton("Plot all X",this);
   connect(plotallxforceButton, SIGNAL(clicked()), this, SLOT(plotallxforce()));
   QPushButton *plotallyforceButton = new QPushButton("Plot all Y",this);
   connect(plotallyforceButton, SIGNAL(clicked()), this, SLOT(plotallyforce()));
   QPushButton *plotallzforceButton = new QPushButton("Plot all Z",this);
   connect(plotallzforceButton, SIGNAL(clicked()), this, SLOT(plotallzforce()));
   // Plot X,Y,Z Buttons
   QPushButton *plotxforceButton = new QPushButton("Plot X",this);
   connect(plotxforceButton, SIGNAL(clicked()), this, SLOT(plotxforce()));
   QPushButton *plotyforceButton = new QPushButton("Plot Y",this);
   connect(plotyforceButton, SIGNAL(clicked()), this, SLOT(plotyforce()));
   QPushButton *plotzforceButton = new QPushButton("Plot Z",this);
   connect(plotzforceButton, SIGNAL(clicked()), this, SLOT(plotzforce()));
   int nmaxpp=run->getNmaxpp();
   // CHECK DELTA OVER Xc,Yc,Zc, XcWF, YcWF, ZcWF
   for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
     if ( fabs(run->getXc()[i]) < DELTA ) run->getXc()[i] = 0;
     if ( fabs(run->getYc()[i]) < DELTA ) run->getYc()[i] = 0;
     if ( fabs(run->getZc()[i]) < DELTA ) run->getZc()[i] = 0;
   }
   for (int i = 0 ; i < nmaxpp ; i++) {
     if ( fabs(run->getXcWF()[i]) < DELTA ) run->getXcWF()[i] = 0;
     if ( fabs(run->getYcWF()[i]) < DELTA ) run->getYcWF()[i] = 0;
     if ( fabs(run->getZcWF()[i]) < DELTA ) run->getZcWF()[i] = 0;
   }
   // Plot X list
   double refx;
   QList<double> xlist;
   QList<double> xwflist;
   refx = run->getXc()[0];
   xlist.append(refx);
   for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
       bool exists = false;
       if ( run->getXc()[i] != refx) {
         for (int j = 0 ; j < xlist.size() ; j++) {
           if ( abs(run->getXc()[i] - xlist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
          }
          if ( exists == false ) {
            refx = run->getXc()[i];
            xlist.append(refx);
          }
       } 
   }
   refx = run->getXcWF()[0];
   xwflist.append(refx);
   for (int i = 0 ; i < nmaxpp ; i++) {
       bool exists = false;
       if ( run->getXcWF()[i] != refx) {
         for (int j = 0 ; j < xwflist.size() ; j++) {
           if ( abs(run->getXcWF()[i] - xwflist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
          }
          if ( exists == false ) {
            refx = run->getXcWF()[i];
            xwflist.append(refx);
          }
       } 
   }
   xnearfieldComboBox = new QComboBox(this);
   xfarfieldComboBox = new QComboBox(this);
   xforceComboBox = new QComboBox(this);
   qSort(xlist.begin(), xlist.end());
   QStringList xliststr;
   for ( int i = 0 ; i < xlist.size(); i++ ) {
      xliststr << QString::number(xlist.at(i),'g',8);
     QLOG_DEBUG() << " COMPARE " << xlist.at(i) << " AND " << xliststr.at(i);
   }
   qSort(xwflist.begin(), xwflist.end());
   QStringList xwfliststr;
   for ( int i = 0 ; i < xwflist.size(); i++ ) {
      xwfliststr << QString::number(xwflist.at(i),'g',8);
     QLOG_DEBUG() << " COMPARE " << xwflist.at(i) << " AND " << xwfliststr.at(i);
   }
   xforceComboBox->addItems(xliststr);
   if (options->getNproche() != 2) {
     xnearfieldComboBox->addItems(xliststr);
     xfarfieldComboBox->addItems(xliststr);
   }
   else if ( options->getNproche() == 2) {
     xnearfieldComboBox->addItems(xwfliststr);
     xfarfieldComboBox->addItems(xwfliststr);
   }
   
   // Plot Y list
   double refy; 
   QList<double> ylist;
   QList<double> ywflist;
   refy = run->getYc()[0];
   ylist.append(refy);
   for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
       bool exists = false;
       if ( run->getYc()[i] != refy) {
         for (int j = 0 ; j < ylist.size() ; j++) {
           if ( abs(run->getYc()[i] - ylist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
         }
         if ( exists == false ) {
           refy = run->getYc()[i];
           ylist.append(refy);
         }
       } 
   }
   refy = run->getYcWF()[0];
   ywflist.append(refy);
   for (int i = 0 ; i < nmaxpp ; i++) {
       bool exists = false;
       if ( run->getYcWF()[i] != refy) {
         for (int j = 0 ; j < ywflist.size() ; j++) {
           if ( abs(run->getYcWF()[i] - ywflist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
         }
         if ( exists == false ) {
           refy = run->getYcWF()[i];
           ywflist.append(refy);
         }
       } 
   }
   ynearfieldComboBox = new QComboBox(this);
   yfarfieldComboBox = new QComboBox(this);
   yforceComboBox = new QComboBox(this);
   qSort(ylist.begin(), ylist.end());
   QStringList yliststr;
   for ( int i = 0 ; i < ylist.size(); i++ )
      yliststr << QString::number(ylist.at(i),'g',8);
   qSort(ywflist.begin(), ywflist.end());
   QStringList ywfliststr;
   for ( int i = 0 ; i < ywflist.size(); i++ )
      ywfliststr << QString::number(ywflist.at(i),'g',8);
   yforceComboBox->addItems(yliststr);
   if (options->getNproche() != 2) {
     ynearfieldComboBox->addItems(yliststr);
     yfarfieldComboBox->addItems(yliststr);
   }
   else if (options->getNproche() == 2) {
     ynearfieldComboBox->addItems(ywfliststr);
     yfarfieldComboBox->addItems(ywfliststr);
   }
  
   // Plot Z list
    double refz;
   QList<double> zlist;
   QList<double> zwflist;
   refz = run->getZc()[0];
   zlist.append(refz);
   for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
       bool exists = false;
       if ( run->getZc()[i] != refz) {
         for (int j = 0 ; j < zlist.size() ; j++) {
           if ( abs(run->getZc()[i] - zlist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
          }
          if ( exists == false ) {
            refz = run->getZc()[i];
            zlist.append(refz);
          }
       } 
   }
   refz = run->getZcWF()[0];
   zwflist.append(refz);
   for (int i = 0 ; i < nmaxpp ; i++) {
       bool exists = false;
       if ( run->getZcWF()[i] != refz) {
         for (int j = 0 ; j < zwflist.size() ; j++) {
           if ( abs(run->getZcWF()[i] - zwflist.at(j)) < DELTA ) {
             exists = true;
             break;
           }
          }
          if ( exists == false ) {
            refz = run->getZcWF()[i];
            zwflist.append(refz);
          }
       } 
   }
   znearfieldComboBox = new QComboBox(this);
   zfarfieldComboBox = new QComboBox(this);
   zforceComboBox = new QComboBox(this);
   qSort(zlist.begin(), zlist.end());
   QStringList zliststr;
   for ( int i = 0 ; i < zlist.size(); i++ ) {
      zliststr << QString::number(zlist.at(i),'g',8);
     QLOG_DEBUG() << " COMPARE " << zlist.at(i) << " AND " << zliststr.at(i);
   }
   qSort(zwflist.begin(), zwflist.end());
   QStringList zwfliststr;
   for ( int i = 0 ; i < zwflist.size(); i++ ) {
      zwfliststr << QString::number(zwflist.at(i),'g',8);
     QLOG_DEBUG() << " COMPARE " << zwflist.at(i) << " AND " << zwfliststr.at(i);
   }
   zforceComboBox->addItems(zliststr);
   if (options->getNproche() != 2) {
     znearfieldComboBox->addItems(zliststr);
     zfarfieldComboBox->addItems(zliststr);
   }
   else if ( options->getNproche() == 2) {
     znearfieldComboBox->addItems(zwfliststr);
     zfarfieldComboBox->addItems(zwfliststr);
   }
   

   // Plot Selection
   // Near Field ComboBox
   nearfieldComboBox = new QComboBox(this);
   QStringList nearfieldList;
   nearfieldList << "Incident field";
    if ( options->getLocalfield() )
      nearfieldList << "Local field";
   if ( options->getMacroscopicfield() )
      nearfieldList << "Macroscopic field";
   nearfieldComboBox->addItems(nearfieldList);
   nearfieldtypeComboBox = new QComboBox(this);
   nearfieldtypeComboBox->addItems((QStringList() 
	<< "Modulus" << "x" << "y" << "z"));
   // Far Field ComboBox
   farfieldComboBox = new QComboBox(this);
   int ntypemic;
   ntypemic = options->getNtypemic();
   int nside;
   nside = options->getNside();
   int maxpp=run->getNmaxpp();
   QLOG_DEBUG () << "Beam:" << options->getNtypemic() << ntypemic;
   QLOG_DEBUG () << "Beam:" << options->getNside() << nside;
   if (ntypemic == 0) {
     if (nside == 0){
       farfieldComboBox->addItems((QStringList() << "Fourier plane: scattered field" <<  "Fourier plane: total field" << "Image plane: scattered field" << "Image plane: total field"));
     }
     else
       farfieldComboBox->addItems((QStringList() << "Fourier plane: scattered field"  << "Image plane: scattered field" ));
       }
   else if (ntypemic == 1 || ntypemic == 2) {
     if (nside == 0){
       farfieldComboBox->addItems((QStringList() <<  "Image plane: scattered field" << "Image plane: total field"));
     }
     else
       farfieldComboBox->addItems((QStringList() <<  "Image plane: scattered field" ));
   }
   farfieldtypeComboBox = new QComboBox(this);
   farfieldtypeComboBox->addItems((QStringList() 
   << "Intensity" << "Modulus" << "x" << "y" << "z"));
 
   // Force ComboBox
   forceComboBox = new QComboBox(this);
   QStringList forceList;
   if ( options->getOpticalforcedensity() )
      forceList << "force";
   if ( options->getOpticaltorquedensity() )
      forceList << "torque";
   forceComboBox->addItems(forceList);
   forcetypeComboBox = new QComboBox(this);
   forcetypeComboBox->addItems((QStringList() << "Modulus" ));

   ///////////////////////////////////////////////////////////////////////////
   // Dipolepsilon Layout
   QGridLayout *dipolepsilonLayout = new QGridLayout();
   dipolepsilonLayout->addWidget(dipoles3DButton,0,0);
   ///////////////////////////////////////////////////////////////////////////
   // Near Field Layout
   QGridLayout *nearfieldLayout = new QGridLayout();
   nearfieldLayout->setColumnStretch(1, 1);
   nearfieldLayout->setColumnStretch(3, 1);
   nearfieldLayout->setColumnStretch(4, 1);
   nearfieldLayout->setColumnStretch(5, 1);
   nearfieldLayout->addWidget(new QLabel("Cross section X:"),0,2);
   nearfieldLayout->addWidget(xnearfieldComboBox,0,3);
   nearfieldLayout->addWidget(plotxnearfieldButton,0,4);
   nearfieldLayout->addWidget(plotallxnearfieldButton,0,5);
   nearfieldLayout->addWidget(new QLabel("Field:"),0,0);
   nearfieldLayout->addWidget(nearfieldComboBox,0,1);
   nearfieldLayout->addWidget(new QLabel("Type:"),1,0);
   nearfieldLayout->addWidget(nearfieldtypeComboBox,1,1);
   nearfieldLayout->addWidget(new QLabel("Cross section Y:"),1,2);
   nearfieldLayout->addWidget(ynearfieldComboBox,1,3);
   nearfieldLayout->addWidget(plotynearfieldButton,1,4);
   nearfieldLayout->addWidget(plotallynearfieldButton,1,5);
   nearfieldLayout->addWidget(new QLabel("Cross section Z:"),2,2);
   nearfieldLayout->addWidget(znearfieldComboBox,2,3);
   nearfieldLayout->addWidget(plotznearfieldButton,2,4);
   nearfieldLayout->addWidget(plotallznearfieldButton,2,5);
   ///////////////////////////////////////////////////////////////////////////
   // Force Layout
   QGridLayout *forceLayout = new QGridLayout();
   forceLayout->setColumnStretch(1, 1);
   forceLayout->setColumnStretch(3, 1);
   forceLayout->setColumnStretch(4, 1);
   forceLayout->setColumnStretch(5, 1);
   if ( options->getOpticalforcedensity() || options->getOpticaltorquedensity() ) {
     forceLayout->addWidget(new QLabel("Cross section X:"),0,2);
     forceLayout->addWidget(xforceComboBox,0,3);
     forceLayout->addWidget(plotxforceButton,0,4);
     forceLayout->addWidget(plotallxforceButton,0,5);
     forceLayout->addWidget(new QLabel("Field:"),0,0);
     forceLayout->addWidget(forceComboBox,0,1);
     forceLayout->addWidget(new QLabel("Type:"),1,0);
     forceLayout->addWidget(forcetypeComboBox,1,1);
     forceLayout->addWidget(new QLabel("Cross section Y:"),1,2);
     forceLayout->addWidget(yforceComboBox,1,3);
     forceLayout->addWidget(plotyforceButton,1,4);
     forceLayout->addWidget(plotallyforceButton,1,5);
     forceLayout->addWidget(new QLabel("Cross section Z:"),2,2);
     forceLayout->addWidget(zforceComboBox,2,3);
     forceLayout->addWidget(plotzforceButton,2,4);
     forceLayout->addWidget(plotallzforceButton,2,5);
   }
   ///////////////////////////////////////////////////////////////////////////
   // Far Field Layout
   QGridLayout *farfieldLayout = new QGridLayout();
   if ( options->getMicroscopy() ) {
     farfieldLayout->setColumnStretch(1, 1);
     farfieldLayout->setColumnStretch(2, 1);
     farfieldLayout->addWidget(new QLabel("Microscopy:"),0,0);
     farfieldLayout->addWidget(farfieldComboBox,0,1);
     farfieldLayout->addWidget(new QLabel("Type:"),1,0);
     farfieldLayout->addWidget(farfieldtypeComboBox,1,1);
     farfieldLayout->addWidget(plotmicroscopyButton,0,2);
   }
   if ( options->getCrosssectionpoynting() ) {
     if ( options->getMicroscopy() )
       farfieldLayout->addWidget(poyntingButton,1,2);
     else 
       farfieldLayout->addWidget(poyntingButton,0,0);
   }

   QGroupBox *dipolepsilonGroupBox = new QGroupBox();
   dipolepsilonGroupBox->setTitle("Check discretization");
   dipolepsilonGroupBox->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 0px;}");
   dipolepsilonGroupBox->setLayout(dipolepsilonLayout);

   QGroupBox *nearfieldGroupBox = new QGroupBox();
   nearfieldGroupBox->setTitle("Near field investigation");
   nearfieldGroupBox->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 0px;}");
   nearfieldGroupBox->setLayout(nearfieldLayout);

   QGroupBox *farfieldGroupBox = new QGroupBox();
   farfieldGroupBox->setTitle("Far field  and microscopy investigation");
   farfieldGroupBox->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 0px;}");
   farfieldGroupBox->setLayout(farfieldLayout);

   QGroupBox *forceGroupBox = new QGroupBox();
   forceGroupBox->setTitle("Force investigation");
   forceGroupBox->setStyleSheet("QGroupBox{border:2px solid gray;border-radius:5px;margin-top: 1ex;} QGroupBox::title{subcontrol-origin: margin;subcontrol-position:top center;padding:0 0px;}");
   forceGroupBox->setLayout(forceLayout);

   QBoxLayout *utilitylayout = new QBoxLayout(QBoxLayout::TopToBottom);
   if ( options->getNearfield() || options->getFarfield() ||  
	options->getForce() || options->getMicroscopy() || options->getDipolepsilon() )
   utilitylayout->addWidget(dipolepsilonGroupBox); // Always present
   if ( options->getNearfield() )
     utilitylayout->addWidget(nearfieldGroupBox);
   if ( options->getFarfield() || options->getMicroscopy())
     utilitylayout->addWidget(farfieldGroupBox);
   if ( options->getOpticalforcedensity() || options->getOpticaltorquedensity() )
     utilitylayout->addWidget(forceGroupBox);
   // Central Layout
   centrallayout->addRow(utilitylayout);
}
void
RunWidget::plotxnearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   plotx(field, xnearfieldComboBox->currentText().toDouble());
}
void
RunWidget::plotynearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   ploty(field, ynearfieldComboBox->currentText().toDouble());
}
void
RunWidget::plotznearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   plotz(field, znearfieldComboBox->currentText().toDouble());
}
void
RunWidget::plotmicroscopy() {
   // Modulus case
   QVector<QwtPoint3D> *data;
   data = new QVector<QwtPoint3D>();
   // Phase case
   QVector<QwtPoint3D> *datapc;
   datapc = new QVector<QwtPoint3D>();
   QWidget *plotfinalwidget = new QWidget();
   QBoxLayout *plotlayout = new QBoxLayout(QBoxLayout::LeftToRight);
   QString type = farfieldtypeComboBox->currentText();
   QString field = farfieldComboBox->currentText();
   int line = 0;
   int col = 0;
   int nfft2d;
   nfft2d = options->getnfft2d();
   QString title,xtitle,ytitle;
   if ( type == "Intensity" ) {
     if ( field == "Image plane: scattered field" ) {
       title = "Image plane: scattered field";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], (
                          norm(run->getEimageX()[i]) +
                          norm(run->getEimageY()[i]) +
                          norm(run->getEimageZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Fourier plane: scattered field" ) {
       title = "Fourier plane: scattered field";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], (
              norm(run->getEfourierX()[i]) +
              norm(run->getEfourierY()[i]) +
              norm(run->getEfourierZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Fourier plane: total field" ) {
       title = "Fourier plane: total field";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], (
              norm(run->getEfourierincX()[i]) +
              norm(run->getEfourierincY()[i]) +
              norm(run->getEfourierincZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Image plane: total field" ) {
       title = "Image plane: total field";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], (
                            norm(run->getEimageincX()[i]) +
                            norm(run->getEimageincY()[i]) +
                            norm(run->getEimageincZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     plot = new PlotRaster(this, data, (int)nfft2d, title, xtitle, ytitle, options->getColors());
     plot->setObjectName("plotraster");
     plotlayout->addWidget(plot);
   }
   else if ( type == "Modulus" ) {
     if ( field == "Image plane: scattered field" ) {
       title = "Image plane: scattered field";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], sqrt(
                          norm(run->getEimageX()[i]) +
                          norm(run->getEimageY()[i]) +
                          norm(run->getEimageZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Fourier plane: scattered field" ) {
       title = "Fourier plane: scattered field";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], sqrt(
              norm(run->getEfourierX()[i]) +
              norm(run->getEfourierY()[i]) +
              norm(run->getEfourierZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Fourier plane: total field" ) {
       title = "Fourier plane: total field";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], sqrt(
              norm(run->getEfourierincX()[i]) +
              norm(run->getEfourierincY()[i]) +
              norm(run->getEfourierincZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     else if ( field == "Image plane: total field" ) {
       title = "Image plane: total field";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
	  data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], sqrt(
                            norm(run->getEimageincX()[i]) +
                            norm(run->getEimageincY()[i]) +
                            norm(run->getEimageincZ()[i]))));
          col++;
          if ( col == nfft2d ) {
             col = 0;
             line++;
          }
       }
     }
     plot = new PlotRaster(this, data, (int)nfft2d, title, xtitle, ytitle, options->getColors());
     plot->setObjectName("plotraster");
     plotlayout->addWidget(plot);
   }
   else { // plot Modulus and phase
     if ( field == "Image plane: scattered field" && type == "x" ) {
       title = "Image plane: scattered field X";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageX()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageX()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: scattered field" && type == "x" ) {
       title = "Fourier plane: scattered field X";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierX()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierX()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: total field" && type == "x" ) {
       title = "Fourier plane: total field X";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierincX()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierincX()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Image plane: total field" && type == "x" ) {
       title = "Image plane: total field X";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageincX()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageincX()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Image plane: scattered field" && type == "y" ) {
       title = "Image plane: scattered field Y";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageY()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageY()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: scattered field" && type == "y" ) {
       title = "Fourier plane: scattered field Y";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierY()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierY()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: total field" && type == "y" ) {
       title = "Fourier plane: total field Y";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierincY()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierincY()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Image plane: total field" && type == "y" ) {
       title = "Image plane: total field";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageincY()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageincY()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Image plane: scattered field" && type == "z" ) {
       title = "Image plane: scattered field Z";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageZ()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageZ()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: scattered field" && type == "z" ) {
       title = "Fourier plane: scattered field Z";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierZ()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierZ()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Fourier plane: total field" && type == "z" ) {
       title = "Fourier plane: total field Z";
       xtitle = "k_x/k_0";
       ytitle = "k_y/k_0";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], abs(run->getEfourierincZ()[i])));
         datapc->push_back(QwtPoint3D(run->getKXY()[line],run->getKXY()[col], arg(run->getEfourierincZ()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     else if ( field == "Image plane: total field" && type == "z" ) {
       title = "Image plane: total field Z";
       xtitle = "x(m)";
       ytitle = "y(m)";
       for ( int i = 0 ; i < nfft2d*nfft2d ; i++ ) {
         data->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], abs(run->getEimageincZ()[i])));
         datapc->push_back(QwtPoint3D(run->getXY()[line],run->getXY()[col], arg(run->getEimageincZ()[i])));
         col++;
         if ( col == nfft2d ) {
            col = 0;
            line++;
         }
       }
     }
     plot = new PlotRaster(this, data, (int)nfft2d, title + " Modulus", xtitle, ytitle, options->getColors());
     plot->setObjectName("plotraster");
     plotlayout->addWidget(plot);
     plotpc = new PlotRaster(this, datapc, (int)nfft2d, title + " phase", xtitle, ytitle, options->getColors());
     plotpc->setObjectName("plotraster");
     plotlayout->addWidget(plotpc);
   }
   plotfinalwidget->setLayout(plotlayout);
   plotwidget->setCurrentIndex(plotwidget->addTab(plotfinalwidget,"Microscopy"));
}

void
RunWidget::plotxforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   plotx(field, xforceComboBox->currentText().toDouble());
   QLOG_DEBUG() << "RunWidget::plotxforce> plotting" << field;
}
void
RunWidget::plotyforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   ploty(field, yforceComboBox->currentText().toDouble());
}
void
RunWidget::plotzforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   plotz(field, zforceComboBox->currentText().toDouble());
}
void
RunWidget::plotallxnearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   plotallx(field, xnearfieldComboBox);
}
void
RunWidget::plotallynearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   plotally(field, ynearfieldComboBox);
}
void
RunWidget::plotallznearfield() {
   QString type = nearfieldtypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = nearfieldComboBox->currentText();
   else
     field = nearfieldComboBox->currentText() + type;
   plotallz(field, znearfieldComboBox);
}
void
RunWidget::plotallxforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   plotallx(field, xforceComboBox);
}
void
RunWidget::plotallyforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   plotally(field, yforceComboBox);
}
void
RunWidget::plotallzforce() {
   QString type = forcetypeComboBox->currentText();
   QString field = "";
   if ( type == "Modulus")
     field = forceComboBox->currentText();
   else
     field = forceComboBox->currentText() + type;
   plotallz(field, zforceComboBox);
}
void
RunWidget::plotallx(QString field, QComboBox *xlist) {
    for (int i = 0 ; i < xlist->count(); i++)
       plotAxislist(field,
		run->getYc(),run->getZc(),run->getXc(),
		run->getYcWF(),run->getZcWF(),run->getXcWF(),
		run->getZc(),run->getZcWF(),xlist->itemText(i).toDouble(),
		":modulus(Y,Z)_", ":phase(Y,Z)_", " (Y,Z)_", "X:");
}
void
RunWidget::plotally(QString field, QComboBox *ylist) {
    for (int i = 0 ; i < ylist->count(); i++)
       plotAxislist(field,
		run->getXc(),run->getZc(),run->getYc(),
		run->getXcWF(),run->getZcWF(),run->getYcWF(),
		run->getZc(),run->getZcWF(),ylist->itemText(i).toDouble(),
		":modulus(X,Z)_", ":phase(X,Z)_", " (X,Z)_", "Y:");
}
void
RunWidget::plotallz(QString field, QComboBox *zlist) {
    for (int i = 0 ; i < zlist->count(); i++)
       plotAxislist(field,
		run->getXc(),run->getYc(),run->getZc(),
		run->getXcWF(),run->getYcWF(),run->getZcWF(),
		run->getYc(),run->getYcWF(),zlist->itemText(i).toDouble(),
		":modulus(X,Y)_", ":phase(X,Y)_", " (X,Y)_", "Z:");
}
void
RunWidget::plotx(QString field, double xvalue) {
    QLOG_DEBUG() << "RunWidget::plotx> " << field << " " << QString::number(xvalue);
    plotAxislist(field,
	run->getYc(),run->getZc(),run->getXc(),
	run->getYcWF(),run->getZcWF(),run->getXcWF(),
	run->getZc(),run->getZcWF(),xvalue,
	":modulus(Y,Z)_", ":phase(Y,Z)_", " (Y,Z)_", "X:");
}
void
RunWidget::ploty(QString field, double yvalue) {
    plotAxislist(field,
	run->getXc(),run->getZc(),run->getYc(),
	run->getXcWF(),run->getZcWF(),run->getYcWF(),
	run->getZc(),run->getZcWF(),yvalue,
	":modulus(X,Z)_", ":phase(X,Z)_", " (X,Z)_", "Y:");
}
void
RunWidget::plotz(QString field, double zvalue) {
    plotAxislist(field,
	run->getXc(),run->getYc(),run->getZc(),
	run->getXcWF(),run->getYcWF(),run->getZcWF(),
	run->getYc(),run->getYcWF(),zvalue,
	":modulus(X,Y)_", ":phase(X,Y)_", " (X,Y)_", "Z:");
}
void
RunWidget::dipoles3Dplot() {
  QVector<QwtPoint3D> *data;
  data = new QVector<QwtPoint3D>();
  QVector<double> *epsilon;
  epsilon = new QVector<double>();
  int nmax = options->getNxm() * options->getNym() * options->getNzm();
  QLOG_DEBUG() << " getObjectSubunits = " << run->getObjectSubunits();
  /*QFile filex("xc.txt");
  QFile filey("yc.txt");
  QFile filez("zc.txt");
  QFile fileps("epsilon.txt");
   if (!filex.open(QIODevice::WriteOnly | QIODevice::Text))
         return;
    if (!filey.open(QIODevice::WriteOnly | QIODevice::Text))
         return;
    if (!filez.open(QIODevice::WriteOnly | QIODevice::Text))
         return;
    if (!fileps.open(QIODevice::WriteOnly | QIODevice::Text))
         return;
  QTextStream outx(&filex);
  QTextStream outz(&filez);
  QTextStream outps(&fileps);*/
 
    for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
       data->push_back(QwtPoint3D(run->getXc()[i],run->getYc()[i],run->getZc()[i]));
       epsilon->push_back(abs(run->getEpsilonField()[i+4*nmax]));
       /*outx << "  " << QString::number(run->getXc()[i],'g',16) << endl;
       outy <<"  " << QString::number(run->getYc()[i],'g',16)<< endl;
       outz << "  " <<QString::number(run->getZc()[i],'g',16)<< endl;
       outps << "  " <<QString::number(epsilon->at(i),'g',16)<< endl;
       */
     }
 
  
  /*filex.close();
  filey.close();
   filez.close();
  fileps.close();*/
  Plot *plot = new Plot(data,epsilon,"Dipoles",0);
  plot->setObjectName("plot3d");
  plotwidget->setCurrentIndex(plotwidget->addTab(plot,"Dipoles with epsilon:"));
}
void
RunWidget::poyntingplot() {
  QVector<QwtPoint3D> *data;
  data = new QVector<QwtPoint3D>();
  QVector<double> *Modulus;
  Modulus = new QVector<double>();
  ntheta = options->getNtheta();
  nphi = options->getNphi();
  for (int i = 0 ; i < (ntheta+1) * nphi ; i++) {
    double r = run->getPoyntingField()[i];
    double theta = run->getThetaField()[i] * PI / 180;
    double phi = run->getPhiField()[i] * PI / 180;
    Modulus->push_back(r);
    QLOG_DEBUG() << " poynting = " << run->getPoyntingField()[i]
		 << " theta = " << run->getThetaField()[i]
		 << " phi = " << run->getPhiField()[i];
    data->push_back(QwtPoint3D(r*sin(theta)*cos(phi),
			       r*sin(theta)*sin(phi),
			       r*cos(theta)));
    QLOG_DEBUG() << " x = " << r*sin(theta)*cos(phi) 
		 << " y = " << r*sin(theta)*sin(phi)
		 << " z = " << r*cos(theta);
  }
  QWidget *plotfinalwidget = new QWidget();
  Plot *plot = new Plot(data,Modulus,"Poynting",nphi);
  plot->setObjectName("plot3d");
  plotwidget->setCurrentIndex(plotwidget->addTab(plot,"Poynting:"));
}
void
RunWidget::cleanupPlots() {
   for (int i = plotwidget->count() - 1; i >= 0 ; i--) {
     QWidget * tmpwidget = plotwidget->widget(i);
     Plot *plot = NULL;
     plot = tmpwidget->findChild<Plot *>("plot3d");
     if ( plot != NULL )
        delete plot;
     PlotRaster *raster = NULL;
     raster = tmpwidget->findChild<PlotRaster *>("plotraster");
     if ( raster != NULL )
        delete raster;
     PlotVector *plotvect = NULL;
     plotvect = tmpwidget->findChild<PlotVector *>("plotvector");
     if ( plotvect != NULL )
        delete plotvect;
     plotwidget->removeTab(i);
     delete tmpwidget;
   }
}
void
RunWidget::colorPlots() {
  QVector<QColor> *colors;
  colors = options->getColors();
  QColor colorlow = QColorDialog::getColor(colors->at(0),NULL,
		    QString("Color Low"),QColorDialog::DontUseNativeDialog);
  if (colorlow.isValid())
    colors->replace(0,colorlow);
  QColor colorhigh = QColorDialog::getColor(colors->at(1),NULL,
		     QString("Color High"),QColorDialog::DontUseNativeDialog);
  if (colorhigh.isValid()) 
    colors->replace(1,colorhigh);
  QLOG_DEBUG() << "Set new colors interval low:" << options->getColors()->at(0) 
		<< " high:" << options->getColors()->at(1);
}
void
RunWidget::savePlots() { 
 
  QWidget *currwidget = plotwidget->currentWidget();
  if (currwidget == NULL) return;
  if ( currwidget->objectName() == "plot3d" ) {
      QLOG_DEBUG() <<"found plot object"<< currwidget->objectName();
      Plot *plot = (Plot*)currwidget;
      QString defaultName = QDir::currentPath() + QDir::separator () + plot->getTitle() + ".ps";
      QString filename = QFileDialog::getSaveFileName(this,tr("Save file"),
			 defaultName,
        		 tr("Postscript files (*.ps);;PDF files (*.pdf);;JPEG files (*.jpg)"),
			 0,QFileDialog::DontUseNativeDialog);
      QStringList fileext = filename.split(".");
       QLOG_DEBUG() << " Saving " << defaultName;
      plot->save(filename,fileext.at(fileext.size() - 1).toUpper());
  }
  else  {
     QList<QWidget*> plotlist;
     plotlist = currwidget->findChildren<QWidget*>();
     for (int j = 0; j < plotlist.size(); j++) {
       QString name = plotlist.at(j)->objectName();
       QLOG_DEBUG() << "RunWidget::savePlots> " << name; 
       if ( name == "plotraster" || name == "plotvector" ) {
         QLOG_DEBUG() <<"found plot object"<< j <<" = "<< name;
         QwtPlot *plot = (QwtPlot*)plotlist.at(j);
         QString defaultName = QDir::currentPath() + QDir::separator () + plot->title().text() + ".ps";
         QLOG_DEBUG() << " Saving " << defaultName;
         QString filename = QFileDialog::getSaveFileName(this,tr("Save file"),
		   	    defaultName,
        		    tr("Postscript files (*.ps);;PDF files (*.pdf);;JPEG files (*.jpg)"),
			    0,QFileDialog::DontUseNativeDialog);
         QwtPlotRenderer renderer;
         renderer.setLayoutFlags(QwtPlotRenderer::FrameWithScales);
         renderer.renderDocument( plot, filename, QSizeF(300,200), 85 );
       }
     }
   }
}
void
RunWidget::printPlots() {
    
  QWidget *currwidget = plotwidget->currentWidget();
  if (currwidget == NULL) return;
  if ( currwidget->objectName() == "plot3d" ) {
    QLOG_DEBUG() <<"found plot object"<< currwidget->objectName();
    Plot *plot = (Plot*)currwidget;
    QPrinter *printer = new QPrinter;
    printer->setDocName ( plot->getTitle() );
    QPrintDialog *printDialog = new QPrintDialog(printer, this);
    if (printDialog->exec() == QDialog::Accepted) {
	    QPainter p(printer);
	    QPixmap pixmap(plot->size());
	    plot->render(&pixmap);
	    p.scale(0.80,0.80);
	    p.drawPixmap(0, 0, pixmap);
    }
  }
  else  {
    QList<QWidget*> plotlist;
    plotlist = currwidget->findChildren<QWidget*>();
    for (int j = 0; j < plotlist.size(); j++) {
      QString name = plotlist.at(j)->objectName();
      QLOG_DEBUG() << "RunWidget::savePlots> " << name; 
      if ( name == "plotraster" || name == "plotvector" ) {
        QLOG_DEBUG() <<"found plot object"<< j <<" = "<< name;
        QwtPlot *plot = (QwtPlot*)plotlist.at(j);

        QPrinter printer( QPrinter::HighResolution );
        //printer.setPaperSize(QSizeF(150, 100), QPrinter::Millimeter);
        QString docName = plot->title().text();
        if ( !docName.isEmpty() ) {
          docName.replace ( QRegExp ( QString::fromLatin1 ( "\n" ) ), tr ( " -- " ) );
          printer.setDocName ( docName );
        }
        printer.setCreator( "Plot" );
        printer.setOrientation( QPrinter::Portrait );
        QPrintDialog dialog( &printer );
        if ( dialog.exec() ) {
          QwtPlotRenderer renderer;
          if ( printer.colorMode() == QPrinter::GrayScale ) {
            renderer.setDiscardFlag( QwtPlotRenderer::DiscardBackground );
            renderer.setDiscardFlag( QwtPlotRenderer::DiscardCanvasBackground );
            renderer.setDiscardFlag( QwtPlotRenderer::DiscardCanvasFrame );
            renderer.setLayoutFlag( QwtPlotRenderer::FrameWithScales );
          }
          renderer.renderTo( plot, printer );
       }
     }
   }
  }
}
void 
RunWidget::plotAxislist(QString field, double *X, double *Y, double *Z, 
			double *XWF, double *YWF, double *ZWF, double *REFCOL, 
			double *REFCOLWF, double refaxis,
			QString title1, QString title2, QString title3, QString title4) {
   QVector<int> *num_colp;
   num_colp = new QVector<int>();
   int colp = 0, ref_colp = 0;
   double refp = 0;
   // Modulus case
   QVector<QwtPoint3D> *data;
   data = new QVector<QwtPoint3D>();
   // Vector field case
   QVector<QwtPointPolar> *datap;
   datap = new QVector<QwtPointPolar>();
   // Phase case
   QVector<QwtPoint3D> *datapc;
   datapc = new QVector<QwtPoint3D>();
   QWidget *plotfinalwidget = new QWidget();
   QBoxLayout *plotlayout = new QBoxLayout(QBoxLayout::LeftToRight);
   int cnt = 0;
   int pos = 0;
   QString xtitle,ytitle;
   if ( title4 == "X:" ) {
     xtitle = "y(m)";
     ytitle = "z(m)";
   }
   else if ( title4 == "Y:" ) {
     xtitle = "x(m)";
     ytitle = "z(m)";
   }
   else if ( title4 == "Z:" ) {
     xtitle = "x(m)";
     ytitle = "y(m)";
   }
   QLOG_DEBUG() << "RunWidget::plotAxislist> NPROCHE " << options->getNproche();
       for (int i = 0 ; i < run->getObjectSubunits() ; i++) {
         if (QString::number(Z[i],'g',8) == QString::number(refaxis,'g',8)) {
           if (pos == 0) {
             refp = REFCOL[i];
             pos = -1;
           }
	   cnt++;
	   if (field == "Incident field" && options->getNproche() != 2)
	     data->push_back(QwtPoint3D(X[i],Y[i], run->getIncidentField()[i]));
	   else if (field == "Local field" && options->getNproche() != 2)
	      data->push_back(QwtPoint3D(X[i],Y[i], run->getLocalField()[i]));
	   else if (field == "Macroscopic field" && options->getNproche() != 2)
	      data->push_back(QwtPoint3D(X[i],Y[i], run->getMacroscopicField()[i]));
	   else if (field == "force") {
	      data->push_back(QwtPoint3D(X[i],Y[i], run->getForceX()[i]));
              datap->push_back(QwtPointPolar(atan2( run->getForceY()[i],
 						    run->getForceX()[i]),
				sqrt(run->getForceX()[i] * run->getForceX()[i] +
 				run->getForceY()[i] * run->getForceY()[i]) ));
           }
	   else if (field == "torque") {
	      data->push_back(QwtPoint3D(X[i],Y[i], run->getTorqueX()[i]));
              datap->push_back(QwtPointPolar(atan2( run->getTorqueY()[i],
						    run->getTorqueX()[i]),
		    		sqrt(run->getTorqueX()[i] * run->getTorqueX()[i] +
                    		run->getTorqueY()[i] * run->getTorqueY()[i])));
           }
           else if (field == "Incident fieldx" && options->getNproche() != 2) {
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getIncidentFieldX()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getIncidentFieldX()[i])));
           }
           else if (field == "Local fieldx" && options->getNproche() != 2) {
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getLocalFieldX()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getLocalFieldX()[i])));
           }
           else if (field == "Macroscopic fieldx" && options->getNproche() != 2){
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getMacroscopicFieldX()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getMacroscopicFieldX()[i])));
           }
           else if (field == "Incident fieldy" && options->getNproche() != 2){
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getIncidentFieldY()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getIncidentFieldY()[i])));
           }
           else if (field == "Local fieldy" && options->getNproche() != 2){
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getLocalFieldY()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getLocalFieldY()[i])));
           }
           else if (field == "Macroscopic fieldy" && options->getNproche() != 2){
	      data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getMacroscopicFieldY()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getMacroscopicFieldY()[i])));
           }
           else if (field == "Incident fieldz" && options->getNproche() != 2){
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getIncidentFieldZ()[i])));
             datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getIncidentFieldZ()[i])));
           }
           else if (field == "Local fieldz" && options->getNproche() != 2){
              data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getLocalFieldZ()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getLocalFieldZ()[i])));
           }
           else if (field == "Macroscopic fieldz" && options->getNproche() != 2){
	      data->push_back(QwtPoint3D(X[i],Y[i], abs(run->getMacroscopicFieldZ()[i])));
              datapc->push_back(QwtPoint3D(X[i],Y[i], arg(run->getMacroscopicFieldZ()[i])));
           }
	   if (abs(refp - REFCOL[i]) < DELTA )
              colp++;
           else {
              num_colp->push_back(colp);
              QLOG_DEBUG() << "Column number:" << colp;
              if (ref_colp < colp) {
                ref_colp = colp;
              }
              refp = REFCOL[i];
              colp = 1;
           }
         }
     }
     if (options->getNproche() == 2 && field != "force" 
				    && field != "torque"
				    ) {
        num_colp->clear();
        cnt = pos = colp = ref_colp = 0;
        int nmax = options->getNxm() * options->getNym() * options->getNzm();
	int nmaxpp=run->getNmaxpp();
        for (int i = 0 ; i < nmaxpp ; i++) {
          if (QString::number(ZWF[i],'g',8) == QString::number(refaxis,'g',8)) {
           if (pos == 0) {
             refp = REFCOLWF[i];
             pos = -1;
           }
	   cnt++;
           if (field == "Incident field")
	      data->push_back(QwtPoint3D(XWF[i],YWF[i], run->getIncidentField()[i]));
	   else if (field == "Local field")
	     data->push_back(QwtPoint3D(XWF[i],YWF[i], run->getLocalField()[i]));
	   else if (field == "Macroscopic field")
	     data->push_back(QwtPoint3D(XWF[i],YWF[i], run->getMacroscopicField()[i]));
           else if (field == "Incident fieldx") {
              data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getIncidentFieldX()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getIncidentFieldX()[i])));
           }
           else if (field == "Local fieldx") {
              data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getLocalFieldX()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getLocalFieldX()[i])));
           }
           else if (field == "Macroscopic fieldx"){
              data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getMacroscopicFieldX()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getMacroscopicFieldX()[i])));
           }
           else if (field == "Incident fieldy"){
              data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getIncidentFieldY()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getIncidentFieldY()[i])));
           }
           else if (field == "Local fieldy"){
              data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getLocalFieldY()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getLocalFieldY()[i])));
           }
           else if (field == "Macroscopic fieldy"){
	      data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getMacroscopicFieldY()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getMacroscopicFieldY()[i])));
           }
           else if (field == "Incident fieldz"){
              data->push_back(QwtPoint3D(XWF[i],YWF[i],abs(run->getIncidentFieldZ()[i])));
             datapc->push_back(QwtPoint3D(XWF[i],YWF[i],arg(run->getIncidentFieldZ()[i])));
           }
           else if (field == "Local fieldz"){
              data->push_back(QwtPoint3D(XWF[i],YWF[i],abs(run->getLocalFieldZ()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i],arg(run->getLocalFieldZ()[i])));
           }
           else if (field == "Macroscopic fieldz"){
	      data->push_back(QwtPoint3D(XWF[i],YWF[i], abs(run->getMacroscopicFieldZ()[i])));
              datapc->push_back(QwtPoint3D(XWF[i],YWF[i], arg(run->getMacroscopicFieldZ()[i])));
           }
	   if (abs(refp - REFCOLWF[i]) < DELTA )
              colp++;
           else {
              num_colp->push_back(colp);
              QLOG_DEBUG() << "Column number:" << colp;
              if (ref_colp < colp) {
                ref_colp = colp;
              }
              refp = REFCOLWF[i];
              colp = 1;
           }
         }
       }
     }
     num_colp->push_back(colp);
     QLOG_DEBUG() << "Column number:" << colp;
     QLOG_DEBUG() << " Points found:" << cnt;
     QLOG_DEBUG() << " Max number of points found in line:" << ref_colp;
     QLOG_DEBUG() << " Number of Points found:" << data->size();
     QLOG_DEBUG() << " Number of lines:" << num_colp->size();
     int num_col3;
     if (field != "force" && field != "torque") {
      // Insert 0 values to complete matrix where necessary
      // First calculate min and max for x,y,z
      double minx=1e308,maxx=-1e308;
      double miny=1e308,maxy=-1e308;
      double minz=1e308,maxz=-1e308;
      for (int i = 0 ; i < data->size(); i++) {
        QwtPoint3D point = data->at(i);
        double tmpx = point.x();
        double tmpy = point.y();
        double tmpz = point.z();
        if (tmpx > maxx) maxx = tmpx;
        if (tmpx < minx) minx = tmpx;
        if (tmpy > maxy) maxy = tmpy;
        if (tmpy < miny) miny = tmpy;
        if (tmpz > maxz) maxz = tmpz;
        if (tmpz < minz) minz = tmpz;
      }
       QLOG_DEBUG() << "Min x:" << QString::number(minx);
       QLOG_DEBUG() << "Max x:" << QString::number(maxx);
       QLOG_DEBUG() << "Min y:" << QString::number(miny);
       QLOG_DEBUG() << "Max y:" << QString::number(maxy);
       QLOG_DEBUG() << "Min field:" << QString::number(minz);
       QLOG_DEBUG() << "Max field:" << QString::number(maxz);
       // Calculate xstepsize and num_col
       for (int k = 0 ; k < num_colp->size() ; k++) {
         if ( num_colp->at(k) > 1 ) {
            int poscalc = 0;
            for (int kk = 0 ; kk < k; kk++)
              poscalc+=num_colp->at(kk);
            xstepsize =  fabs(data->at(poscalc).x() - data->at(poscalc+1).x());
            num_col = (maxx - minx) / xstepsize + 1;
            QLOG_DEBUG() << "StepSize X:" << QString::number(xstepsize) 
			<< " Number of colums:" << QString::number(num_col);
            break;
         }
       }
      QLOG_DEBUG() << "1:StepSize X:" << QString::number(xstepsize) 
			<< " Number of colums:" << num_col;
      // Complete matrix data for raster
      pos = 0;
      QLOG_DEBUG() << "num_colp->size()= = " << num_colp->size();
      for (int k = 0 ; k < num_colp->size() ; k++) {
       double curx = minx;
       int tmppos = pos;
       
       QLOG_DEBUG() << "tmppos = " << tmppos << "num_colp(" << k << ")=" << num_colp->at(k);
       double tmp = data->at(tmppos).x();
       while ( curx < tmp - xstepsize / 10 ) {
          QLOG_DEBUG() << "insert right curx = " << QString::number(curx) 
			<< "tmp = " << QString::number(tmp)
			<< " at pos = " << pos << " data size = " << data->size();
          data->insert(pos,QwtPoint3D(curx,data->at(tmppos).y(),-100));
          if (field.contains("x") || field.contains("y") || field.contains("z")) 
            datapc->insert(pos,QwtPoint3D(curx,data->at(tmppos).y(),-100));
          curx = curx + xstepsize;
          QLOG_DEBUG() << "incremented curx = " << QString::number(curx) ;
          pos++;
       }
       curx = data->at(pos+num_colp->at(k) - 1).x() + xstepsize;
       QLOG_DEBUG() << "second loop curx = " << QString::number(curx);
       pos= pos + num_colp->at(k);
       while ( curx < maxx + xstepsize / 10 ) {
          QLOG_DEBUG() << "insert left curx = " << QString::number(curx) 
			<< " at pos = " << pos;
          data->insert(pos,QwtPoint3D(curx,data->at(tmppos).y(),-100));
          if (field.contains("x") || field.contains("y") || field.contains("z")) 
            datapc->insert(pos,QwtPoint3D(curx,data->at(tmppos).y(),-100));
          curx = curx + xstepsize;
          pos++;
       }
       for (int ll=0;ll<pos;ll++)
        QLOG_DEBUG() << "newdata(" << ll << ")=" << data->at(ll).x();
      }
     }
     QLOG_DEBUG() << " Total Number of Points:" << data->size();
     if (field != "force" && field != "torque") {
       
        QLOG_DEBUG() << "Before Calling PlotRaster StepSize X:" << QString::number(xstepsize) 
			<< " Number of colums:" << (int)qRound(num_col);
        plot = new PlotRaster(this, data, (int) (int)qRound(num_col), field+title1+title4+QString::number(refaxis),
			      xtitle, ytitle, options->getColors());
        plot->setObjectName("plotraster");
        plotlayout->addWidget(plot);
        if (field.contains("x") || field.contains("y") || field.contains("z")) {
          plotpc = new PlotRaster(this, datapc, (int)qRound(num_col), field+title2+title4+QString::number(refaxis),
			          xtitle, ytitle, options->getColors());
          plotpc->setObjectName("plotraster");
          plotlayout->addWidget(plotpc);
        }
     }
     else {
        // absalize datap radius
        double max = 0;
        for (int i = 0 ; i < datap->size(); i++)
         if ( max < datap->at(i).radius()) max = datap->at(i).radius();
        for (int i = 0 ; i < datap->size(); i++)
          datap->replace(i, QwtPointPolar(datap->at(i).azimuth(),datap->at(i).radius()/max));
        plotv = new PlotVector(this, data, datap, ref_colp, field+title3+title4+QString::number(refaxis), xtitle, ytitle);
        plotv->setObjectName("plotvector");
        plotlayout->addWidget(plotv);
     }
     plotfinalwidget->setLayout(plotlayout);
     plotwidget->setCurrentIndex(plotwidget->addTab(plotfinalwidget,title4 + QString::number(refaxis)));
   
}
