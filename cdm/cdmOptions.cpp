#include "cdmOptions.h"

Options::Options()
{
  colors = new QVector<QColor>();
  colors->push_back(Qt::blue);
  colors->push_back(Qt::yellow);

  beamList = (QStringList() << "Linear plane wave" <<  "Circular plane wave"
	                    << "Multiplane wave" << "Antenna"
			    << "Circular Gaussian" << "Linear Gaussian"
	                    << "Circular Gaussian (FFT)" << "Linear Gaussian (FFT)" 
			    << "Circular Gaussian (para)" << "Linear Gaussian (para)" 
		            << "Arbitrary wave (file)");

  objectList = (QStringList() << "sphere" << "inhomogeneous sphere" << "random spheres (length)"<< "random spheres (meshsize)" << "cube"
		<<  "cuboid (length)" << "cuboid (meshsize)"  << "inhomogeneous cuboid (length)" << "inhomogeneous cuboid (meshsize)"
		<< "ellipsoid" << "multiple spheres" 
		<<  "cylinder" << "concentric spheres" << "arbitrary");

  anisotropyList = (QStringList() << "iso" << "ani");

  materialList = (QStringList() << "xx" << "Au" << "Ag" << "Al" << "Ir" << "Pt" << "Cu" 
				<< "Glass" << "Water" << "W" << "Si" << "Mo" << "Ni" 
				<< "Os" << "Rh");

  methodeitList = (QStringList() << "GPBICG1" << "GPBICG2" << "GPBICGsafe"<< "GPBICGplus" << "GPBICGAR1" << "GPBICGAR2"
				 << "QMRCLA" << "TFQMR" << "CG" << "BICGSTAB" << "QMRBICGSTAB1" 
 		                 << "QMRBICGSTAB2" << "CORS" << "GPBICOR" << "BICGSTARPLUS" );

  polarizabilityList = (QStringList() << "RR" <<  "GB" << "LA" << "LR" << "CM" << "PS");

  quadList = (QStringList() << "0" << "1" <<  "2" << "3" << "4" << "5");

  nfft2dList = (QStringList() << "64" << "128" << "256" <<  "512"  << "768" << "1024" << "1536" << "2048" << "3072" << "4096" << "8192" << "16384");
  
  nrigList = (QStringList() << "rigorous" << "renormalized Born" << "Born"
	                    << "Born order 1" << "renormalized Rytov" << "Rytov"
	                    << "Beam propagation" << "renormalized Beam propagation" );
#ifdef USE_HDF5
  nmatlabList = (QStringList() << "Save in ascii file" << "Do not save file" << "Save in HDF5 file");
#else
  nmatlabList = (QStringList() << "Save in ascii file" << "Do not save file");
#endif
  ntypemicList = (QStringList() << "Holographic" << "Brightfield" << "Darkfield & phase");
  nsideList = (QStringList() << "kz>0 (transmission)" << "kz<0 (reflexion)");
    
  rangeofstudyList = (QStringList() << "Object" << "Cube around object" << "Wide field");

  initOptions();
}
Options::~Options() {
  QLOG_DEBUG ( ) << "Deleting Options";
}
void Options::initOptions(){ 
  sphereradius.resize(MAX_OBJECT_NUMBER);
  sphereradius.fill(0,MAX_OBJECT_NUMBER);
  positionx.resize(MAX_OBJECT_NUMBER);
  positionx.fill(0,MAX_OBJECT_NUMBER);
  for (int i = 0 ; i < MAX_OBJECT_NUMBER; i++)
   QLOG_DEBUG () << "positionx[" << i << "]=" << positionx.at(i);
  positiony.resize(MAX_OBJECT_NUMBER);
  positiony.fill(0,MAX_OBJECT_NUMBER);
  positionz.resize(MAX_OBJECT_NUMBER);
  positionz.fill(0,MAX_OBJECT_NUMBER);
  epsilon.resize(MAX_OBJECT_NUMBER);
  epsilon.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon11.resize(MAX_OBJECT_NUMBER);
  epsilon11.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  for (int i = 0 ; i < MAX_OBJECT_NUMBER; i++)
   QLOG_DEBUG () << "epsilon11[" << i << "]=" << real(epsilon11.at(i)) << " " << imag(epsilon11.at(i));
  epsilon12.resize(MAX_OBJECT_NUMBER);
  epsilon12.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon13.resize(MAX_OBJECT_NUMBER);
  epsilon13.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon21.resize(MAX_OBJECT_NUMBER);
  epsilon21.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon22.resize(MAX_OBJECT_NUMBER);
  epsilon22.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon23.resize(MAX_OBJECT_NUMBER);
  epsilon23.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon31.resize(MAX_OBJECT_NUMBER);
  epsilon31.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon32.resize(MAX_OBJECT_NUMBER);
  epsilon32.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  epsilon33.resize(MAX_OBJECT_NUMBER);
  epsilon33.fill(dcmplx(0,0),MAX_OBJECT_NUMBER);
  material.resize(MAX_OBJECT_NUMBER);
  material.fill("xx",MAX_OBJECT_NUMBER);
  polarizationTM = 0;
  polarizationTE = 0;
  polarizationRL = 0;
  incidenceangle_theta_z = 0;
  incidenceangle_phi_x = 0;
  P0 = 0;
  W0 = 0;
  demiaxea = 0;
  demiaxeb = 0;
  demiaxec = 0;
  thetaobj = 0;
  phiobj = 0;
  psiobj = 0;
  hauteur = 0;
  density = 0;
  cubeside = 0;
  cubesidex = 0;
  cubesidey = 0;
  cubesidez = 0;
  objectnumber = 1;
  thetam.resize(MAX_WAVEMULTI_NUMBER);
  thetam.fill(0,MAX_WAVEMULTI_NUMBER);
  phim.resize(MAX_WAVEMULTI_NUMBER);
  phim.fill(0,MAX_WAVEMULTI_NUMBER);
  ppm.resize(MAX_WAVEMULTI_NUMBER);
  ppm.fill(0,MAX_WAVEMULTI_NUMBER);
  ssm.resize(MAX_WAVEMULTI_NUMBER);
  ssm.fill(0,MAX_WAVEMULTI_NUMBER);
  E0m.resize(MAX_WAVEMULTI_NUMBER);
  E0m.fill(dcmplx(0,0),MAX_WAVEMULTI_NUMBER);
  wavemultinumber = 1;
  nxm = nym = nzm = ntheta = nphi = 0;
  na = nainc = 0.9;
  ntypemic= 0;
  nside=0;
  nxmp = nymp = nzmp = 0;
  nxx = nyy = nzz = 0;
  meshsize = 50;
  xgaus = ygaus = zgaus = 0;
  dipolepsilon = 0;
  farfield = 0;
  nearfield = 0;
  force = 0;
  localfield = 0;
  macroscopicfield = 0;
  crosssection = 0;
  crosssectionpoynting = 0;
  nenergie = 0;
  quickdiffract = 0;
  nrig = 0;
  nread = 0;
  nmatlab = 0;
  advancedinterface = 0;
  microscopy = 0;
  microscopyFFT = 0;
  opticalforce = 0;
  opticalforcedensity = 0;
  opticaltorque = 0;
  opticaltorquedensity = 0;
  nproche = 0;
  beam = "Linear plane wave";
  object = "sphere";
  methodeit = "GPBICG1";
  fileh5 = "ifdda.h5";
  polarizability = "RR";
  quad = 0;
  nfft2d = 128;
  discretization = 0;
  wavelength = 0;
  gross=0;
  zlens=0;
}
void Options::initDb() {
  dbpath = "options.db3";
  QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE",dbpath);
  db.setDatabaseName(dbpath);  
  if ( !db.open() ) {
    QLOG_ERROR() << "Options::initDb::No Db open!";
    return;
  }
  // Create options_tbl table
  QSqlQuery query(QSqlDatabase::database(dbpath));
  query.exec("create table options_tbl "
	     "(name varchar(255) not null primary key, "
	     "description varchar(255), "
	     "date varchar(255), "
	     "wavelength double, "
             "P0 double, "
             "W0 double, "
	     "beam varchar(255), "
             "wavemultinumber integer, "
	     "object varchar(255), "
             "objectnumber integer, "
	     "anisotropy varchar(255), "
	     "material varchar(255), "
	     "discretization integer, "
	     "tolerance double)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("insert into options_tbl values ('new','',datetime('now'),632.8,1,6328,'Linear plane wave',1,'sphere',1,'iso','xx',10,1e-4)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("create table beam_tbl "
	     "(parent varchar(255) not null primary key, "
             "beam varchar(255), "
             "beamfile varchar(255), "
	     "incidenceangle_theta_z double, "
	     "incidenceangle_phi_x double, "
	     "polarizationTM double, "
	     "polarizationTE double, "
	     "polarizationRL double, "
             "xgaus double, "
	     "ygaus double, "
	     "zgaus double)");
  query.exec("insert into beam_tbl values ('new','','',0,0,0,0,0,0,0,0)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("create table object_tbl "
	     "(parent varchar(255) not null, "
             "objectnum integer not null, "
             "object varchar(255), "
             "objectfile varchar(255), "
             "spherecoherencelength double, "
             "spherestandardev double, "
             "sphereseed int, "
	     "sphereradius double, "
	     "density double, "
	     "cubeside double, "
	     "cubesidex double, "
	     "cubesidey double, "
	     "cubesidez double, "
	     "positionx double, "
	     "positiony double, "
	     "positionz double, "
             "demiaxea double, "
	     "demiaxeb double, "
	     "demiaxec double, "
             "thetaobj double, "
	     "phiobj double, "
	     "psiobj double, "
             "hauteur double, "
             "primary key (parent, objectnum))");
  query.exec("insert into object_tbl values ('new',0,'sphere','',0,0,0,100,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("create table epsilon_tbl "
	     "(parent varchar(255) not null, "
             "objectnum integer not null, "
             "material varchar(255), "
   	     "anisotropy varchar(255), "
	     "epsilonr double, "
	     "epsiloni double, "
	     "epsilon11r double, "
 	     "epsilon11i double, "
	     "epsilon12r double, "
	     "epsilon12i double, "
	     "epsilon13r double, "
	     "epsilon13i double, "
	     "epsilon21r double, "
	     "epsilon21i double, "
	     "epsilon22r double, "
	     "epsilon22i double, "
	     "epsilon23r double, "
	     "epsilon23i double, "
	     "epsilon31r double, "
	     "epsilon31i double, "
	     "epsilon32r double, "
             "epsilon32i double, "
             "epsilon33r double, "
	     "epsilon33i double, "
             "primary key (parent, objectnum))");
  query.exec("insert into epsilon_tbl values ('new',0,'xx','iso',1.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("create table wavemulti_tbl "
	     "(parent varchar(255) not null, "
             "wavemultinum integer not null, "
	     "thetam double, "
	     "phim double, "
	     "ppm double, "
 	     "ssm double, "
	     "E0mr double, "
	     "E0mi double, "
             "primary key (parent, wavemultinum))");
  query.exec("insert into wavemulti_tbl values ('new',0,0,0,0,0,0,0)");
  QLOG_DEBUG () << query.lastError().text();
  query.exec("create table run_tbl "
	     "(parent varchar(255) not null primary key, "
             "methode varchar(255), "
             "polarizability varchar(255), "
             "quad integer, "
             "nrig integer, "
             "nread integer, "
             "filereread varchar(255), "
	     "nmatlab integer, "
             "fileh5 varchar(255), "
             "dipolepsilon integer, "
             "farfield integer, "
	     "nearfield integer, "
             "force integer, "
	     "localfield integer, "
	     "macroscopicfield integer, "
             "nenergie integer, "
	     "crosssection integer, "
	     "crosssectionpoynting integer, "
             "quickdiffract integer, "
             "microscopy integer, "
             "microscopyFFT integer, "
	     "opticalforce integer, "
	     "opticalforcedensity integer, "
   	     "opticaltorque integer, "
	     "opticaltorquedensity integer, "
	     "nproche integer, "
             "nxx integer, "
             "nyy integer, "
             "nzz integer, "
             "nxmp integer, "
             "nymp integer, "
             "nzmp integer, "
	     "ntheta integer, "
	     "nphi integer, "
             "na double, "
             "nainc double, "	     
             "gross double, "
	     "zlens double, "
	     "ntypemic integer,"
	     "nside integer,"
	     "meshsize double, "
	     "nfft2d integer, "
	     "advancedinterface integer) ");
  query.exec("insert into run_tbl values ('new','GPBICG1','RR',0,0,0,'',0,'ifdda.h5',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,72,0.9,0.9,100,0,0,0,25,128,0)");
}
void 
Options::saveDb(QString name, QString description){
  QSqlQuery query(QSqlDatabase::database(dbpath));
  // create material string (multi object)
  QString materialStr;
  QLOG_DEBUG() << "this->getObjectNumber()" << this->getObjectNumber();
  QLOG_DEBUG() << "this->getObjectNumber()" << objectnumber;
  for ( int i = 0 ; i < this->getObjectNumber(); i++ )
    materialStr = materialStr + this->getMaterial()[i] + " ";

  QLOG_DEBUG() << "this->getWaveMultiNumber()" << this->getWaveMultiNumber();

  this->setDescription(description);
  QString queryStr = "insert or replace into options_tbl values ('"+ 
	name + "', '"+ 
	description + "', datetime('now'), " +
        QString::number(this->getWavelength()) + ", " +
        QString::number(this->getP0()) + ", " +
        QString::number(this->getW0()) + ", '" +
        this->getBeam() + "', '" +
        QString::number(this->getWaveMultiNumber()) +  "', '" +
        this->getObject() +  "', '" +
        QString::number(this->getObjectNumber()) +  "', '" +
        this->getAnisotropy() +  "', '" +
        materialStr +  "', " +
        QString::number(this->getDiscretization()) + ", " +
        QString::number(this->getTolerance()) + ")";
  query.exec(queryStr);
  QLOG_DEBUG () << "OPTIONS_TBL query " << queryStr;
  QLOG_DEBUG () << "OPTIONS_TBL " << query.lastError().text();

  queryStr = "insert or replace into beam_tbl values ('"+ 
	name + "', '"+ 
        this->getBeam() + "', '" +
        this->getBeamFile() + "', " +
        QString::number(this->getIncidenceangle_theta_z()) +  ", " +
        QString::number(this->getIncidenceangle_phi_x()) +  ", " +
        QString::number(this->getPolarizationTM()) +  ", " +
	QString::number(this->getPolarizationTE()) +  ", " +
        QString::number(this->getPolarizationRL()) +  ", " +
        QString::number(this->getXgaus()) +  ", " +
	QString::number(this->getYgaus()) +  ", " +
        QString::number(this->getZgaus()) + ")";
  query.exec(queryStr);
  QLOG_DEBUG () << "BEAM_TBL query " << queryStr;
  QLOG_DEBUG () << "BEAM_TBL " << query.lastError().text();

  for (int i = 0 ; i < this->getObjectNumber(); i++) {
   QString queryStr = "insert or replace into object_tbl values ('"+ 
	name + "', '" +
        QString::number(i) + "', '" +
        this->getObject() + "', '" +
        this->getObjectFile() + "', " +
        QString::number(this->getSpherecoherencelength()) +  ", " +
        QString::number(this->getSpherestandardev()) +  ", " +
        QString::number(this->getSphereseed()) +  ", " +
        QString::number(this->getSphereradius().at(i)) +  ", " +
        QString::number(this->getDensity()) +  ", " +
        QString::number(this->getCubeside()) +  ", " +
        QString::number(this->getCubesidex()) +  ", " +
        QString::number(this->getCubesidey()) +  ", " +
        QString::number(this->getCubesidez()) +  ", " +
        QString::number(this->getPositionx().at(i)) +  ", " +
	QString::number(this->getPositiony().at(i)) +  ", " +
        QString::number(this->getPositionz().at(i)) +  ", " +
        QString::number(this->getDemiaxea()) +  ", " +
	QString::number(this->getDemiaxeb()) +  ", " +
        QString::number(this->getDemiaxec()) +  ", " +
        QString::number(this->getThetaobj()) +  ", " +
	QString::number(this->getPhiobj()) +  ", " +
        QString::number(this->getPsiobj()) +  ", " +
	QString::number(this->getHauteur()) + ")";
   query.exec(queryStr);
   QLOG_DEBUG () << "OBJECT_TBL query " << queryStr;
   QLOG_DEBUG () << "OBJECT_TBL " << query.lastError().text();
  }
  for (int i = 0 ; i < this->getObjectNumber(); i++) {
   QString queryStr = "insert or replace into epsilon_tbl values ('"+ 
	name + "', '" +
        QString::number(i) + "', '" +
        this->getMaterial()[i] + "', '" +
        this->getAnisotropy() + "', " +
        QString::number(real(this->getEpsilon().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon().at(i))) +  ", " +
        QString::number(real(this->getEpsilon11().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon11().at(i))) +  ", " +
        QString::number(real(this->getEpsilon12().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon12().at(i))) +  ", " +
	QString::number(real(this->getEpsilon13().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon13().at(i))) +  ", " +
        QString::number(real(this->getEpsilon21().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon21().at(i))) +  ", " +
        QString::number(real(this->getEpsilon22().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon22().at(i))) +  ", " +
        QString::number(real(this->getEpsilon23().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon23().at(i))) +  ", " +
	QString::number(real(this->getEpsilon31().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon31().at(i))) +  ", " +
        QString::number(real(this->getEpsilon32().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon32().at(i))) +  ", " +
        QString::number(real(this->getEpsilon33().at(i))) +  ", " +
	QString::number(imag(this->getEpsilon33().at(i))) + ")";
   query.exec(queryStr);
   QLOG_DEBUG () << "EPSILON_TBL query " << queryStr;
   QLOG_DEBUG () << "EPSILON_TBL " << query.lastError().text();
  }
for (int i = 0 ; i < this->getWaveMultiNumber(); i++) {
   QString queryStr = "insert or replace into wavemulti_tbl values ('"+ 
	name + "', " +
        QString::number(i) + ", " +
        QString::number(this->getThetam().at(i)) +  ", " +
	QString::number(this->getPhim().at(i)) +  ", " +
        QString::number(this->getPpm().at(i)) +  ", " +
        QString::number(this->getSsm().at(i)) +  ", " +
        QString::number(real(this->getE0m().at(i))) +  ", " +
	QString::number(imag(this->getE0m().at(i))) + ")";
   query.exec(queryStr);
   QLOG_DEBUG () << "WAVEMULTI_TBL query " << queryStr;
   QLOG_DEBUG () << "WAVEMULTI_TBL " << query.lastError().text();
  }
  queryStr = "insert or replace into run_tbl values ('"+ 
	name + "', '"+ 
        this->getMethodeit() +  "', '" +
        this->getPolarizability() +  "', " +
        QString::number(this->getQuad()) +  ", " +
        QString::number(this->getNrig()) +  ", " +
        QString::number(this->getNread()) +  ", '" +
        this->getFilereread() +  "', " +
        QString::number(this->getNmatlab()) +  ", '" +
        this->getH5File() +  "', " +
        QString::number(this->getDipolepsilon()) +  ", " +
        QString::number(this->getFarfield()) +  ", " +
        QString::number(this->getNearfield()) +  ", " +
        QString::number(this->getForce()) +  ", " +
        QString::number(this->getLocalfield()) +  ", " +
        QString::number(this->getMacroscopicfield()) +  ", " +
        QString::number(this->getNenergie()) +  ", " +
        QString::number(this->getCrosssection()) +  ", " +
	QString::number(this->getCrosssectionpoynting()) +  ", " +
        QString::number(this->getQuickdiffract()) +  ", " +
        QString::number(this->getMicroscopy()) +  ", " +
        QString::number(this->getMicroscopyFFT()) +  ", " +
        QString::number(this->getOpticalforce()) +  ", " +
        QString::number(this->getOpticalforcedensity()) +  ", " +
        QString::number(this->getOpticaltorque()) +  ", " +
	QString::number(this->getOpticaltorquedensity()) +  ", " +
        QString::number(this->getNproche()) +  ", " +
        QString::number(this->getNxx()) +  ", " +
        QString::number(this->getNyy()) +  ", " +
        QString::number(this->getNzz()) +  ", " +
        QString::number(this->getNxmp()) +  ", " +
        QString::number(this->getNymp()) +  ", " +
        QString::number(this->getNzmp()) +  ", " +    
        QString::number(this->getNtheta()) +  ", " +
        QString::number(this->getNphi()) +  ", " +
        QString::number(this->getNA()) +  ", " +
        QString::number(this->getNAinc()) +  ", " + 
        QString::number(this->getGross()) +  ", " +
        QString::number(this->getZlens()) +  ", " +
        QString::number(this->getNtypemic()) +  ", " +    
        QString::number(this->getNside()) +  ", " +    
        QString::number(this->getMeshsize()) +  ", " +
        QString::number(this->getnfft2d()) +  ", " +
        QString::number(this->getAdvancedinterface()) + ")";
  query.exec(queryStr);
  QLOG_DEBUG () << "RUN_TBL query " << queryStr;
  QLOG_DEBUG () << "RUN_TBL" << query.lastError().text();
}
void 
Options::loadDb(QString name){
  
  if ( name == "" ) name = "new";
  QSqlQuery query(QSqlDatabase::database(dbpath));
  query.exec("select name,description,wavelength,P0,W0,beam,wavemultinumber,object,objectnumber,"
             "anisotropy,discretization,tolerance "
	     "from options_tbl where name='" + name + "'");
  QLOG_DEBUG () << query.lastError().text();
    
  while (query.next()) {
    QLOG_DEBUG () << "SELECT OPTIONS_TBL " << query.value(7).toString();
    QLOG_DEBUG () << "SELECT OPTIONS_TBL beam " << query.value(5).toString();
    this->setName(query.value(0).toString());
    this->setDescription(query.value(1).toString());
    this->setWavelength(query.value(2).toDouble());
    this->setP0(query.value(3).toDouble());
    this->setW0(query.value(4).toDouble());
    this->setBeam(query.value(5).toString());
    this->setWaveMultiNumber(query.value(6).toInt());
    this->setObject(query.value(7).toString());
    this->setObjectNumber(query.value(8).toInt());
    this->setAnisotropy(query.value(9).toString());
    this->setDiscretization(query.value(10).toInt());
    this->setTolerance(query.value(11).toDouble());
  }
  query.exec("select beamfile,incidenceangle_theta_z,incidenceangle_phi_x,"
	     "polarizationTM,polarizationTE,polarizationRL,xgaus,ygaus,zgaus "
	     "from beam_tbl where parent='" + name + "'");
  QLOG_DEBUG () << "SELECT OPTIONS_TBL " << query.lastError().text();
  while (query.next()) {
    QLOG_DEBUG () << "SELECT OPTIONS_TBL " << query.value(4).toDouble();
    this->setBeamFile(query.value(0).toString());
    this->setIncidenceangle_theta_z(query.value(1).toDouble());
    this->setIncidenceangle_phi_x(query.value(2).toDouble());
    this->setPolarizationTM(query.value(3).toDouble());
    this->setPolarizationTE(query.value(4).toDouble());
    this->setPolarizationRL(query.value(5).toDouble());
    this->setXgaus(query.value(6).toDouble());
    this->setYgaus(query.value(7).toDouble());
    this->setZgaus(query.value(8).toDouble());
  }
  QLOG_DEBUG() << "Options::loadDb>> Object Number:" <<  this->getObjectNumber();
 
  for (int i = 0 ; i < this->getObjectNumber(); i++) {
   QString queryStr = "select objectfile,spherecoherencelength,spherestandardev,sphereseed,sphereradius,"
                  "density,cubeside,cubesidex,cubesidey,cubesidez,"
                  "positionx,positiony,positionz,demiaxea,demiaxeb,demiaxec,"
                  "thetaobj,phiobj,psiobj,hauteur "
                  "from object_tbl where parent='" + name + "' "
                  "and objectnum=" + QString::number(i);
   QLOG_DEBUG() << "Options::loadDb>> query :" <<  queryStr ;
   query.exec(queryStr);
   QLOG_DEBUG () << "SELECT OBJECT_TBL " << this->getObjectNumber() << " " << query.lastError().text();
   while (query.next()) {
     QLOG_DEBUG() << "SPHERE RADIUS " <<  query.value(4).toDouble();
     this->setObjectFile(query.value(0).toString());
     this->setSpherecoherencelength(query.value(1).toDouble());
     this->setSpherestandardev(query.value(2).toDouble());
     this->setSphereseed(query.value(3).toInt());
     sphereradius.replace(i,query.value(4).toDouble());
     this->setDensity(query.value(5).toDouble());
     this->setCubeside(query.value(6).toDouble());
     this->setCubesidex(query.value(7).toDouble());
     this->setCubesidey(query.value(8).toDouble());
     this->setCubesidez(query.value(9).toDouble());
     positionx.replace(i,query.value(10).toDouble());
     positiony.replace(i,query.value(11).toDouble());
     positionz.replace(i,query.value(12).toDouble());
     this->setDemiaxea(query.value(13).toDouble());
     this->setDemiaxeb(query.value(14).toDouble());
     this->setDemiaxec(query.value(15).toDouble());
     this->setThetaobj(query.value(16).toDouble());
     this->setPhiobj(query.value(17).toDouble());
     this->setPsiobj(query.value(18).toDouble());
     this->setHauteur(query.value(19).toDouble());
   }
  }
  for (int i = 0 ; i < this->getWaveMultiNumber(); i++) {
    QString queryStr = "select thetam,phim,ppm,ssm,E0mr,E0mi "
	       "from wavemulti_tbl where parent='" + name + "' "
               "and wavemultinum=" + QString::number(i);
    QLOG_DEBUG() << "Options::loadDb>> query :" <<  queryStr ;
    query.exec(queryStr);
    QLOG_DEBUG () << "SELECT WAVEMULTI_TBL" << query.lastError().text();
    while (query.next()) {
      QLOG_DEBUG () << "SELECT WAVEMULTI_TBL" << query.value(0) << " " << query.value(1);
      QLOG_DEBUG () << "SELECT WAVEMULTI_TBL" << query.value(2) << " " << query.value(3);
      thetam.replace(i,query.value(0).toDouble());
      phim.replace(i,query.value(1).toDouble());
      ppm.replace(i,query.value(2).toDouble());
      ssm.replace(i,query.value(3).toDouble());
      E0m.replace(i,dcmplx(query.value(4).toDouble(),query.value(5).toDouble()));
    }
  }
  for (int i = 0 ; i < this->getObjectNumber(); i++) {
    QString queryStr = "select epsilonr,epsiloni,epsilon11r,epsilon11i,epsilon12r,epsilon12i,"
	       "epsilon13r,epsilon13i,epsilon21r,epsilon21i,epsilon22r,epsilon22i,"
	       "epsilon23r,epsilon23i,epsilon31r,epsilon31i,epsilon32r,epsilon32i,"
	       "epsilon33r,epsilon33i,material "
	       "from epsilon_tbl where parent='" + name + "' "
               "and objectnum=" + QString::number(i);
    QLOG_DEBUG() << "Options::loadDb>> query :" <<  queryStr ;
    query.exec(queryStr);
    QLOG_DEBUG () << "SELECT EPSILON_TBL" << query.lastError().text();
    while (query.next()) {
      QLOG_DEBUG () << "SELECT EPSILON_TBL" << query.value(0) << " " << query.value(1);
      QLOG_DEBUG () << "SELECT EPSILON_TBL" << query.value(2) << " " << query.value(3);
      epsilon.replace(i,dcmplx(query.value(0).toDouble(),query.value(1).toDouble()));
      epsilon11.replace(i,dcmplx(query.value(2).toDouble(),query.value(3).toDouble()));
      epsilon12.replace(i,dcmplx(query.value(4).toDouble(),query.value(5).toDouble()));
      epsilon13.replace(i,dcmplx(query.value(6).toDouble(),query.value(7).toDouble()));
      epsilon21.replace(i,dcmplx(query.value(8).toDouble(),query.value(9).toDouble()));
      epsilon22.replace(i,dcmplx(query.value(10).toDouble(),query.value(11).toDouble()));
      epsilon23.replace(i,dcmplx(query.value(12).toDouble(),query.value(13).toDouble()));
      epsilon31.replace(i,dcmplx(query.value(14).toDouble(),query.value(15).toDouble()));
      epsilon32.replace(i,dcmplx(query.value(16).toDouble(),query.value(17).toDouble()));
      epsilon33.replace(i,dcmplx(query.value(18).toDouble(),query.value(19).toDouble()));
      material.replace(i,query.value(20).toString());
    }
  }
  query.exec("select methode,polarizability,quad,nrig,nread,filereread,nmatlab,fileh5,dipolepsilon,farfield,nearfield,force,"
	     "localfield,macroscopicfield,nenergie,crosssection,crosssectionpoynting,quickdiffract,microscopy,"
	     "microscopyFFT,opticalforce,opticalforcedensity,opticaltorque,opticaltorquedensity,nproche," 
	     "nxx,nyy,nzz,nxmp,nymp,nzmp,ntheta,nphi,na,nainc,gross,zlens,ntypemic,nside,meshsize,nfft2d,advancedinterface from run_tbl where parent='" + name + "'");
  QLOG_DEBUG () << "SELECT RUN_TBL" << query.lastError().text();
  QLOG_DEBUG () << "SELECT RUN_TBL" << query.value(37).toInt();
  while (query.next()) {
    this->setMethodeit(query.value(0).toString());
    this->setPolarizability(query.value(1).toString());
    this->setQuad(query.value(2).toInt());
    this->setNrig(query.value(3).toInt());
    this->setNread(query.value(4).toInt());
    this->setFilereread(query.value(5).toString());
    this->setNmatlab(query.value(6).toInt());
    this->setH5File(query.value(7).toString());
    this->setDipolepsilon(query.value(8).toInt());
    this->setFarfield(query.value(9).toInt());
    this->setNearfield(query.value(10).toInt());
    this->setForce(query.value(11).toInt());
    this->setLocalfield(query.value(12).toInt());
    this->setMacroscopicfield(query.value(13).toInt());
    this->setNenergie(query.value(14).toInt());
    this->setCrosssection(query.value(15).toInt());
    this->setCrosssectionpoynting(query.value(16).toInt());
    this->setQuickdiffract(query.value(17).toInt());
    this->setMicroscopy(query.value(18).toInt());
    this->setMicroscopyFFT(query.value(19).toInt());
    this->setOpticalforce(query.value(20).toInt());
    this->setOpticalforcedensity(query.value(21).toInt());
    this->setOpticaltorque(query.value(22).toInt());
    this->setOpticaltorquedensity(query.value(23).toInt());
    this->setNproche(query.value(24).toInt());
    this->setNxx(query.value(25).toInt());
    this->setNyy(query.value(26).toInt());
    this->setNzz(query.value(27).toInt());
    this->setNxmp(query.value(28).toInt());
    this->setNymp(query.value(29).toInt());
    this->setNzmp(query.value(30).toInt());
    this->setNtheta(query.value(31).toInt());
    this->setNphi(query.value(32).toInt());
    this->setNA(query.value(33).toDouble());
    this->setNAinc(query.value(34).toDouble());    
    this->setGross(query.value(35).toDouble());
    this->setZlens(query.value(36).toDouble());    
    this->setNtypemic(query.value(37).toInt());
    this->setNside(query.value(38).toInt());
    this->setMeshsize(query.value(39).toDouble());
    this->setnfft2d(query.value(40).toInt());
    this->setAdvancedinterface(query.value(41).toInt());
    QLOG_DEBUG () << "SELECT TYPEMIC" << query.value(37).toInt();
  }
}
void 
Options::removeDb(QString name){
  QSqlQuery query(QSqlDatabase::database(dbpath));
  query.exec("delete from options_tbl where name='"+ name + "'");
  QLOG_DEBUG () << "DELETE FROM OPTIONS_TBL" << query.lastError().text();
  query.exec("delete from beam_tbl where parent='"+ name + "'");
  QLOG_DEBUG () << "DELETE FROM BEAM_TBL" << query.lastError().text();
  query.exec("delete from object_tbl where parent='"+ name + "'");
  QLOG_DEBUG () << "DELETE FROM OBJECT_TBL" << query.lastError().text();
  query.exec("delete from epsilon_tbl where parent='"+ name + "'");
  QLOG_DEBUG () << "DELETE FROM EPSILON_TBL" << query.lastError().text();
  query.exec("delete from run_tbl where parent='"+ name + "'");
  QLOG_DEBUG () << "DELETE FROM RUN_TBL" << query.lastError().text();

}
void 
Options::setName(QString _runname){
  runname = _runname;
}
void 
Options::setDescription(QString _description){
  description = _description;
}
void 
Options::setWavelength(double _wavelength){
  wavelength = _wavelength;
}
void 
Options::setBeam(QString _beam){
  beam = _beam;
}
void 
Options::setBeamFile(QString _beamfile){
  beamfile = _beamfile;
}
void 
Options::setObject(QString _object){
  object = _object;
}
void 
Options::setObjectFile(QString _objectfile){
  objectfile = _objectfile;
}
void 
Options::setAnisotropy(QString _anisotropy){
  anisotropy = _anisotropy;
}
void 
Options::setMaterial(QVector<QString> _material){
  for (int i = 0 ; i < objectnumber; i++)
    material.replace(i,_material.at(i));
}
void 
Options::setFilereread(QString _filereread){
  filereread = _filereread;
}
void 
Options::setDiscretization(int _discretization){
  discretization = _discretization;
}
void 
Options::setTolerance(double _tolerance){
  tolerance = _tolerance;
}
void 
Options::setMethodeit(QString _methodeit){
  methodeit = _methodeit;
}
void 
Options::setH5File(QString _fileh5){
  fileh5 = _fileh5;
}

void 
Options::setPolarizability(QString _polarizability){
  polarizability = _polarizability;
}
void 
Options::setQuad(int _quad){
  quad = _quad;
}
// ellipsoid
void 
Options::setDemiaxea(double _demiaxea) {
 demiaxea = _demiaxea;
}
void 
Options::setDemiaxeb(double _demiaxeb) {
 demiaxeb = _demiaxeb;
}
void 
Options::setDemiaxec(double _demiaxec) {
 demiaxec = _demiaxec;
}
void 
Options::setThetaobj(double _thetaobj) {
 thetaobj = _thetaobj;
}
void 
Options::setPhiobj(double _phiobj) {
 phiobj = _phiobj;
}
void 
Options::setPsiobj(double _psiobj) {
 psiobj = _psiobj;
}
// planewave
void 
Options::setIncidenceangle_theta_z(double _incidenceangle_theta_z) {
 incidenceangle_theta_z = _incidenceangle_theta_z;
}
void 
Options::setIncidenceangle_phi_x(double _incidenceangle_phi_x) {
 incidenceangle_phi_x = _incidenceangle_phi_x;
}
void 
Options::setPolarizationTM(double _polarizationTM) {
 polarizationTM = _polarizationTM;
}
void 
Options::setPolarizationTE(double _polarizationTE) {
 polarizationTE = _polarizationTE;
}
void 
Options::setPolarizationRL(double _polarizationRL) {
 polarizationRL = _polarizationRL;
}
void 
Options::setXgaus(double _xgaus) {
 xgaus = _xgaus;
}
void 
Options::setYgaus(double _ygaus) {
 ygaus = _ygaus;
}
void 
Options::setZgaus(double _zgaus) {
 zgaus = _zgaus;
}
void 
Options::setP0(double _P0) {
 P0 = _P0;
}
void 
Options::setW0(double _W0) {
 W0 = _W0;
}
// wave multi
void 
Options::setThetam(QVector<double> _thetam) {
 for (int i = 0 ; i < wavemultinumber; i++)
   thetam.replace(i,_thetam.at(i));
}
void 
Options::setPhim(QVector<double> _phim) {
 for (int i = 0 ; i < wavemultinumber; i++)
   phim.replace(i,_phim.at(i));
}
void 
Options::setPpm(QVector<double> _ppm) {
 for (int i = 0 ; i < wavemultinumber; i++)
   ppm.replace(i,_ppm.at(i));
}
void 
Options::setSsm(QVector<double> _ssm) {
 for (int i = 0 ; i < wavemultinumber; i++)
   ssm.replace(i,_ssm.at(i));
}
void 
Options::setE0m(QVector<dcmplx> _E0m) {
  for (int i = 0 ; i < wavemultinumber; i++)
   E0m.replace(i,_E0m.at(i));
}
// sphere (including multiple)
void 
Options::setDensity(double _density) {
  density = _density;
}
void 
Options::setCubeside(double _cubeside) {
  cubeside = _cubeside;
}
void 
Options::setCubesidex(double _cubesidex) {
  cubesidex = _cubesidex;
}
void 
Options::setCubesidey(double _cubesidey) {
  cubesidey = _cubesidey;
}
void 
Options::setCubesidez(double _cubesidez) {
  cubesidez = _cubesidez;
}
void 
Options::setHauteur(double _hauteur) {
  hauteur = _hauteur;
}
void 
Options::setObjectNumber(int _objectnumber) {
  objectnumber = _objectnumber;
}
void 
Options::setWaveMultiNumber(int _wavemultinumber) {
  wavemultinumber = _wavemultinumber;
}
void
Options::setSpherecoherencelength(double _lc) {
   lc = _lc;
}
void
Options::setSpherestandardev(double _hc) {
   hc = _hc;
}
void
Options::setSphereseed(int _sphereseed) {
   sphereseed = _sphereseed;
}
void 
Options::setSphereradius(QVector<double> _sphereradius) {
  for (int i = 0 ; i < objectnumber; i++)
   sphereradius.replace(i,_sphereradius.at(i));
}
void 
Options::setPositionx(QVector<double> _positionx) {
  for (int i = 0 ; i < objectnumber; i++)
   positionx.replace(i,_positionx.at(i));
}
void 
Options::setPositiony(QVector<double> _positiony) {
  for (int i = 0 ; i < objectnumber; i++)
   positiony.replace(i,_positiony.at(i));
}
void 
Options::setPositionz(QVector<double> _positionz) {
  for (int i = 0 ; i < objectnumber; i++)
   positionz.replace(i,_positionz.at(i));
}
void 
Options::setEpsilon(QVector<dcmplx> _epsilon) {
  for (int i = 0 ; i < objectnumber; i++)
    epsilon.replace(i,_epsilon.at(i));
}
void 
Options::setEpsilon11(QVector<dcmplx> _epsilon11) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon11.replace(i,_epsilon11.at(i));
}
void 
Options::setEpsilon12(QVector<dcmplx> _epsilon12) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon12.replace(i,_epsilon12.at(i));
}
void 
Options::setEpsilon13(QVector<dcmplx> _epsilon13) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon13.replace(i,_epsilon13.at(i));
}
void 
Options::setEpsilon21(QVector<dcmplx> _epsilon21) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon21.replace(i,_epsilon21.at(i));
}
void 
Options::setEpsilon22(QVector<dcmplx> _epsilon22) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon22.replace(i,_epsilon22.at(i));
}
void 
Options::setEpsilon23(QVector<dcmplx> _epsilon23) {
  for (int i = 0 ; i < objectnumber; i++)
   epsilon23.replace(i,_epsilon23.at(i));
}
void Options::setEpsilon31(QVector<dcmplx> _epsilon31){
  for (int i = 0 ; i < objectnumber; i++)
   epsilon31.replace(i,_epsilon31.at(i));
}
void 
Options::setEpsilon32(QVector<dcmplx> _epsilon32){
  for (int i = 0 ; i < objectnumber; i++)
   epsilon32.replace(i,_epsilon32.at(i));
}
void 
Options::setEpsilon33(QVector<dcmplx> _epsilon33){
  for (int i = 0 ; i < objectnumber; i++)
   epsilon33.replace(i,_epsilon33.at(i));
}
void 
Options::setNread(int _nread){
  nread = _nread;
}
void
Options::setNmatlab(int _nmatlab){
  nmatlab  = _nmatlab;
}
void
Options::setAdvancedinterface(int _advancedinterface){
  advancedinterface  = _advancedinterface;
}
void 
Options::setDipolepsilon(int _dipolepsilon){
  dipolepsilon = _dipolepsilon;
}
void 
Options::setFarfield(int _farfield){
  farfield = _farfield;
}
void 
Options::setNearfield(int _nearfield){
  nearfield = _nearfield;
}
void 
Options::setForce(int _force){
  force = _force;
}
void 
Options::setLocalfield(int _localfield){
  localfield = _localfield;
}
void 
Options::setMacroscopicfield(int _macroscopicfield){
  macroscopicfield = _macroscopicfield;
}
void 
Options::setNenergie(int _nenergie){
  nenergie = _nenergie;
}
void 
Options::setCrosssection(int _crosssection){
  crosssection = _crosssection;
}
void 
Options::setCrosssectionpoynting(int _crosssectionpoynting){
  crosssectionpoynting = _crosssectionpoynting;
}
void 
Options::setQuickdiffract(int _quickdiffract){
  quickdiffract = _quickdiffract;
}
void 
Options::setNrig(int _nrig){
  nrig = _nrig;
}
void 
Options::setMicroscopy(int _microscopy){
  microscopy = _microscopy;
}
void 
Options::setMicroscopyFFT(int _microscopyFFT){
  microscopyFFT = _microscopyFFT;
}
void 
Options::setOpticalforce(int _opticalforce){
  opticalforce = _opticalforce;
}
void 
Options::setOpticalforcedensity(int _opticalforcedensity){
  opticalforcedensity = _opticalforcedensity;
}
void 
Options::setOpticaltorque(int _opticaltorque){
  opticaltorque = _opticaltorque;
}
void 
Options::setOpticaltorquedensity(int _opticaltorquedensity){
  opticaltorquedensity = _opticaltorquedensity;
}
void 
Options::setNproche(int _nproche){
  nproche = _nproche;
}
void 
Options::setNxx( int _nxx){
  nxx = _nxx;
}
void 
Options::setNyy( int _nyy){
  nyy = _nyy;
}
void 
Options::setNzz( int _nzz){
  nzz = _nzz;
}
void 
Options::setNxm( int _nxm){
  nxm = _nxm;
}
void 
Options::setNym( int _nym){
  nym = _nym;
}
void 
Options::setNzm( int _nzm){
  nzm = _nzm;
}
void
Options::setNxmp( int _nxmp){
  nxmp = _nxmp;
}
void 
Options::setNymp( int _nymp){
  nymp = _nymp;
}
void 
Options::setNzmp( int _nzmp){
  nzmp = _nzmp;
}
void 
Options::setNtheta( int _ntheta){
  ntheta = _ntheta;
}
void 
Options::setNphi( int _nphi){
  nphi = _nphi;
}
void 
Options::setNA( double _na){
  na = _na;
}
void 
Options::setNAinc( double _nainc){
  nainc = _nainc;
}
void 
Options::setGross( double _gross){
  gross = _gross;
}
void 
Options::setZlens( double _zlens){
  zlens = _zlens;
}
void 
Options::setNtypemic(int _ntypemic){
  ntypemic = _ntypemic;
}
void 
Options::setNside(int _nside){
  nside = _nside;
}

void 
Options::setMeshsize( double _meshsize){
  meshsize = _meshsize;
}
void 
Options::setnfft2d( int _nfft2d){
  nfft2d = _nfft2d;
}
QString
Options::getName(){
  return runname;
}
QString
Options::getDescription(){
  return description;
}
QString 
Options::getFilereread(){
  return filereread;
}
double 
Options::getWavelength(){
  return wavelength;
}
QString 
Options::getBeam(){
  return beam;
}
QString 
Options::getBeamFile(){
  return beamfile;
}
QString 
Options::getObject(){
  return object;
}
QString 
Options::getObjectFile(){
  return objectfile;
}
QString 
Options::getAnisotropy(){
  return anisotropy;
}
QVector<QString> 
Options::getMaterial(){
  return material;
}
QString 
Options::getMethodeit(){
  return methodeit;
}
QString 
Options::getH5File(){
  return fileh5;
}
QString 
Options::getPolarizability(){
  return polarizability;
}
int 
Options::getQuad(){
  return quad;
}
int 
Options::getDiscretization(){
  return discretization;
}
double 
Options::getTolerance(){
  return tolerance;
}
// ellipsoid
double 
Options::getDemiaxea(){
  return demiaxea;
}
double 
Options::getDemiaxeb(){
  return demiaxeb;
}
double 
Options::getDemiaxec(){
  return demiaxec;
}
double 
Options::getThetaobj(){
  return thetaobj;
}
double 
Options::getPhiobj(){
  return phiobj;
} 
double 
Options::getPsiobj(){
  return psiobj;
} 
// planewave
double 
Options::getIncidenceangle_theta_z(){
  return incidenceangle_theta_z;
} 
double 
Options::getIncidenceangle_phi_x(){
  return incidenceangle_phi_x;
} 
double 
Options::getPolarizationTM(){
  return polarizationTM;
} 
double
Options::getPolarizationTE(){
  return polarizationTE;
} 
double 
Options::getPolarizationRL(){
  return polarizationRL;
}
double 
Options::getXgaus(){
  return xgaus;
}
double 
Options::getYgaus(){
  return ygaus;
}
double 
Options::getZgaus(){
  return zgaus;
}
double 
Options::getP0(){
  return P0;
} 
double 
Options::getW0(){
  return W0;
}
// multi wave
QVector<double>
Options::getThetam(){
  return thetam;
}
QVector<double>
Options::getPhim(){
  return phim;
}
QVector<double> 
Options::getPpm(){
  return ppm;
}
QVector<double> 
Options::getSsm(){
  return ssm;
}
QVector<dcmplx>
Options::getE0m(){
  return E0m;
}

// sphere (includes multiple)
double
Options::getDensity(){
  return density;
} 
double
Options::getCubeside(){
  return cubeside;
} 
double
Options::getCubesidex(){
  return cubesidex;
} 
double
Options::getCubesidey(){
  return cubesidey;
} 
double
Options::getCubesidez(){
  return cubesidez;
} 
double 
Options::getHauteur(){
  return hauteur;
} 
int 
Options::getObjectNumber(){
  return objectnumber;
}
int 
Options::getWaveMultiNumber(){
  return wavemultinumber;
}
double
Options::getSpherecoherencelength(){
  return lc;
}
double
Options::getSpherestandardev(){
  return hc;
}
int
Options::getSphereseed(){
  return sphereseed;
}
QVector<double>
Options::getSphereradius(){
  return sphereradius;
}
QVector<double> 
Options::getPositionx(){
  return positionx;
} 
QVector<double> 
Options::getPositiony(){
  return positiony;
} 
QVector<double>
Options::getPositionz(){
  return positionz;
}
QVector<dcmplx>
Options::getEpsilon() {
  return epsilon;
}
QVector<dcmplx>
Options::getEpsilon11(){
  return epsilon11;
}
QVector<dcmplx>
Options::getEpsilon12(){
  return epsilon12;
}
QVector<dcmplx>
Options::getEpsilon13(){
  return epsilon13;
}
QVector<dcmplx>
Options::getEpsilon21(){
  return epsilon21;
}
QVector<dcmplx> 
Options::getEpsilon22(){
  return epsilon22;
}
QVector<dcmplx> 
Options::getEpsilon23(){
  return epsilon23;
}
QVector<dcmplx>
Options::getEpsilon31(){
  return epsilon31;
}
QVector<dcmplx> 
Options::getEpsilon32(){
  return epsilon32;
}
QVector<dcmplx> 
Options::getEpsilon33(){
  return epsilon33;
}
// Options check box parameters
int    
Options::getNread(){
  return nread;
}
int    
Options::getNmatlab(){
  return nmatlab;
}
int    
Options::getAdvancedinterface(){
  return advancedinterface;
}
int   
Options::getDipolepsilon(){
  return dipolepsilon;
}
int   
Options::getFarfield(){
  return farfield;
}
int   
Options::getNearfield(){
  return nearfield;
}
int   
Options::getForce(){
  return force;
}
int   
Options::getLocalfield(){
  return localfield;
}
int   
Options::getMacroscopicfield(){
  return macroscopicfield;
}
int    
Options::getNenergie(){
  return nenergie;
}
int    
Options::getCrosssection(){
  return crosssection;
}
int    
Options::getCrosssectionpoynting(){
  return crosssectionpoynting;
}
int    
Options::getQuickdiffract(){
  return quickdiffract;
}
int   
Options::getNrig(){
  return nrig;
}
int    
Options::getMicroscopy(){
  return microscopy;
}
int    
Options::getMicroscopyFFT(){
  return microscopyFFT;
}
int   
Options::getOpticalforce(){
  return opticalforce;
}
int    
Options::getOpticalforcedensity(){
  return opticalforcedensity;
}
int    
Options::getOpticaltorque(){
  return opticaltorque;
}
int    
Options::getOpticaltorquedensity(){
  return opticaltorquedensity;
}
// Options parameters
int    
Options::getNproche(){
  return nproche;
}
int    
Options::getNxx(){
  return nxx;
}
int    
Options::getNyy(){
  return nyy;
}
int    
Options::getNzz(){
  return nzz;
}
int    
Options::getNxm(){
  return nxm;
}
int    
Options::getNym(){
  return nym;
}
int    
Options::getNzm(){
  return nzm;
}
int    
Options::getNxmp(){
  return nxmp;
}
int    
Options::getNymp(){
  return nymp;
}
int    
Options::getNzmp(){
  return nzmp;
}
int    
Options::getNtheta(){
  return ntheta;
}
int    
Options::getNphi(){
  return nphi;
}
int  
Options::getnfft2d(){
  return nfft2d;
}
double  
Options::getNA(){
  return na;
}
double  
Options::getNAinc(){
  return nainc;
}
double  
Options::getGross(){
  return gross;
}
double  
Options::getZlens(){
  return zlens;
}
int  
Options::getNside(){
  return nside;
}
int  
Options::getNtypemic(){
  return ntypemic;
}
double 
Options::getMeshsize(){
  return meshsize;
} 
QVector<QColor>* 
Options::getColors(){
  return colors;
}

