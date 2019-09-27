#ifndef OPTIONS_H
#define OPTIONS_H

#include <QtSql>
#include <QtGui>
#include "QsLog.h"

#include <complex>

using namespace std;

typedef complex<double> dcmplx;

#define MAX_OBJECT_NUMBER 20

#define MAX_WAVEMULTI_NUMBER 10

class Options
{
   public:
   Options();
   ~Options();

   void initDb();
   void saveDb(QString name, QString description);
   void loadDb(QString name);
   void removeDb(QString name);
   void initOptions();

   QString  getName();
   QString  getDescription();
   QString  getFilereread();
   QString  getH5File();
   double   getWavelength();
   QString  getBeam();
   QString  getBeamFile();
   QString  getObject();
   QString  getObjectFile();
   QString  getAnisotropy();
   QVector<QString>  getMaterial();
   QString  getMethodeit();
   QString  getPolarizability();
   int      getQuad();
   QString  getMethode();
   int      getDiscretization();
   double   getTolerance();   

   // ellipsoid
   double getDemiaxea();
   double getDemiaxeb();
   double getDemiaxec();
   double getThetaobj();
   double getPhiobj();
   double getPsiobj();

   // beam shape
   double getIncidenceangle_theta_z();
   double getIncidenceangle_phi_x();
   double getPolarizationTM();
   double getPolarizationTE();
   double getPolarizationRL();
   double getXgaus();
   double getYgaus();
   double getZgaus();
   double getP0();
   double getW0();

   // wave multi
   QVector<double> getThetam();
   QVector<double> getPhim();
   QVector<double> getPpm();
   QVector<double> getSsm();
   QVector<dcmplx> getE0m();

   // sphere (includes multiple)
   double   getDensity();
   double   getCubeside();
   double   getCubesidex();
   double   getCubesidey();
   double   getCubesidez();
   double   getHauteur();
   int      getObjectNumber();
   int      getWaveMultiNumber();
   int      getSphereseed();
   double   getSpherecoherencelength();
   double   getSpherestandardev();
   QVector<double> getSphereradius();
   QVector<double> getPositionx();
   QVector<double> getPositiony();
   QVector<double> getPositionz();
   QVector<dcmplx> getEpsilon();
   QVector<dcmplx> getEpsilon11();
   QVector<dcmplx> getEpsilon12();
   QVector<dcmplx> getEpsilon13();
   QVector<dcmplx> getEpsilon21();
   QVector<dcmplx> getEpsilon22();
   QVector<dcmplx> getEpsilon23();
   QVector<dcmplx> getEpsilon31();
   QVector<dcmplx> getEpsilon32();
   QVector<dcmplx> getEpsilon33();

   int  getNread();
   int  getNmatlab();
   int  getAdvancedinterface();
   int  getDipolepsilon();
   int  getFarfield();
   int  getNearfield();
   int  getForce();
   int  getLocalfield();
   int  getMacroscopicfield();
   int  getCrosssection();
   int  getCrosssectionpoynting();
   int  getNenergie();
   int  getQuickdiffract();
   int  getNrig();
   int  getMicroscopy();
   int  getMicroscopyFFT();
   int  getOpticalforce();
   int  getOpticalforcedensity();
   int  getOpticaltorque();
   int  getOpticaltorquedensity();
   int  getNproche(); 

   int  getNxx();
   int  getNyy();
   int  getNzz();   
   int  getNxm();
   int  getNym();
   int  getNzm();
   int  getNxmp();
   int  getNymp();
   int  getNzmp();
   int  getNtheta();
   int  getNphi();
   int  getnfft2d();
   double getNA();
   double getNAinc();
   double getGross();
   double getZlens();   
   int getNtypemic();
   int getNside();
   double getMeshsize();
   
   QVector<QColor>* getColors();

   void setName(QString _runname);
   void setDescription(QString _description);
   void setFilereread(QString _filereread);
   void setH5File(QString _fileh5);
   void setWavelength(double _wavelength);
   void setBeam(QString _beam);
   void setBeamFile(QString _beamfile);
   void setObject(QString _object);
   void setObjectFile(QString _objectfile);
   void setAnisotropy(QString _anisotropy);
   void setMaterial(QVector<QString> _material);
   void setMethodeit(QString _methodeit);
   void setPolarizability(QString _polarizability);
   void setQuad(int _quad);
   void setMethode(QString _methode);
   void setDiscretization(int _discretization);
   void setTolerance(double _tolerance);

   // ellipsoid
   void setDemiaxea( double _demiaxea);
   void setDemiaxeb( double _demiaxeb);
   void setDemiaxec( double _demiaxec);
   void setThetaobj( double _thetaobj);
   void setPhiobj( double _psiobj);
   void setPsiobj( double _psiobj);

   // beam shape
   void setIncidenceangle_theta_z(double _incidenceangle_theta_z);
   void setIncidenceangle_phi_x(double _incidenceangle_phi_x);
   void setPolarizationTM(double _polarizationTM);
   void setPolarizationTE(double _polarizationTE);
   void setPolarizationRL(double _polarizationRL);
   void setXgaus(double _xgaus);
   void setYgaus(double _ygaus);
   void setZgaus(double _zgaus);
   void setP0(double P0);
   void setW0(double W0);

   // wave multi
   void setThetam(QVector<double> _thetam);
   void setPhim(QVector<double> _phim);
   void setPpm(QVector<double> _ppm);
   void setSsm(QVector<double> _ssm);
   void setE0m(QVector<dcmplx> _E0m);

   // objects (includes multiple)
   void setDensity(double _density);
   void setCubeside(double _cubeside);
   void setCubesidex(double _cubesidex);
   void setCubesidey(double _cubesidey);
   void setCubesidez(double _cubesidez);
   void setHauteur(double _hauteur);
   void setObjectNumber(int _objectnumber);
   void setWaveMultiNumber(int _wavemultinumber);
   void setSphereseed(int _sphereseed);
   void setSpherecoherencelength(double _lc);
   void setSpherestandardev(double _hc);
   void setSphereradius(QVector<double> _sphereradius);
   void setPositionx(QVector<double> _positionx);
   void setPositiony(QVector<double> _positiony);
   void setPositionz(QVector<double> _positionz);
   void setEpsilon(QVector<dcmplx> _epsilon);
   void setEpsilon11(QVector<dcmplx> _epsilon11);
   void setEpsilon12(QVector<dcmplx> _epsilon12);
   void setEpsilon13(QVector<dcmplx> _epsilon13);
   void setEpsilon21(QVector<dcmplx> _epsilon21);
   void setEpsilon22(QVector<dcmplx> _epsilon22);
   void setEpsilon23(QVector<dcmplx> _epsilon23);
   void setEpsilon31(QVector<dcmplx> _epsilon31);
   void setEpsilon32(QVector<dcmplx> _epsilon32);
   void setEpsilon33(QVector<dcmplx> _epsilon33);

   void setNread(int _nread);
   void setNmatlab(int _nmatlab);
   void setAdvancedinterface(int _advancedinterface);
   void setDipolepsilon(int _dipolepsilon);
   void setFarfield(int _farfield);
   void setNearfield(int _nearfield);
   void setForce(int _force);
   void setLocalfield(int _localfield);
   void setMacroscopicfield(int _macroscopicfield);
   void setCrosssection(int _crosssection);
   void setCrosssectionpoynting(int _crosssectionpyonting);
   void setNenergie(int _nenergie);
   void setQuickdiffract(int _quickdiffract);
   void setNrig(int _nrig);
   void setMicroscopy(int _microscopy);
   void setMicroscopyFFT(int _microscopyFFT);
   void setOpticalforce(int _opticalforce);
   void setOpticalforcedensity(int _opticalforcedensity);
   void setOpticaltorque(int _opticaltorque);
   void setOpticaltorquedensity(int _opticaltorquedensity);
   void setNproche(int _nproche);

   void setNxx(int _nxx);
   void setNyy(int _nyy);
   void setNzz(int _nzz);
   void setNxm(int _nxm);
   void setNym(int _nym);
   void setNzm(int _nzm);
   void setNxmp(int _nxmp);
   void setNymp(int _nymp);
   void setNzmp(int _nzmp);
   void setNtheta(int _ntheta);
   void setNphi(int _nphi);
   void setNA(double _na);
   void setNAinc(double _nainc);
   void setGross(double _gross);
   void setZlens(double _zlens);   
   void setNtypemic(int _ntypemic);
   void setNside(int _nside);
   void setMeshsize(double _meshsize);
   void setnfft2d(int _nfft2d);
 
   QStringList beamList;
   QStringList objectList;
   QStringList anisotropyList;
   QStringList materialList;
   QStringList methodeitList;
   QStringList quadList;
   QStringList nfft2dList;
   QStringList polarizabilityList;
   QStringList nrigList;
   QStringList nmatlabList;
   QStringList ntypemicList;
   QStringList nsideList;
   QStringList rangeofstudyList;

  private:

    QString dbpath;
    QString runname;
    QString description;
    QString filereread;
    double  wavelength;
    QString beam;
    QString beamfile;
    QString object;
    QString objectfile;
    QString anisotropy;
    QString fileh5;
    QString methodeit;
    QString polarizability;
    int     quad;
    int     discretization;
    double  tolerance;
    
    // ellipsoid
    double demiaxea;
    double demiaxeb;
    double demiaxec;
    double thetaobj;
    double phiobj;
    double psiobj;

    // planewave
    double incidenceangle_theta_z;
    double incidenceangle_phi_x;
    double polarizationTM;
    double polarizationTE;
    double polarizationRL;
    double xgaus;
    double ygaus;
    double zgaus;
    double P0;
    double W0;

    // wave multi
    QVector<double> thetam;
    QVector<double> phim;
    QVector<double> ppm;
    QVector<double> ssm;
    QVector<dcmplx> E0m;

    // sphere (include multiple)
    double density;
    double cubeside;
    double cubesidex;
    double cubesidey;
    double cubesidez;
    double hauteur;
    int    objectnumber;
    int    wavemultinumber;
    int    sphereseed;
    double hc;
    double lc;
    QVector<double> sphereradius;
    QVector<double> positionx;
    QVector<double> positiony;
    QVector<double> positionz;
    QVector<dcmplx> epsilon;
    QVector<dcmplx> epsilon11;
    QVector<dcmplx> epsilon12;
    QVector<dcmplx> epsilon13;
    QVector<dcmplx> epsilon21;
    QVector<dcmplx> epsilon22;
    QVector<dcmplx> epsilon23;
    QVector<dcmplx> epsilon31;
    QVector<dcmplx> epsilon32;
    QVector<dcmplx> epsilon33;
  
    QVector<QString> material;

    int    nread;
    int    nmatlab;
    int    advancedinterface;
    int    dipolepsilon;
    int    farfield;
    int    nearfield;
    int    force;
    int    localfield;
    int    macroscopicfield;
    int    crosssection;
    int    nenergie;
    int    crosssectionpoynting;
    int    quickdiffract;
    int    nrig;
    int    microscopy;
    int    microscopyFFT;
    int    opticalforce;
    int    opticalforcedensity;
    int    opticaltorque;
    int    opticaltorquedensity;
    int    nproche;
    int    nxx, nyy, nzz;
    int    nxm, nym, nzm;
    int    nxmp, nymp, nzmp;
    int    ntheta, nphi;
    int    nfft2d;
    double na,nainc,gross,zlens,meshsize;
    int ntypemic,nside;
    
    QVector<QColor> *colors;

};
#endif
