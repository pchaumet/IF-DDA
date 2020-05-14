#ifndef OPTIONSWIDGET_H
#define OPTIONSWIDGET_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif

#include "QsLog.h"
#include "cdmOptions.h"
#include "cdmRunWidget.h"
#include "cdmBeamConfigDialog.h"
#include "cdmObjectConfigDialog.h"
#include "cdmEpsilonConfigDialog.h"

class OptionsWidget : public QWidget
{
  Q_OBJECT
    
    public:
  OptionsWidget(QMainWindow *_mainwindow = 0, Options *_options = 0);
  ~OptionsWidget();

   QString  getFilereread();
   QString  getH5File();
   double   getWavelength();
   double   getP0();
   double   getW0();
   QString  getBeam();
   QString  getObject();
   int      getObjectNumber();
   int      getWaveMultiNumber();
   QString  getAnisotropy();
   int      getDiscretization();
   double   getTolerance();   
   QString  getMethodeit();
   QString  getPolarizability();
   int      getQuad();
   int      getnfft2d();

   QMainWindow *mainwindow;
   Options  *options;
   BeamConfigDialog *beamconfigdlg;
   ObjectConfigDialog *objectconfigdlg;
   EpsilonConfigDialog *epsilonconfigdlg;

   void setWavelength(double _wavelength);
   void setP0(double _P0);
   void setW0(double _W0);
   void setBeam(QString _beam);
   void setObject(QString _beam);
   void setObjectNumber(int _objectnumber);
   void setWaveMultiNumber(int _wavemultinumber);
   void setAnisotropy(QString _anisotropy);
   void setDiscretization(int _discretization);
   void setTolerance(double _tolerance);
   void setMethodeit(QString _methodeit);
   void setH5file(QString _fileh5);
   void setPolarizability(QString _polarizability);
   void setQuad(int _quad);
   void setnfft2d(int _nfft2d);
   void update();
   void updateFarfield();
   void updateAdvancedinterface();
   void updateMicroscopy();
   void updateNearfield();
   void updateForce();
   void updateOptions();

    signals:
      void updateOptionsWindow();

    public slots:
      void execute();
      void saveName();
      void saveAsName();
      void finish();
      void cancel();
      void handleObjectSelectionChanged(int index);
      void handleObjectNumberSelectionChanged(int index);
      void handleBeamSelectionChanged(int index);
      void handleWaveMultiNumberSelectionChanged(int index);
      void handleRangeofstudySelectionChanged(int index);  
      void configureBeam();
      void configureObject();
      void configureEpsilon();
      void nreadCheckBoxStateChanged(int state);
      void advancedinterfaceCheckBoxStateChanged(int state);
      void farfieldCheckBoxStateChanged(int state);
      void nenergieCheckBoxStateChanged(int state);
      void crosssectionCheckBoxStateChanged(int state);
      void microscopyFFTCheckBoxStateChanged(int state);      
      void crosssectionpoyntingCheckBoxStateChanged(int state);
      void microscopyCheckBoxStateChanged(int state);
      void ntypemicStateChanged(int index);
      void nsideStateChanged(int index);
      void nmatlabStateChanged(int index);
      void forceCheckBoxStateChanged(int state);
      void opticalforceCheckBoxStateChanged(int state);
      void opticalforcedensityCheckBoxStateChanged(int state);
      void opticaltorqueCheckBoxStateChanged(int state);
      void opticaltorquedensityCheckBoxStateChanged(int state);
      void nearfieldCheckBoxStateChanged(int state);
      void localfieldCheckBoxStateChanged(int state);
      void macroscopicfieldCheckBoxStateChanged(int state);
      void dipolepsilonCheckBoxStateChanged(int state);

    private:

    QString         name, description;
    QFormLayout     *layout;
    QPushButton     *executeButton, *saveButton;
    QLabel          *emptynenergieLabel;
    QLabel          *emptycrosssectionLabel, *emptycrosssectionpoyntingLabel;
    QLabel          *emptylocalfieldLabel, *emptymacroscopicfieldLabel, *emptyrangeofstudyLabel;
    QLabel          *emptyopticalforceLabel, *emptyopticalforcedensityLabel;
    QLabel          *emptyopticaltorqueLabel, *emptyopticaltorquedensityLabel;
    QLabel          *emptymicroscopyLabel, *emptynaLabel, *emptygrossLabel, *emptymicroscopyFFTLabel;
    QLabel          *emptyzlensLabel, *emptynaincLabel;
    QLabel          *wavelengthLabel;
    QLineEdit       *wavelength;
    QLabel          *beamLabel;
    QLineEdit       *P0;
    QLabel          *P0Label;
    QLineEdit       *W0;
    QLabel          *W0Label;
    QComboBox       *beam;
    QPushButton     *beamButton;
    QLabel          *objectLabel;
    QComboBox       *object;
    QLabel          *objectnumberLabel;
    QSpinBox        *objectnumber;
    QLabel          *wavemultinumberLabel;
    QSpinBox        *wavemultinumber;
    QPushButton     *objectButton;
    QLabel          *anisotropyLabel;
    QComboBox       *anisotropy;
    QPushButton     *epsilonButton;
    QLabel          *discretizationLabel;
    QLineEdit       *discretization;
    QLabel          *toleranceLabel;
    QLineEdit       *tolerance;
    QLabel          *methodeitLabel;
    QComboBox       *methodeit;
    QLabel          *polarizabilityLabel;
    QComboBox       *polarizability;
    QLabel          *quadLabel;
    QComboBox       *quad;
    QLabel          *nfft2dLabel;
    QComboBox       *nfft2d;
    QLabel          *nmatlabLabel;
    QComboBox       *nmatlab;
    QLabel          *nrigLabel;
    QComboBox       *nrig;
    QLabel          *ntypemicLabel;
    QComboBox       *ntypemic;
    QLabel          *nsideLabel;
    QComboBox       *nside;
    QLabel          *localfieldLabel;
    QCheckBox       *localfield;
    QLabel          *nreadLabel;
    QCheckBox       *nread;
    QLabel          *filerereadLabel;
    QLineEdit       *filereread;
    QLabel          *fileh5Label;
    QLineEdit       *fileh5;
    QLabel          *advancedinterfaceLabel;
    QCheckBox       *advancedinterface;
    QLabel          *macroscopicfieldLabel;
    QCheckBox       *macroscopicfield;
    QLabel          *nenergieLabel;
    QCheckBox       *nenergie;
    QLabel          *crosssectionLabel;
    QCheckBox       *crosssection;
    QLabel          *crosssectionpoyntingLabel;
    QCheckBox       *crosssectionpoynting;
    QLabel          *quickdiffractLabel;
    QCheckBox       *quickdiffract;
    QLabel          *microscopyLabel;
    QCheckBox       *microscopy;
    QLabel          *microscopyFFTLabel;
    QCheckBox       *microscopyFFT;
    QLabel          *opticalforceLabel;
    QCheckBox       *opticalforce;
    QLabel          *opticalforcedensityLabel;
    QCheckBox       *opticalforcedensity;
    QLabel          *opticaltorqueLabel;
    QCheckBox       *opticaltorque;
    QLabel          *opticaltorquedensityLabel;
    QCheckBox       *opticaltorquedensity;
    QLabel          *rangeofstudyLabel;
    QComboBox       *rangeofstudy;
    QCheckBox       *dipolepsilon, *farfield, *force, *nearfield ;
    QLabel          *nthetaLabel, *nphiLabel, *naLabel, *naincLabel, *grossLabel, *zlensLabel;
    QLabel          *nxmpLabel, *nympLabel, *nzmpLabel, *emptynxmpLabel, *emptynympLabel, *emptynzmpLabel;
    QLineEdit       *nxm, *nym, *nzm, *ntheta, *nphi, *na, *gross, *zlens, *nainc, *meshsize;
    QLineEdit       *nxmp, *nymp, *nzmp, *nxx, *nyy, *nzz;
    QDialog         *cfgWindow;
    QLineEdit       *nameLineEdit, *descriptionLineEdit;
    
};
#endif
