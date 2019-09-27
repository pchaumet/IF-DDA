#ifndef BEAMCONFIGDIALOG_H
#define BEAMCONFIGDIALOG_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif
#include "QsLog.h"
#include "cdmOptions.h"

class BeamConfigDialog : public QDialog
{
  
  Q_OBJECT
    
    public:
  
  BeamConfigDialog( QWidget *_parentwidget = 0, Qt::WindowFlags fl = Qt::SubWindow,
		    Options *_options = 0);
  virtual ~BeamConfigDialog();

 public slots:
  void ok();

 protected:
  void closeEvent(QCloseEvent *event);
  void showWarning(QString message);

 private:

  void         updateOptions();

  Options      *options;
  QLineEdit    *beamfile;
  QLineEdit    *incidenceangle_theta_z;
  QLineEdit    *incidenceangle_phi_x;
  QLineEdit    *polarizationTM;
  QLineEdit    *polarizationTE;
  double        polarizationRL;
  QRadioButton *positiveRadio;
  QRadioButton *negativeRadio;
  QLineEdit    *xgaus;
  QLineEdit    *ygaus;
  QLineEdit    *zgaus;
  QVector<QLineEdit*>   thetam;
  QVector<QLineEdit*>   phim;
  QVector<QLineEdit*>   ppm;
  QVector<QLineEdit*>   ssm;
  QVector<QLineEdit*>   E0mr;
  QVector<QLineEdit*>   E0mi;
  QTabWidget   *tabwavemulticonfig;
  QPushButton  *okButton;
};



#endif
