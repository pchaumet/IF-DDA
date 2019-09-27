#ifndef OBJECTCONFIGDIALOG_H
#define OBJECTCONFIGDIALOG_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif

#include "QsLog.h"
#include "cdmOptions.h"

class ObjectConfigDialog : public QDialog 
{
  
  Q_OBJECT
    
    public:
  
  ObjectConfigDialog( QWidget *parentwidget = 0, Qt::WindowFlags fl = Qt::SubWindow,
		      Options *_options = 0);
  virtual ~ObjectConfigDialog();

 public slots:

  void ok();

 protected:
  void closeEvent(QCloseEvent *event);
  void showWarning(QString message);

 private:

  void                  updateOptions();

  Options               *options;
  QVector<QLineEdit*>    objectfile;
  QVector<QLineEdit*>    sphereradius;
  QVector<QLineEdit*>    spherecoherencelength;
  QVector<QLineEdit*>    spherestandardev;
  QVector<QLineEdit*>    sphereseed;
  QVector<QLineEdit*>    density;
  QVector<QLineEdit*>    cubeside;
  QVector<QLineEdit*>    cubesidex;
  QVector<QLineEdit*>    cubesidey;
  QVector<QLineEdit*>    cubesidez;
  QVector<QLineEdit*>    positionx;
  QVector<QLineEdit*>    positiony;
  QVector<QLineEdit*>    positionz;
  QVector<QLineEdit*>    demiaxea;
  QVector<QLineEdit*>    demiaxeb;
  QVector<QLineEdit*>    demiaxec;
  QVector<QLineEdit*>    thetaobj;
  QVector<QLineEdit*>    phiobj;
  QVector<QLineEdit*>    psiobj;
  QVector<QLineEdit*>    hauteur;
  QVector<QLineEdit*>    meshsize;
  QVector<QLineEdit*>    nxx;
  QVector<QLineEdit*>    nyy;
  QVector<QLineEdit*>    nzz;

  
  QTabWidget            *tabobjectconfig;
  QPushButton           *okButton;
};



#endif
