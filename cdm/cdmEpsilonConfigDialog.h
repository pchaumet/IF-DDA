#ifndef EPSILONCONFIGDIALOG_H
#define EPSILONCONFIGDIALOG_H

#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif
#include "QsLog.h"
#include "cdmOptions.h"

class EpsilonConfigDialog : public QDialog 
{
  
  Q_OBJECT
    
    public:
  
  EpsilonConfigDialog(QWidget *parent = 0, Qt::WindowFlags fl = Qt::SubWindow,
		      Options *_options = 0);
  virtual ~EpsilonConfigDialog();

 public slots:

  void ok();
  void onMaterialChanged(int index);

  signals:

  void handleSelectionChanged(int index);

 protected:
  void closeEvent(QCloseEvent *event);
  void showWarning(QString message);

 private:

  QVBoxLayout *layout;

  void         updateOptions();

  Options      *options;

  QVector<QLineEdit*>    epsilonr;
  QVector<QLineEdit*>    epsiloni;
  QVector<QLineEdit*>    epsilon11r;
  QVector<QLineEdit*>    epsilon11i;
  QVector<QLineEdit*>    epsilon12r;
  QVector<QLineEdit*>    epsilon12i;
  QVector<QLineEdit*>    epsilon13r;
  QVector<QLineEdit*>    epsilon13i;
  QVector<QLineEdit*>    epsilon21r;
  QVector<QLineEdit*>    epsilon21i;
  QVector<QLineEdit*>    epsilon22r;
  QVector<QLineEdit*>    epsilon22i;
  QVector<QLineEdit*>    epsilon23r;
  QVector<QLineEdit*>    epsilon23i;
  QVector<QLineEdit*>    epsilon31r;
  QVector<QLineEdit*>    epsilon31i;
  QVector<QLineEdit*>    epsilon32r;
  QVector<QLineEdit*>    epsilon32i;
  QVector<QLineEdit*>    epsilon33r;
  QVector<QLineEdit*>    epsilon33i;
  QVector<QComboBox*>    material;

  QTabWidget             *tabepsilonconfig;
  QPushButton            *okButton;
};



#endif
