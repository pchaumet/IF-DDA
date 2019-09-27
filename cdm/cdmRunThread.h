#ifndef RUNTHREAD_H
#define RUNTHREAD_H

#include <QtSql>
#include <QtGui>
#include "QsLog.h"
#include "cdmOptions.h"
#include "cdmlib.h"

class RunThread : public QThread
{
  Q_OBJECT
    
    public:
  RunThread(QObject* parent = 0, Options *_options);
  ~RunThread();

  QMutex *mutex;

  void stop();
  

 signals:

 protected:
  virtual void run();

 private:
  bool    suspend;
  Options *options;
 
};

#endif // RUNTHREAD_H
