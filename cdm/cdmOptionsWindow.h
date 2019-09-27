#ifndef OPTIONSWINDOW_H
#define OPTIONSWINDOW_H

#include <QMessageBox>
#include <QtSql>
#include <QtGui>
#if QT_VERSION < 0x050000
#else
#include <QtWidgets>
#endif


#include "cdmOptionsWidget.h"
#include "QsLog.h"

class ManifestModel: public QSqlTableModel
{
    public:
    ManifestModel(QObject * parent = 0, QSqlDatabase db = QSqlDatabase() ):
    QSqlTableModel(parent, db)
    { }
     
    ~ManifestModel() { }
     
    Qt::ItemFlags flags ( const QModelIndex & index ) const
    {
      if ((index.column() == 0 && index.data() == "") || index.column() == 1) 
         return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
      else
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    }
};
class ColorDelegate: public QItemDelegate
{
    public:
    ColorDelegate(QObject *parent = 0) : QItemDelegate(parent) {}

    public:
    virtual void paint(QPainter *painter, const QStyleOptionViewItem &option, 
			const QModelIndex &index) const
    {
       QItemDelegate::paint(painter, option, index);
       if(option.state & QStyle::State_Selected)
        {
         drawBackground(painter, option, index);
        }
    }
    protected:
    virtual void drawBackground(QPainter *painter, const QStyleOptionViewItem &option, 
			const QModelIndex &index) const
    {
       if ((index.column() == 0 && index.data() == "") || index.column() == 1) {
        painter->fillRect(option.rect, Qt::lightGray);
        painter->save();
        QPen pen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap, Qt::MiterJoin);
        int w = pen.width()/2;
        painter->setPen(pen);
        painter->drawRect(option.rect.adjusted(w,w,-w,-w));
        painter->restore();
      }
   }
};

class OptionsWindow : public QMainWindow 
{
  
  Q_OBJECT
    
    public:
  
  OptionsWindow( QMainWindow* _parent = 0, Options *_options = 0, Qt::WindowFlags fl = Qt::Window);
  virtual ~OptionsWindow();
 
  void  initOptionsTbl();

  signals:
   void updateOptionsWidget();

 public slots:
   void save();
   void remove();
   void load();
   void tofile();
   void updatewindow();

 protected:
  void closeEvent(QCloseEvent *event);
  void showWarning(QString message);

 private:

  QMainWindow     *parent;
  OptionsWidget   *optionswidget;
  Options         *options;

  QString         dbpath;

  QLabel          *optionstitle;
  ManifestModel   *optionstable;
  QTableView      *optionsview;

  QVBoxLayout*    vboxlayout;
   
  QPushButton*    saveButton;
  QPushButton*    removeButton;
  QPushButton*    loadButton;
  QPushButton*    exportButton;
};
#endif // OPTIONSWINDOW_H
