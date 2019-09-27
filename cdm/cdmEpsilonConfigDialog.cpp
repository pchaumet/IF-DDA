#include "cdmEpsilonConfigDialog.h"


EpsilonConfigDialog::EpsilonConfigDialog(QWidget *parent,Qt::WindowFlags fl,
					 Options *_options)
: QDialog(parent,fl)
{
   this->setWindowTitle(tr("Epsilon properties"));
   options = _options;

   layout = new QVBoxLayout(this);

   okButton = new QPushButton("Ok");
   connect(okButton, SIGNAL(clicked()), this, SLOT(ok()));
   layout->addWidget(okButton);

   tabepsilonconfig = new QTabWidget();
   layout->addWidget(tabepsilonconfig);

   QSignalMapper *signalMapper = new QSignalMapper(this);

   QLOG_DEBUG() << "EpsilonConfigDialog::EpsilonConfigDialog> Object number: " 
	       << options->getObjectNumber();

   for (int i = 0 ; i < options->getObjectNumber(); i++) {
      
       QFormLayout *epsilonconfiglayout = new QFormLayout();
     
       if (options->getAnisotropy() == "iso") {
         epsilonr.push_back(new QLineEdit(QString::number(real(options->getEpsilon().at(i)))));
         epsiloni.push_back(new QLineEdit(QString::number(imag(options->getEpsilon().at(i)))));
         QComboBox *materialBox = new QComboBox(this);
         connect( materialBox, SIGNAL(currentIndexChanged(int)), signalMapper, SLOT(map()));
         signalMapper->setMapping(materialBox, i);
         materialBox->addItems(options->materialList);
         materialBox->setCurrentIndex(0);
         material.push_back(materialBox);
         
         QBoxLayout   *epsilonlayout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilonlayout->addWidget(epsilonr.at(i));
         epsilonlayout->addWidget(epsiloni.at(i));
         epsilonlayout->addWidget(material.at(i));
    
         epsilonconfiglayout->addRow("Epsilon:",epsilonlayout);
       }
       else if (options->getAnisotropy() == "ani") {
         epsilon11r.push_back(new QLineEdit(QString::number(real(options->getEpsilon11().at(i)))));
         epsilon11i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon11().at(i)))));
         epsilon12r.push_back(new QLineEdit(QString::number(real(options->getEpsilon12().at(i)))));
         epsilon12i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon12().at(i)))));
         epsilon13r.push_back(new QLineEdit(QString::number(real(options->getEpsilon13().at(i)))));
         epsilon13i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon13().at(i)))));
         epsilon21r.push_back(new QLineEdit(QString::number(real(options->getEpsilon21().at(i)))));
         epsilon21i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon21().at(i)))));
         epsilon22r.push_back(new QLineEdit(QString::number(real(options->getEpsilon22().at(i)))));
         epsilon22i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon22().at(i)))));
         epsilon23r.push_back(new QLineEdit(QString::number(real(options->getEpsilon23().at(i)))));
         epsilon23i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon23().at(i)))));
         epsilon31r.push_back(new QLineEdit(QString::number(real(options->getEpsilon31().at(i)))));
         epsilon31i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon31().at(i)))));
         epsilon32r.push_back(new QLineEdit(QString::number(real(options->getEpsilon32().at(i)))));
         epsilon32i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon32().at(i)))));
         epsilon33r.push_back(new QLineEdit(QString::number(real(options->getEpsilon33().at(i)))));
         epsilon33i.push_back(new QLineEdit(QString::number(imag(options->getEpsilon33().at(i)))));
       
         QBoxLayout   *epsilon11layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon11layout->addWidget(epsilon11r.at(i));
         epsilon11layout->addWidget(epsilon11i.at(i));
         QBoxLayout   *epsilon12layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon12layout->addWidget(epsilon12r.at(i));
         epsilon12layout->addWidget(epsilon12i.at(i));
         QBoxLayout   *epsilon13layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon13layout->addWidget(epsilon13r.at(i));
         epsilon13layout->addWidget(epsilon13i.at(i));
         QBoxLayout   *epsilon21layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon21layout->addWidget(epsilon21r.at(i));
         epsilon21layout->addWidget(epsilon21i.at(i));
         QBoxLayout   *epsilon22layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon22layout->addWidget(epsilon22r.at(i));
         epsilon22layout->addWidget(epsilon22i.at(i));
         QBoxLayout   *epsilon23layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon23layout->addWidget(epsilon23r.at(i));
         epsilon23layout->addWidget(epsilon23i.at(i));
         QBoxLayout   *epsilon31layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon31layout->addWidget(epsilon31r.at(i));
         epsilon31layout->addWidget(epsilon31i.at(i));
         QBoxLayout   *epsilon32layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon32layout->addWidget(epsilon32r.at(i));
         epsilon32layout->addWidget(epsilon32i.at(i));
         QBoxLayout   *epsilon33layout = new QBoxLayout(QBoxLayout::LeftToRight);
         epsilon33layout->addWidget(epsilon33r.at(i));
         epsilon33layout->addWidget(epsilon33i.at(i));
   
         epsilonconfiglayout->addRow("Epsilon11:",epsilon11layout);
         epsilonconfiglayout->addRow("Epsilon12:",epsilon12layout);
         epsilonconfiglayout->addRow("Epsilon13:",epsilon13layout);
         epsilonconfiglayout->addRow("Epsilon21:",epsilon21layout);
         epsilonconfiglayout->addRow("Epsilon22:",epsilon22layout);
         epsilonconfiglayout->addRow("Epsilon23:",epsilon23layout);
         epsilonconfiglayout->addRow("Epsilon31:",epsilon31layout);
         epsilonconfiglayout->addRow("Epsilon32:",epsilon32layout);
         epsilonconfiglayout->addRow("Epsilon33:",epsilon33layout);
       }

       QWidget *epsilonconfigwidget = new QWidget();
       epsilonconfigwidget->setLayout(epsilonconfiglayout);
       tabepsilonconfig->setCurrentIndex(tabepsilonconfig->addTab(epsilonconfigwidget,
		"object " + 
		QString::number(i+1)));
   }

   if (options->getAnisotropy() == "iso") {
      connect(signalMapper, SIGNAL(mapped(int)), this, SIGNAL(handleSelectionChanged(int)));
      connect(this, SIGNAL(handleSelectionChanged(int)), this, SLOT(onMaterialChanged(int)));
     // init combo index (triggers the handle)
     for (int i = 0 ; i < options->getObjectNumber(); i++) {
        QComboBox *materialBox = material.at(i);
        QLOG_DEBUG() << "EpsilonConfigDialog::EpsilonConfigDialog> Object " << i+1 
		    << " material: " <<  options->getMaterial().at(i);
        materialBox->setCurrentIndex(materialBox->findText(options->getMaterial().at(i)));
     }
   }
}
EpsilonConfigDialog::~EpsilonConfigDialog()
{
  QLOG_DEBUG ( ) << "Deleting EpsilonConfigDialog";
}
void 
EpsilonConfigDialog::onMaterialChanged(int index) {
 QLOG_DEBUG() << "EpsilonConfigDialog::onMaterialChanged";
 QComboBox *materialBox = material.at(index);
 if ( materialBox->currentText() != "xx" )  {
   QLineEdit *epsilonrEdit = epsilonr.at(index);
   QLineEdit *epsiloniEdit = epsiloni.at(index);
   epsilonrEdit->setEnabled( false );
   epsiloniEdit->setEnabled( false );
 }
 else {
   QLineEdit *epsilonrEdit = epsilonr.at(index);
   QLineEdit *epsiloniEdit = epsiloni.at(index);
   epsilonrEdit->setEnabled( true );
   epsiloniEdit->setEnabled( true );
 }
}
void
EpsilonConfigDialog::ok()
{
  this->updateOptions();
  this->close();
}
void
EpsilonConfigDialog::updateOptions()
{
   QVector<dcmplx> _epsilon;
   QVector<dcmplx> _epsilon11;
   QVector<dcmplx> _epsilon12;
   QVector<dcmplx> _epsilon13;
   QVector<dcmplx> _epsilon21;
   QVector<dcmplx> _epsilon22;
   QVector<dcmplx> _epsilon23;
   QVector<dcmplx> _epsilon31;
   QVector<dcmplx> _epsilon32;
   QVector<dcmplx> _epsilon33;
   QVector<QString> _material;

    for (int i = 0 ; i < options->getObjectNumber(); i++) {
      if (options->getAnisotropy() == "iso") {
       _epsilon.push_back(dcmplx(epsilonr.at(i)->text().toDouble(),epsiloni.at(i)->text().toDouble()));
       QLOG_DEBUG() << " material object " << i << " : " <<  material.at(i)->currentText();
       _material.push_back(material.at(i)->currentText());
      }
      else if (options->getAnisotropy() == "ani") {
       _epsilon11.push_back(dcmplx(epsilon11r.at(i)->text().toDouble(),epsilon11i.at(i)->text().toDouble()));
       _epsilon12.push_back(dcmplx(epsilon12r.at(i)->text().toDouble(),epsilon12i.at(i)->text().toDouble()));
       _epsilon13.push_back(dcmplx(epsilon13r.at(i)->text().toDouble(),epsilon13i.at(i)->text().toDouble()));
       _epsilon21.push_back(dcmplx(epsilon21r.at(i)->text().toDouble(),epsilon21i.at(i)->text().toDouble()));
       _epsilon22.push_back(dcmplx(epsilon22r.at(i)->text().toDouble(),epsilon22i.at(i)->text().toDouble()));
       _epsilon23.push_back(dcmplx(epsilon23r.at(i)->text().toDouble(),epsilon23i.at(i)->text().toDouble()));
       _epsilon31.push_back(dcmplx(epsilon31r.at(i)->text().toDouble(),epsilon31i.at(i)->text().toDouble()));
       _epsilon32.push_back(dcmplx(epsilon32r.at(i)->text().toDouble(),epsilon32i.at(i)->text().toDouble()));
       _epsilon33.push_back(dcmplx(epsilon33r.at(i)->text().toDouble(),epsilon33i.at(i)->text().toDouble()));
     }
   }
   if (options->getAnisotropy() == "iso") {
     options->setEpsilon(_epsilon);
     options->setMaterial(_material);
   }
   else if (options->getAnisotropy() == "ani") {
     options->setEpsilon11(_epsilon11);
     options->setEpsilon12(_epsilon12);
     options->setEpsilon13(_epsilon13);
     options->setEpsilon21(_epsilon21);
     options->setEpsilon22(_epsilon22);
     options->setEpsilon23(_epsilon23);
     options->setEpsilon31(_epsilon31);
     options->setEpsilon32(_epsilon32);
     options->setEpsilon33(_epsilon33);
   } 
}
void 
EpsilonConfigDialog::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing EpsilonConfigDialog";
  this->close();
}
void 
EpsilonConfigDialog::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}
