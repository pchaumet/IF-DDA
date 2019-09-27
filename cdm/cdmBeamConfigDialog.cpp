#include "cdmBeamConfigDialog.h"

BeamConfigDialog::BeamConfigDialog(QWidget *parent,Qt::WindowFlags fl,
				   Options *_options)
: QDialog(parent,fl)
{
   options = _options;
   beamfile = new QLineEdit(options->getBeamFile());
   incidenceangle_theta_z = new QLineEdit(QString::number(options->getIncidenceangle_theta_z()));
   incidenceangle_phi_x= new QLineEdit(QString::number(options->getIncidenceangle_phi_x()));
   polarizationTM= new QLineEdit(QString::number(options->getPolarizationTM()));
   polarizationTE= new QLineEdit(QString::number(options->getPolarizationTE()));
   polarizationRL= options->getPolarizationRL();
   positiveRadio = new QRadioButton("Right");
   negativeRadio = new QRadioButton("Left");
   QBoxLayout *polarizationlayout = new QBoxLayout(QBoxLayout::LeftToRight);
   polarizationlayout->addWidget(positiveRadio);
   polarizationlayout->addWidget(negativeRadio);
   if ( polarizationRL > 0 )
     positiveRadio->setChecked(true);
   else
     negativeRadio->setChecked(true);
   xgaus = new QLineEdit(QString::number(options->getXgaus()));
   ygaus = new QLineEdit(QString::number(options->getYgaus()));
   zgaus = new QLineEdit(QString::number(options->getZgaus()));

   tabwavemulticonfig = new QTabWidget();

   okButton = new QPushButton("Ok");
   connect(okButton, SIGNAL(clicked()), this, SLOT(ok()));

   QFormLayout *layout = new QFormLayout(this);
   
   if (options->getBeam() == "Circular plane wave") {
    layout->addRow(tr("Incidence angle (theta with respect to z):"),incidenceangle_theta_z);
    layout->addRow(tr("Incidence angle (phi with respect to x):"),incidenceangle_phi_x);
    layout->addRow(tr("polarization:"),polarizationlayout);
   }
   else if (options->getBeam() == "Linear plane wave") {
    layout->addRow(tr("Incidence angle (theta with respect to z):"),incidenceangle_theta_z);
    layout->addRow(tr("Incidence angle (phi with respect to x):"),incidenceangle_phi_x);
    layout->addRow(tr("Polarization TM (1) TE (0):"),polarizationTM);
    //layout->addRow(tr("polarization (TE with respect to the plane x,y):"),polarizationTE);
   }
   else if (options->getBeam() == "Multiplane wave") {

    layout->addWidget(tabwavemulticonfig);
  
    for (int i = 0 ; i < options->getWaveMultiNumber(); i++) {

       QFormLayout *wavemulticonfiglayout = new QFormLayout();

       thetam.push_back(new QLineEdit(QString::number(options->getThetam().at(i))));
       phim.push_back(new QLineEdit(QString::number(options->getPhim().at(i))));
       ppm.push_back(new QLineEdit(QString::number(options->getPpm().at(i))));
       //       ssm.push_back(new QLineEdit(QString::number(options->getSsm().at(i))));
       E0mr.push_back(new QLineEdit(QString::number(real(options->getE0m().at(i)))));
       E0mi.push_back(new QLineEdit(QString::number(imag(options->getE0m().at(i)))));
         
       QBoxLayout *thetamlayout = new QBoxLayout(QBoxLayout::LeftToRight);
       thetamlayout->addWidget(thetam.at(i));
       QBoxLayout *phimlayout = new QBoxLayout(QBoxLayout::LeftToRight);
       phimlayout->addWidget(phim.at(i));
       QBoxLayout *ppmlayout = new QBoxLayout(QBoxLayout::LeftToRight);
       ppmlayout->addWidget(ppm.at(i));
       //       QBoxLayout *ssmlayout = new QBoxLayout(QBoxLayout::LeftToRight);
       // ssmlayout->addWidget(ssm.at(i));
       QBoxLayout *E0mlayout = new QBoxLayout(QBoxLayout::LeftToRight);
       E0mlayout->addWidget(E0mr.at(i));
       E0mlayout->addWidget(E0mi.at(i));
    
       wavemulticonfiglayout->addRow("Incidence angle (theta with respect to z):",thetamlayout);
       wavemulticonfiglayout->addRow("Incidence angle (phi with respect to x):",phimlayout);
       wavemulticonfiglayout->addRow("Polarization TM (1) TE (0):",ppmlayout);
       // wavemulticonfiglayout->addRow("Polarisation S:",ssmlayout);
       wavemulticonfiglayout->addRow("Amplitude (complex):",E0mlayout);
       
       QWidget *wavemulticonfigwidget = new QWidget();
       wavemulticonfigwidget->setLayout(wavemulticonfiglayout);
       tabwavemulticonfig->setCurrentIndex(tabwavemulticonfig->addTab(wavemulticonfigwidget,
		"plane wave " + 
		QString::number(i+1)));
      }
   }
   else if (options->getBeam() == "Linear Gaussian" || options->getBeam() == "Linear Gaussian (para)" 
	                                        || options->getBeam() == "Linear Gaussian (FFT)" ) {
    layout->addRow(tr("Incidence angle (theta with respect to z):"),incidenceangle_theta_z);
    layout->addRow(tr("Incidence angle (phi with respect to x):"),incidenceangle_phi_x);
    layout->addRow(tr("polarization TM (1) TE (0):"),polarizationTM);
    layout->addRow(tr("Gaussian X (nm):"),xgaus);
    layout->addRow(tr("Gaussian Y (nm):"),ygaus);
    layout->addRow(tr("Gaussian Z (nm):"),zgaus);
   }
   else if (options->getBeam() == "Antenna" ) {
     layout->addRow(tr("Orientation angle (theta with respect to z):"),incidenceangle_theta_z);
    layout->addRow(tr("Orientation angle (phi with respect to x):"),incidenceangle_phi_x);
    layout->addRow(tr("position X (nm):"),xgaus);
    layout->addRow(tr("position Y (nm):"),ygaus);
    layout->addRow(tr("position Z (nm):"),zgaus);
   }
   else if (options->getBeam() == "Circular Gaussian" || options->getBeam() == "Circular Gaussian (para)"
	                                          ||  options->getBeam() == "Circular Gaussian (FFT)" ) {
    layout->addRow(tr("Incidence angle (theta with respect to z):"),incidenceangle_theta_z);
    layout->addRow(tr("Incidence angle (phi with respect to x):"),incidenceangle_phi_x);
    layout->addRow(tr("polarization:"),polarizationlayout);
    layout->addRow(tr("Gaussian X:"),xgaus);
    layout->addRow(tr("Gaussian Y:"),ygaus);
    layout->addRow(tr("Gaussian Z:"),zgaus);
   }
   else if (options->getBeam() == "Arbitrary wave (file)") {
    layout->addRow(tr("File path:"),beamfile);
   }
   else {
    layout->addRow(tr(""),new QLabel("  Beam type not defined yet  "));
   }
   layout->addRow("",okButton);
   this->setWindowTitle(tr("Beam properties"));
}
BeamConfigDialog::~BeamConfigDialog()
{
  QLOG_DEBUG ( ) << "Deleting BeamConfigDialog";
}
void
BeamConfigDialog::ok()
{
  this->updateOptions();
  this->close();
}
void
BeamConfigDialog::updateOptions()
{
   options->setBeamFile(beamfile->text());
   options->setIncidenceangle_theta_z(incidenceangle_theta_z->text().toDouble());
   options->setIncidenceangle_phi_x(incidenceangle_phi_x->text().toDouble());
   options->setPolarizationTM(polarizationTM->text().toDouble());
   //options->setPolarizationTE(polarizationTE->text().toDouble());
   if ( positiveRadio->isChecked() == true )
    polarizationRL = 1;
   else if ( negativeRadio->isChecked() == true )
    polarizationRL = -1;
   options->setPolarizationRL(polarizationRL);
   options->setXgaus(xgaus->text().toDouble());
   options->setYgaus(ygaus->text().toDouble());
   options->setZgaus(zgaus->text().toDouble());
   
   if (options->getBeam() == "Multiplane wave") {
	   QVector<double> _thetam;
	   QVector<double> _phim;
	   QVector<double> _ppm;
	   // QVector<double> _ssm;
	   QVector<dcmplx> _E0m;
	   for (int i = 0 ; i < options->getWaveMultiNumber(); i++) {
	    _thetam.push_back(thetam.at(i)->text().toDouble());
	    _phim.push_back(phim.at(i)->text().toDouble());
	    _ppm.push_back(ppm.at(i)->text().toDouble());
	    //_ssm.push_back(ssm.at(i)->text().toDouble());
	    _E0m.push_back(dcmplx(E0mr.at(i)->text().toDouble(),E0mi.at(i)->text().toDouble()));
	   }
	   options->setThetam(_thetam);
	   options->setPhim(_phim);
	   options->setPpm(_ppm);
	   //	   options->setSsm(_ssm);
	   options->setE0m(_E0m);
  }
}
void 
BeamConfigDialog::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing BeamConfigDialog";
  this->close();
}
void 
BeamConfigDialog::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}
