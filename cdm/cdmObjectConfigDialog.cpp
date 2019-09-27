#include "cdmObjectConfigDialog.h"

ObjectConfigDialog::ObjectConfigDialog( QWidget *parent,Qt::WindowFlags fl,
					Options *_options)
: QDialog(parent,fl)
{
   this->setWindowTitle(tr("Object properties"));
   options = _options;

   QVBoxLayout *layout = new QVBoxLayout(this);

   okButton = new QPushButton("Ok");
   connect(okButton, SIGNAL(clicked()), this, SLOT(ok()));
   layout->addWidget(okButton);

   tabobjectconfig = new QTabWidget();  
   layout->addWidget(tabobjectconfig);
   
   if (options->getObject() != "multiple spheres" && 
       options->getObject() != "concentric spheres")
	 options->setObjectNumber(1);
   QLOG_DEBUG() << "ObjectConfigDialog::ObjectConfigDialog> Object Number = " 
  	       << options->getObjectNumber()
               << "ObjectConfigDialog::ObjectConfigDialog> Object radius = "
               <<  options->getSphereradius().at(0);
   for (int i = 0 ; i < options->getObjectNumber(); i++) {
     QFormLayout *objectconfiglayout = new QFormLayout();
     QWidget *objectconfigwidget = new QWidget();
     objectconfigwidget->setLayout(objectconfiglayout);
     tabobjectconfig->setCurrentIndex(tabobjectconfig->addTab(objectconfigwidget,
		"object " + 
		QString::number(i+1)));
     objectfile.push_back(new QLineEdit(options->getObjectFile()));
     sphereradius.push_back(new QLineEdit(QString::number(options->getSphereradius().at(i))));
     spherecoherencelength.push_back(new QLineEdit(QString::number(options->getSpherecoherencelength())));
     spherestandardev.push_back(new QLineEdit(QString::number(options->getSpherestandardev())));
     sphereseed.push_back(new QLineEdit(QString::number(options->getSphereseed())));
     density.push_back(new QLineEdit(QString::number(options->getDensity())));
     cubeside.push_back(new QLineEdit(QString::number(options->getCubeside())));
     cubesidex.push_back(new QLineEdit(QString::number(options->getCubesidex())));
     cubesidey.push_back(new QLineEdit(QString::number(options->getCubesidey())));
     cubesidez.push_back(new QLineEdit(QString::number(options->getCubesidez())));
     positionx.push_back(new QLineEdit(QString::number(options->getPositionx().at(i))));
     positiony.push_back(new QLineEdit(QString::number(options->getPositiony().at(i))));
     positionz.push_back(new QLineEdit(QString::number(options->getPositionz().at(i))));
     demiaxea.push_back(new QLineEdit(QString::number(options->getDemiaxea())));
     demiaxeb.push_back(new QLineEdit(QString::number(options->getDemiaxeb())));
     demiaxec.push_back(new QLineEdit(QString::number(options->getDemiaxec())));
     psiobj.push_back(new QLineEdit(QString::number(options->getPsiobj())));
     thetaobj.push_back(new QLineEdit(QString::number(options->getThetaobj())));
     phiobj.push_back(new QLineEdit(QString::number(options->getPhiobj())));
     hauteur.push_back(new QLineEdit(QString::number(options->getHauteur())));
     nxx.push_back(new QLineEdit(QString::number(options->getNxx())));
     nyy.push_back(new QLineEdit(QString::number(options->getNyy())));
     nzz.push_back(new QLineEdit(QString::number(options->getNzz())));
     meshsize.push_back(new QLineEdit(QString::number(options->getMeshsize())));
     
     if (options->getObject() == "sphere" || options->getObject() == "multiple spheres") {
      objectconfiglayout->addRow(tr("Sphere radius (nm):"),sphereradius.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));	
     }
     else if (options->getObject() == "inhomogeneous sphere") {
      objectconfiglayout->addRow(tr("Sphere radius (nm):"),sphereradius.at(i));
      objectconfiglayout->addRow(tr("Sphere seed :"),sphereseed.at(i));
      objectconfiglayout->addRow(tr("Coherence length (nm):"),spherecoherencelength.at(i));
      objectconfiglayout->addRow(tr("Standard deviation (nm):"),spherestandardev.at(i));
     }
     else if (options->getObject() == "random spheres (length)") {
      objectconfiglayout->addRow(tr("Cuboid side x (nm):"),cubesidex.at(i));
      objectconfiglayout->addRow(tr("Cuboid side y (nm):"),cubesidey.at(i));
      objectconfiglayout->addRow(tr("Cuboid side z (nm):"),cubesidez.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("Sphere radius (nm):"),sphereradius.at(i));
      objectconfiglayout->addRow(tr("Sphere seed :"),sphereseed.at(i));
      objectconfiglayout->addRow(tr("Density :"),density.at(i));
     }
     else if (options->getObject() == "random spheres (meshsize)") {
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("number of subunit X:"),nxx.at(i));
      objectconfiglayout->addRow(tr("number of subunit Y:"),nyy.at(i));
      objectconfiglayout->addRow(tr("number of subunit Z:"),nzz.at(i));
      objectconfiglayout->addRow(tr("meshsize (nm):"),meshsize.at(i));
      objectconfiglayout->addRow(tr("Sphere radius (nm):"),sphereradius.at(i));
      objectconfiglayout->addRow(tr("Sphere seed :"),sphereseed.at(i));
      objectconfiglayout->addRow(tr("Density :"),density.at(i));
     }
     else if (options->getObject() == "concentric spheres") {
      objectconfiglayout->addRow(tr("Sphere radius (nm):"),sphereradius.at(i));
      if (i == 0) {
        objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
        objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
        objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      }
     }
     else if (options->getObject() == "cube") {
      objectconfiglayout->addRow(tr("Cube side (nm):"),cubeside.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));	
     }
     else if (options->getObject() == "cuboid (length)") {
      objectconfiglayout->addRow(tr("Cuboid side x (nm):"),cubesidex.at(i));
      objectconfiglayout->addRow(tr("Cuboid side y (nm):"),cubesidey.at(i));
      objectconfiglayout->addRow(tr("Cuboid side z (nm):"),cubesidez.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("Psi(Z):"),psiobj.at(i));
      objectconfiglayout->addRow(tr("Theta(X):"),thetaobj.at(i));
      objectconfiglayout->addRow(tr("Phi(Z):"),phiobj.at(i));
     }
     else if (options->getObject() == "cuboid (meshsize)") {
      objectconfiglayout->addRow(tr("Position X:"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y:"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z:"),positionz.at(i));
      objectconfiglayout->addRow(tr("number of subunit X:"),nxx.at(i));
      objectconfiglayout->addRow(tr("number of subunit Y:"),nyy.at(i));
      objectconfiglayout->addRow(tr("number of subunit Z:"),nzz.at(i));
      objectconfiglayout->addRow(tr("meshsize (nm):"),meshsize.at(i));
     }
     else if (options->getObject() == "inhomogeneous cuboid (length)") {
      objectconfiglayout->addRow(tr("Cuboid side x (nm):"),cubesidex.at(i));
      objectconfiglayout->addRow(tr("Cuboid side y (nm):"),cubesidey.at(i));
      objectconfiglayout->addRow(tr("Cuboid side z (nm):"),cubesidez.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("Cuboid seed :"),sphereseed.at(i));
      objectconfiglayout->addRow(tr("Coherence length (nm):"),spherecoherencelength.at(i));
      objectconfiglayout->addRow(tr("Standard deviation:"),spherestandardev.at(i));
     }
     else if (options->getObject() == "inhomogeneous cuboid (meshsize)") {
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("number of subunit X:"),nxx.at(i));
      objectconfiglayout->addRow(tr("number of subunit Y:"),nyy.at(i));
      objectconfiglayout->addRow(tr("number of subunit Z:"),nzz.at(i));
      objectconfiglayout->addRow(tr("meshsize (nm):"),meshsize.at(i));
      objectconfiglayout->addRow(tr("Cuboid seed :"),sphereseed.at(i));
      objectconfiglayout->addRow(tr("Coherence length (nm):"),spherecoherencelength.at(i));
      objectconfiglayout->addRow(tr("Standard deviation:"),spherestandardev.at(i));
     }
     else if (options->getObject() == "ellipsoid") {
      objectconfiglayout->addRow(tr("Half axe A (nm):"),demiaxea.at(i));
      objectconfiglayout->addRow(tr("Half axe B (nm):"),demiaxeb.at(i));
      objectconfiglayout->addRow(tr("Half axe C (nm):"),demiaxec.at(i));	
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("Psi(Z):"),psiobj.at(i));	
      objectconfiglayout->addRow(tr("Theta(X):"),thetaobj.at(i));
      objectconfiglayout->addRow(tr("Phi(Z):"),phiobj.at(i));	
     }
     else if (options->getObject() == "cylinder") {
      objectconfiglayout->addRow(tr("cylinder radius (nm):"),sphereradius.at(i));
      objectconfiglayout->addRow(tr("Cylinder height (nm):"),hauteur.at(i));
      objectconfiglayout->addRow(tr("Position X (nm):"),positionx.at(i));
      objectconfiglayout->addRow(tr("Position Y (nm):"),positiony.at(i));
      objectconfiglayout->addRow(tr("Position Z (nm):"),positionz.at(i));
      objectconfiglayout->addRow(tr("Psi(Z):"),psiobj.at(i));	
      objectconfiglayout->addRow(tr("Theta(X):"),thetaobj.at(i));
      objectconfiglayout->addRow(tr("Phi(Z):"),phiobj.at(i));
     }
     else if (options->getObject() == "arbitrary") {
       objectconfiglayout->addRow(tr("File path:"),objectfile.at(i));
     }
   }
}
ObjectConfigDialog::~ObjectConfigDialog()
{
  QLOG_DEBUG ( ) << "Deleting ObjectConfigDialog";
}
void
ObjectConfigDialog::ok()
{
  this->updateOptions();
  this->close();
}
void
ObjectConfigDialog::updateOptions()
{
   QVector<double> _sphereradius;
   QVector<double> _positionx;
   QVector<double> _positiony;
   QVector<double> _positionz;

   QLOG_DEBUG() << "ObjectConfigDialog::updateOptions> Object Number = " 
	       << options->getObjectNumber();

   for (int i = 0 ; i < options->getObjectNumber(); i++) {
     if (options->getObject() == "concentric spheres") {
         _sphereradius.push_back(sphereradius.at(i)->text().toDouble());
         _positionx.push_back(positionx.at(0)->text().toDouble());
         _positiony.push_back(positiony.at(0)->text().toDouble());
         _positionz.push_back(positionz.at(0)->text().toDouble());
     }
     else {
        options->setSphereseed(sphereseed.at(i)->text().toInt());
        options->setSpherecoherencelength(spherecoherencelength.at(i)->text().toDouble());
        options->setSpherestandardev(spherestandardev.at(i)->text().toDouble());
	_sphereradius.push_back(sphereradius.at(i)->text().toDouble());
	_positionx.push_back(positionx.at(i)->text().toDouble());
	_positiony.push_back(positiony.at(i)->text().toDouble());
	_positionz.push_back(positionz.at(i)->text().toDouble());
        options->setObjectFile(objectfile.at(i)->text());
	options->setDensity(density.at(i)->text().toDouble());
	options->setCubeside(cubeside.at(i)->text().toDouble());
	options->setCubesidex(cubesidex.at(i)->text().toDouble());
	options->setCubesidey(cubesidey.at(i)->text().toDouble());
	options->setCubesidez(cubesidez.at(i)->text().toDouble());
	options->setThetaobj(thetaobj.at(i)->text().toDouble());
	options->setPhiobj(phiobj.at(i)->text().toDouble());
	options->setPsiobj(psiobj.at(i)->text().toDouble());
	options->setDemiaxea(demiaxea.at(i)->text().toDouble());
	options->setDemiaxeb(demiaxeb.at(i)->text().toDouble());
	options->setDemiaxec(demiaxec.at(i)->text().toDouble());
	options->setHauteur(hauteur.at(i)->text().toDouble());
	options->setMeshsize(meshsize.at(i)->text().toDouble());
	options->setNxx(nxx.at(i)->text().toInt());
	options->setNyy(nyy.at(i)->text().toInt());
	options->setNzz(nzz.at(i)->text().toInt());
     }
   }
   options->setSphereradius(_sphereradius);
   options->setPositionx(_positionx);
   options->setPositiony(_positiony);
   options->setPositionz(_positionz);
}
void 
ObjectConfigDialog::closeEvent(QCloseEvent* event)
{
  event->accept();  
  QLOG_DEBUG ( ) << "Closing ObjectConfigDialog";
  this->close();
}
void 
ObjectConfigDialog::showWarning(QString message) {
  QMessageBox::warning(this, "Error:", message);
}
