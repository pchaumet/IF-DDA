#include "cdmPlot.h"

class MyQwtSymbol: public QwtSymbol
{
public:
    MyQwtSymbol( const QBrush &brush , QwtPointPolar datap)
    {
        QPen pen( Qt::black, 1.5);
        pen.setJoinStyle( Qt::RoundJoin );
        pen.setCosmetic( true );
        QPainterPath path = createArrow( datap );
        //setPinPoint( QPointF( 0.0, 0.0 ) );
        setPen( pen );
        setBrush( brush );
        setPath( path );    
    }

private:
    QPainterPath createArrow( QwtPointPolar datap ) const
    {
        QPainterPath path; 
        const double h = datap.radius();
        const double a = datap.radius() / 10;
        QLOG_DEBUG() << "MyQwtSymbol::createArrow:radius:" << 
				QString::number(datap.radius());
        QLOG_DEBUG() << "MyQwtSymbol::createArrow:azimuth:" << 
				QString::number(datap.azimuth());
        
        path.lineTo( h - 3 * a , 0 );
        path.moveTo( h, 0 );
        path.lineTo( h - 3 * a, a );
        path.lineTo( h - 3 * a, -a );
        path.lineTo( h, 0 );
      
        QTransform transform;
        transform.rotate( datap.azimuth() * (-180.) / M_PI );
        path = transform.map( path );

        return path;
    }
};
class MyQwtPlotCurve: public QwtPlotCurve
{
public:
    MyQwtPlotCurve( QVector<QwtPoint3D> *_data, QVector<QwtPointPolar> *_datap )
    {
     QLOG_DEBUG() << "MyQwtPlotCurve:MyQwtPlotCurve" ;
          
     data = _data;
     datap = _datap;
     // Init symbols
     setSymbol(new MyQwtSymbol(Qt::magenta,datap->at(0)));
    }
    void drawSymbols (QPainter *painter,
		      const QwtSymbol &symbol,
		      const QwtScaleMap &xMap,
		      const QwtScaleMap &yMap,
		      const QRectF &canvasRect,
		      int   from,
		      int   to)	const 
    {
    //QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols";
    //QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:from" << from << " to:" << to;
    QBrush brush = Qt::magenta;
    QwtPointMapper mapper;
    mapper.setFlag( QwtPointMapper::RoundPoints, 
       QwtPainter::roundingAlignment( painter ) );
    mapper.setFlag( QwtPointMapper::WeedOutPoints, 
       testPaintAttribute( QwtPlotCurve::FilterPoints ) );
    mapper.setBoundingRect( canvasRect );
  
    
    //QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:data->size()=" << data->size();
    // Reconstruct symbols
    // Calculate average distance
    double distance = 0;
     for ( int j = 0 ; j < data->size() ; j++)
        distance = distance + fabs(data->at(j).x());
     distance/=data->size();
     //QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:Average distance =" 
      //            << distance;
     
    for (int i = 0 ; i < data->size(); i++) {
     QVector<double> xvect, yvect,fxvect, fyvect;
     xvect.push_back( data->at(i).x());
     yvect.push_back( data->at(i).y());     
     QwtPointArrayData *dataarr = new QwtPointArrayData(xvect,yvect);
     const QPolygonF points = mapper.toPointsF( xMap, yMap, dataarr, 0, 0);
     
     QLOG_DEBUG() << "distance* radius " << QString::number(datap->at(i).radius() * 
              distance);
     if (datap->at(i).radius() != 0)
       fxvect.push_back( data->at(i).x() + datap->at(i).radius() * distance / 5 + distance / 10);
     else
       fxvect.push_back( data->at(i).x() );
     fyvect.push_back( data->at(i).y());
     QwtPointArrayData *dataarrp = new QwtPointArrayData(fxvect,fyvect);
     const QPolygonF fpoints = mapper.toPointsF( xMap, yMap, dataarrp, 0, 0);
     if ( fpoints.size() > 0 ) {
	QwtPointPolar dataptrans(datap->at(i).azimuth(),
				fpoints.at(0).x() - points.at(0).x());
	QLOG_DEBUG() << "MyQwtPlotCurve::dataptrans radius()=" 
		     << dataptrans.radius();
	QLOG_DEBUG() << "MyQwtPlotCurve::dataptrans::azimuth()=" 
                     << dataptrans.azimuth();
	MyQwtSymbol mysymbol(brush,dataptrans);
	if ( points.size() > 0 )
	  mysymbol.drawSymbols( painter, points );
       }
       
       QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:data->at(" << i << ").x=" 
                   << data->at(i).x();
       QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:data->at(" << i << ").y=" 
                   << data->at(i).y();
       QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:datap->at(" << i << ").radius=" 
                   << datap->at(i).radius();
       QLOG_DEBUG() << "MyQwtPlotCurve::drawSymbols:datap->at(" << i << ").azimuth=" 
                   << datap->at(i).azimuth();
     
       QLOG_DEBUG() << "MyQwtPlotCurve::xvect=" << xvect.at(0);
       QLOG_DEBUG() << "MyQwtPlotCurve::yvect=" << yvect.at(0);
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::points.x()=" << points.at(0).x();
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::points.y()=" << points.at(0).y();
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::points.size()=" << points.size();
       QLOG_DEBUG() << "MyQwtPlotCurve::fxvect=" << fxvect.at(0);
       QLOG_DEBUG() << "MyQwtPlotCurve::fyvect=" << fyvect.at(0);       
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::fpoints.x()=" << fpoints.at(0).x();
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::fpoints.y()=" << fpoints.at(0).y();
       QLOG_DEBUG() << "MyQwtPlotCurve::QPolygonF::fpoints.size()=" << fpoints.size(); 
    }
   }
   ~MyQwtPlotCurve() {
      delete data;
      delete datap;
   }
private:
    QVector<QwtPoint3D>    *data;
    QVector<QwtPointPolar> *datap;
    
};
class MyQwtPlotZoomer: public QwtPlotZoomer
{
public:
    MyQwtPlotZoomer( QWidget *canvas ):QwtPlotZoomer( canvas )
    {
        setTrackerMode( AlwaysOn );
    }

    virtual QwtText trackerTextF( const QPointF &pos ) const
    {
        QColor bg( Qt::white );
        bg.setAlpha( 160 );
        QwtText text("(" + QString::number(pos.x(),'E',2) + 
		"," + QString::number(pos.y(),'E',2) + ") ");
        text.setBackgroundBrush( QBrush( bg ) );
        return text;
    }
};

class MyQwtScaleDraw: public QwtScaleDraw
{
public:
    MyQwtScaleDraw():QwtScaleDraw()
    {
    }

    virtual QwtText label( double value ) const
    {
        QwtText text(QString::number(value,'E',3));
        return text;
    }
};

class MyQwtMatrixRasterData: public QwtMatrixRasterData
{
public:
    MyQwtMatrixRasterData():QwtMatrixRasterData()
    {      
    }

};
class MyQwtLinearColorMap: public QwtLinearColorMap
{
public:
    MyQwtLinearColorMap(QVector<QColor> *colors):QwtLinearColorMap()
    {
       setColorInterval(colors->at(0), colors->at(1));
       setMode( QwtLinearColorMap::ScaledColors);
    }
};
PlotRaster::PlotRaster( QWidget *parent , QVector<QwtPoint3D> *_data, 
			const int &_num_col, QString _title, QString _xtitle, QString _ytitle,
			QVector<QColor> *_colors):
    QwtPlot( parent )
{

    QLOG_DEBUG() << "PlotRaster::PlotRaster" ;
   
    data = _data;
    title = _title;
    xtitle = _xtitle;
    ytitle = _ytitle;
    num_col = _num_col;
    colors = _colors;
    this->setTitle( title );
    d_spectrogram = new QwtPlotSpectrogram();
    d_spectrogram->setRenderThreadCount( 0 ); // use system specific thread count
    MyQwtLinearColorMap *colormap =  new MyQwtLinearColorMap(colors);
    d_spectrogram->setColorMap( colormap );
   
    d_spectrogram->setCachePolicy( QwtPlotRasterItem::PaintCache );
    MyQwtMatrixRasterData *rasterdata = new MyQwtMatrixRasterData();
    QVector<double> zvalues;
    for (int i = 0 ; i < data->size(); i++)
       zvalues.push_back(data->at(i).z());
    
    rasterdata->setValueMatrix(zvalues,num_col);
    double minx=1e308,maxx=-1e308;
    double miny=1e308,maxy=-1e308;
    double minz=1e308,maxz=-1e308;
    for (int i = 0 ; i < data->size(); i++) {
       QwtPoint3D point = data->at(i);
       double tmpx = point.x();
       double tmpy = point.y();
       double tmpz = point.z();
       if (tmpx > maxx) maxx = tmpx;
       if (tmpx < minx) minx = tmpx;
       if (tmpy > maxy) maxy = tmpy;
       if (tmpy < miny) miny = tmpy;
       if (tmpz == -100) continue;
       if (tmpz > maxz) maxz = tmpz;
       if (tmpz < minz) minz = tmpz;
    }
    if (data->size() == 0) {
       minx = maxx = miny = maxy = minz = maxz = 0;
    }
    QLOG_INFO() << "PlotRaster: Number of columns:" << num_col;
    QLOG_INFO() << "PlotRaster: Data size:" << data->size();
    QLOG_INFO() << "PlotRaster: Num columns:" << rasterdata->numColumns();
    QLOG_INFO() << "PlotRaster: Num rows:" << rasterdata->numRows();
    QLOG_INFO() << "Min x:" << QString::number(minx);
    QLOG_INFO() << "Max x:" << QString::number(maxx);
    QLOG_INFO() << "Min y:" << QString::number(miny);
    QLOG_INFO() << "Max y:" << QString::number(maxy);
    QLOG_INFO() << "Min field:" << QString::number(minz);
    QLOG_INFO() << "Max field:" << QString::number(maxz);
    QwtInterval rangex = QwtInterval(minx,maxx);
    QwtInterval rangey = QwtInterval(miny,maxy);
    QwtInterval rangez;
    if (minz == maxz) 
      rangez = QwtInterval(minz*0.5,maxz*1.5);
    else
      rangez = QwtInterval(minz,maxz);
    rasterdata->setInterval( Qt::XAxis, rangex );
    rasterdata->setInterval( Qt::YAxis, rangey );
    rasterdata->setInterval( Qt::ZAxis, rangez );
    d_spectrogram->setData( rasterdata );
    d_spectrogram->attach( this );

    // Axis range
    this->setAxisScale(QwtPlot::xBottom,minx,maxx,0);
    this->setAxisScale(QwtPlot::yLeft,miny,maxy,0);
    this->setAxisTitle(QwtPlot::xBottom, xtitle);
    this->setAxisTitle(QwtPlot::yLeft, ytitle);
    // A color bar on the right axis
    QwtScaleWidget *rightAxis = this->axisWidget( QwtPlot::yRight );
    rightAxis->setColorBarEnabled( true );
    MyQwtLinearColorMap *colormap2 =  new MyQwtLinearColorMap(colors);  
    rightAxis->setColorMap( rangez, colormap2 );
    rightAxis->setSpacing(2);    
    rightAxis->setScaleDraw(new MyQwtScaleDraw());
    setAxisScale( QwtPlot::yRight, minz, maxz, 0 );
    
    enableAxis( QwtPlot::yRight );
    plotLayout()->setAlignCanvasToScales( true );
    replot();

    // LeftButton for the zooming
    // MidButton for the panning
    // RightButton: zoom out by 1
    // Ctrl+RighButton: zoom out to full size

    QwtPlotZoomer* zoomer = new MyQwtPlotZoomer( canvas() );
    zoomer->setMousePattern( QwtEventPattern::MouseSelect2,
        Qt::RightButton, Qt::ControlModifier );
    zoomer->setMousePattern( QwtEventPattern::MouseSelect3,
        Qt::RightButton );

    QwtPlotPanner *panner = new QwtPlotPanner( canvas() );
    panner->setAxisEnabled( QwtPlot::yRight, false );
    panner->setMouseButton( Qt::MidButton );

    // Avoid jumping when labels with more/less digits
    // appear/disappear when scrolling vertically

    const QFontMetrics fm( axisWidget( QwtPlot::yLeft )->font() );
    QwtScaleDraw *sd = axisScaleDraw( QwtPlot::yLeft );
    sd->setMinimumExtent( fm.width( "100.00" ) );

    const QColor c( Qt::darkBlue );
    zoomer->setRubberBandPen( c );
    zoomer->setTrackerPen( c );
}
PlotRaster::~PlotRaster() {
    QLOG_DEBUG() << "PlotRaster::~PlotRaster:: Deleting Plot";
    delete d_spectrogram;
}
PlotVector::PlotVector( QWidget *parent , QVector<QwtPoint3D> *_data, 
			QVector<QwtPointPolar> *_datap, int _num_col, QString _title,
			QString _xtitle, QString _ytitle):
    QwtPlot( parent )
{
    QLOG_DEBUG() << "PlotVector::PlotVector" ;

    data = _data;
    datap = _datap;
    title = _title;
    xtitle = _xtitle;
    ytitle = _ytitle;
    num_col = _num_col;
    d_curve = new MyQwtPlotCurve(data, datap);
    d_curve->setRenderThreadCount( 0 ); // use system specific thread count
    d_curve->setStyle(QwtPlotCurve::NoCurve);
    QPolygonF points;
    for (int i = 0 ; i < data->size(); i++) {
        points << QPointF(data->at(i).x(),data->at(i).y());
    }
    d_curve->setSamples( points );
    this->setTitle( title );
    d_curve->attach( this );
    double minx=1e308,maxx=-1e308;
    double miny=1e308,maxy=-1e308;
    double minz=1e308,maxz=-1e308;
    for (int i = 0 ; i < data->size(); i++) {
       QwtPoint3D point = data->at(i);
       double tmpx = point.x();
       double tmpy = point.y();
       double tmpz = point.z();
       if (tmpx > maxx) maxx = tmpx;
       if (tmpx < minx) minx = tmpx;
       if (tmpy > maxy) maxy = tmpy;
       if (tmpy < miny) miny = tmpy;
       if (tmpz == -100) continue;
       if (tmpz > maxz) maxz = tmpz;
       if (tmpz < minz) minz = tmpz;
    }
    if (data->size() == 0) {
       minx = maxx = miny = maxy = minz = maxz = 0;
    }
     double xabsmax;
    if (fabs(minx) > fabs(maxx))
       xabsmax = fabs(maxx);
    else
       xabsmax = fabs(minx);
    double yabsmax;
    if (fabs(miny) > fabs(maxy))
       yabsmax = fabs(maxy);
    else
       yabsmax = fabs(miny);
    QLOG_DEBUG() << " xabsmax = " << QString::number(xabsmax,'g',8);
    QLOG_DEBUG() << " yabsmax = " << QString::number(yabsmax,'g',8);
    
    double xmin = minx - xabsmax/1.5;
    double xmax = maxx + xabsmax/1.5;
    double ymin = miny - yabsmax/1.5;
    double ymax = maxy + yabsmax/1.5;
   
    this->setAxisScale( QwtPlot::xBottom, xmin, xmax, 0 );
    this->setAxisScale( QwtPlot::yLeft, ymin, ymax, 0 );
    this->setAxisTitle(QwtPlot::xBottom, xtitle);
    this->setAxisTitle(QwtPlot::yLeft, ytitle);
    plotLayout()->setAlignCanvasToScales( true );
    replot();
    // LeftButton for the zooming
    // MidButton for the panning
    // RightButton: zoom out by 1
    // Ctrl+RighButton: zoom out to full size

    QwtPlotZoomer* zoomer = new MyQwtPlotZoomer( this->canvas() );
    zoomer->setZoomBase(QRectF(QPointF(xmin,ymax),QPointF(xmax,ymin))); 
    zoomer->setMousePattern( QwtEventPattern::MouseSelect2,
        Qt::RightButton, Qt::ControlModifier );
    zoomer->setMousePattern( QwtEventPattern::MouseSelect3,
        Qt::RightButton );

    QwtPlotPanner *panner = new QwtPlotPanner( canvas() );
    panner->setAxisEnabled( QwtPlot::yRight, false );
    panner->setMouseButton( Qt::MidButton );

    // Avoid jumping when labels with more/less digits
    // appear/disappear when scrolling vertically

    const QFontMetrics fm( this->axisWidget( QwtPlot::yLeft )->font() );
    QwtScaleDraw *sd = axisScaleDraw( QwtPlot::yLeft );
    sd->setMinimumExtent( fm.width( "100.00" ) );

    const QColor c( Qt::darkBlue );
    zoomer->setRubberBandPen( c );
    zoomer->setTrackerPen( c );
}
PlotVector::~PlotVector() {
    QLOG_DEBUG() << "PlotVector::~PlotVector:: Deleting Plot";
    delete d_curve;
}
MyColor::MyColor(TripleVector  *tripledata , unsigned size)

{
        QLOG_DEBUG() << "MyColor::MyColor> Create MyColor";
        tripledata_ = tripledata;
        size_ = size;
        //reset(size_);
}
MyColor::~MyColor() {
}
void 
MyColor::reset(unsigned size)
{
	colors_ = ColorVector(size);
	RGBA elem;

	double dsize = size;
	
	for (unsigned int i=0; i!=size; ++i)
	{
		elem.r = i / dsize;
		elem.g = i / dsize / 4;
		elem.b = 1 - i/dsize;
		elem.a = 1.0;
		colors_[i] = elem;
	}
}

/**
	Assigns a new ColorVector (Also overwrites the constructors size argument)
*/
void 
MyColor::setColorVector(ColorVector const& cv)
{
    colors_ = cv;
}

void 
MyColor::setAlpha(double a)
{
	if (a<0 || a>1)
		return;
	
	RGBA elem;

	for (unsigned int i=0; i!=colors_.size(); ++i)
	{
		elem = colors_[i];
		elem.a = a;
		colors_[i] = elem;
	}	
}	

RGBA
MyColor::rgba(double x, double y, double z) const
{
    int index = 0;
         //printf("tripledata size=%d data size=%d\n",tripledata_->size(),size_);
        if (tripledata_->size() != size_)
            return colors_[0];;
        for (int i = 0 ; i < tripledata_->size(); i++) {
       Triple elem = tripledata_->at(i);
       if ( elem.x == x &&  elem.y == y && elem.z == z ){
             index = i;
             break;
           }
        }
        if (index < 0)
       index = 0;
        if ((unsigned int)index > colors_.size() - 1)
       index = (int)(colors_.size() - 1);
       // printf("x=%E,y=%E,z=%E index=%d r=%f,g=%f,b=%f,size=%d\n",x,y,z,index,colors_[index].r,colors_[index].g,colors_[index].b, tripledata_->size());
    return colors_[index];
}
void 
MyColor::update(const Qwt3D::Plot3D &val)
{
}

Plot::Plot( QVector<QwtPoint3D> *_data, QVector<double> *_colormap, QString _title, int _param ) {

    data = _data;  
    colormap = _colormap;
    title = _title;
    param = _param;

    tripleVect = new TripleVector();
    cellVect = new CellVector();
    colorvect = new ColorVector();
   
    //
    // Fill data
    //
    for (int j = 0; j < data->size(); j++) {
          QLOG_DEBUG() << " data size = " << data->size() << " j = " << j;
	  Triple triple;
          if (colormap->at(j) == 1) continue;
	  triple.x = data->at(j).x();
	  triple.y = data->at(j).y();
	  triple.z = data->at(j).z();
	  tripleVect->push_back(triple);
	  // quadrangles
	 if ( title == "Poynting" ) {
	  int count = param - 1;
	  if (j%(count+1) == 0 ) {
		  Cell c(4);
	  	  for (int k=0; k < count+1; ++k) {
		     if (k+j < data->size())
	    		c[0] = k+j;
		     if (k%count+1+j < data->size())
	    		c[1] = k%count+1+j;
		     if (count+k%count+2+j < data->size())
	    		c[2] = count+k%count+2+j;
		     if (count+k+j+1 < data->size())
	    		c[3] = count+k+j+1;
	    	     cellVect->push_back(c);
	  	  }
	  }
	}
        else if ( title == "Dipoles" ) {
	   Cell c;
           cellVect->push_back(c);
        }
    }
    (this)->createDataset(*tripleVect,*cellVect);
    
   if ( title == "Dipoles" ) {
     this->setPlotStyle(Qwt3D::POINTS);
   }
   if ( title == "Poynting" )
    this->setPlotStyle(Qwt3D::FILLED);

    //
    // Fill color
    //
    
    if ( colormap != NULL ) {
       double maxcolormap = -1.0E-300;
       double mincolormap = 1.0E300;
      //
      // find colormap maximum,minimum value
      //
      for (int i = 0 ; i < data->size() ; i++ ) {
       if ( colormap->at(i) > maxcolormap)
          maxcolormap = colormap->at(i);
       if ( colormap->at(i) < mincolormap)
          mincolormap = colormap->at(i);
      }
      QLOG_DEBUG() << "Plot::Plot> maxcolormap=" << QString::number(maxcolormap);
      for (int i = 0 ; i < data->size() ; i++ ) {
         if (colormap->at(i) == 1) continue;
         if (title == "Dipoles")
           colorvect->push_back(RGBA((colormap->at(i))/(maxcolormap),0.,0.,0.8));
         else if (title == "Poynting") 
           colorvect->push_back(RGBA(colormap->at(i)/(maxcolormap),0.,0.,0.8));
      }
      color = new MyColor(tripleVect,tripleVect->size());
      color->setColorVector(*colorvect);
      this->setDataColor(*color);
      this->updateNormals();
      this->legend()->setLimits(mincolormap, maxcolormap);
      this->legend()->setAutoScale(true);
      this->showColorLegend(true);
    }
    this->setTitle(title);
    for (unsigned i = 0 ; i != this->coordinates()->axes.size(); ++i) {
      this->coordinates()->axes[i].setMajors(7);
      this->coordinates()->axes[i].setMinors(4);
    }
    this->coordinates()->axes[X1].setLabelString("x");
    this->coordinates()->axes[Y1].setLabelString("y");
    this->coordinates()->axes[Z1].setLabelString("z");
    //this->setCoordinateStyle(BOX);
    this->updateData();
    this->updateGL();
  }
Plot::~Plot() {
    QLOG_DEBUG() << "Plot::~Plot:: Deleting Plot";
    delete data;
    delete colormap;
    delete cellVect;
    delete tripleVect;
    delete colorvect;
}
