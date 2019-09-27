#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_symbol.h>
#include <qwt_plot_spectrogram.h>
#include <qprinter.h>
#include <qprintdialog.h>
#include <qwt_color_map.h>
#include <qwt_scale_widget.h>
#include <qwt_scale_draw.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_layout.h>
#include <qwt_plot_renderer.h>
#include <qwt_text_engine.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_point_mapper.h>
#include <qwt_painter.h>
#include <qwt_series_data.h>
#include <qwt_point_data.h>
#include <QBuffer>
#include <QBoxLayout>
#include <qwt3d_plot3d.h>
#include <qwt3d_meshplot.h>
#include <qwt3d_color.h>
#include <qwt3d_function.h>
#include "QsLog.h"

using namespace Qwt3D;

class MyColor : public Color
{
public:
        explicit MyColor(TripleVector  *tripledata, unsigned size = 100);
        ~MyColor();
	Qwt3D::RGBA rgba(double x, double y, double z) const;
	void setColorVector(Qwt3D::ColorVector const& cv);
	void reset(unsigned size=100);
	void setAlpha(double a);
        Color* clone() const {return new MyColor(*this);}
	Qwt3D::ColorVector& createVector(Qwt3D::ColorVector& vec) {vec = colors_; return vec;}
        void update(const Plot3D& val);

protected:
	Qwt3D::ColorVector colors_;
        TripleVector *tripledata_;
        unsigned int size_;
};

class Plot : public MeshPlot
{
  Q_OBJECT

  public:
      Plot( QVector<QwtPoint3D> *_data, QVector<double> *colormap, QString _title, int _param);
      ~Plot();
      QString getTitle() { return title;};
   private:
      MyColor *color;
      QVector<double> *colormap;
      QVector<QwtPoint3D> *data;
      TripleVector *tripleVect;
      CellVector *cellVect;
      ColorVector *colorvect;
      QString title;
      int param;
};

class PlotRaster: public QwtPlot
{
    Q_OBJECT

public:
    PlotRaster( QWidget * = NULL, QVector<QwtPoint3D> *_data = NULL, 
                const int &_num_col = 0, QString _title = 0, QString _xtitle = 0,
		QString _ytitle = 0, QVector<QColor> *_colors = NULL);
    ~PlotRaster();
    
private:
    QwtPlotSpectrogram *d_spectrogram;
    QVector<QwtPoint3D> *data;
    QVector<QColor> *colors;
    int num_col;
    QString title;
    QString xtitle;
    QString ytitle;
};

class PlotVector: public QwtPlot
{
    Q_OBJECT

public:
    PlotVector( QWidget * = NULL, QVector<QwtPoint3D> *_data = NULL, 
                QVector<QwtPointPolar> *_datap = NULL,
                int _num_col = 0, QString _title = 0, QString _xtitle = 0,
		QString _ytitle = 0);
    ~PlotVector();

private:
    QwtPlotCurve *d_curve;
    QVector<QwtPoint3D> *data;
    QVector<QwtPointPolar> *datap;
    int num_col;
    QString title;
    QString xtitle;
    QString ytitle;
};

