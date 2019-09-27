#include "cdmRun.h"

Run::Run(QString _runname)
{
  runname = _runname;

  xc = NULL;
  yc = NULL;
  zc = NULL;
  kxy = NULL;
  xy = NULL;
  xcwf = NULL;
  ycwf = NULL;
  zcwf = NULL;
  incidentfield = NULL;
  localfield = NULL;
  macroscopicfield = NULL;
  thetafield = NULL;
  phifield = NULL;
  poyntingfield = NULL;  
  forcex = NULL;
  forcey = NULL;
  forcez = NULL;
  forcexmulti = NULL;
  forceymulti = NULL;
  forcezmulti = NULL;
  torquex = NULL;
  torquey = NULL;
  torquez = NULL;
  torquexmulti = NULL;
  torqueymulti = NULL;
  torquezmulti = NULL;
  incidentfieldx = NULL;
  localfieldx = NULL;
  macroscopicfieldx = NULL;
  incidentfieldy = NULL;
  localfieldy = NULL;
  macroscopicfieldy = NULL;
  incidentfieldz = NULL;
  localfieldz = NULL;
  macroscopicfieldz = NULL;
  polarisafield = NULL;
  epsilonfield = NULL;
  eimagex = NULL;
  eimagey = NULL;
  eimagez = NULL;
  efourierx = NULL;
  efouriery = NULL;
  efourierz = NULL;
  efourierincx = NULL;
  efourierincy = NULL;
  efourierincz = NULL;
  eimageincx = NULL;
  eimageincy = NULL;
  eimageincz = NULL;
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  FF = NULL;
  FF0 = NULL;
  FFloc = NULL;
  xr = NULL;
  xi = NULL;
  wrk = NULL;
  FFTTENSORxx = NULL;
  FFTTENSORxy = NULL;
  FFTTENSORxz = NULL;
  FFTTENSORyy = NULL;
  FFTTENSORyz = NULL;
  FFTTENSORzz = NULL;
  vectx = NULL;
  vecty = NULL;
  vectz = NULL;
  Ediffkzpos = NULL;
  Ediffkzneg = NULL;
  Tabdip  = NULL;
  Tabmulti = NULL;
}
Run::~Run() {
  cleanVectorsMemory();
}
int
Run::checkAvailableMemorySize() {

  QLOG_INFO() << "Run debut mem ";
#ifdef OS
#if OS == LINUX
  //    long long pages = sysconf(_SC_AVPHYS_PAGES);
    long long pages = sysconf(_SC_PHYS_PAGES);
    long long page_size = sysconf(_SC_PAGE_SIZE);
    int mem_size = (int)((pages * page_size)/1000000L);
    return mem_size;
#elif OS == WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (unsigned long long) (status.ullAvailPhys);
#endif
#endif
   return 0;
}
int
Run::allocateVectorsMemory(int nmax, int ntheta, int nphi, int nfft2d, int obj_num) {
  
  QLOG_INFO() << "Run::allocateVectorsMemory> "; 

  int mem_used = 0;

  QLOG_INFO() << "Run::allocateVectorsMemory::Memory allocation dcmplx:" << (int)sizeof(dcmplx);
  QLOG_INFO() << "Run::allocateVectorsMemory::Memory allocation double:" << (int)sizeof(double);
  xc = (double*) malloc(sizeof(double)*nmax);
  if (xc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  yc = (double*) malloc(sizeof(double)*nmax);
  if (yc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(yc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  zc = (double*) malloc(sizeof(double)*nmax);
  if (zc == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(zc,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  xcwf = (double*) malloc(sizeof(double)*nmax);
  if (xcwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xcwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  ycwf = (double*) malloc(sizeof(double)*nmax);
  if (ycwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(ycwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  zcwf = (double*) malloc(sizeof(double)*nmax);
  if (zcwf == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(zcwf,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  incidentfield = (double*) malloc(sizeof(double)*nmax);
  if (incidentfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  localfield = (double*) malloc(sizeof(double)*nmax);
  if (localfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  macroscopicfield = (double*) malloc(sizeof(double)*nmax);
  if (macroscopicfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  } 
  memset(macroscopicfield,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  thetafield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (thetafield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(thetafield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  phifield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (phifield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(phifield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  poyntingfield = (double*) malloc(sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  if (poyntingfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(poyntingfield,0,sizeof(double)*max((ntheta+1)*nphi,nfft2d*nfft2d));
  mem_used+=max((ntheta+1)*nphi,nfft2d*nfft2d)*sizeof(double) / 1000000L;
  forcex = (double*) malloc(sizeof(double)*nmax);
  if (forcex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcex,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  forcey = (double*) malloc(sizeof(double)*nmax);
  if (forcey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcey,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  forcez = (double*) malloc(sizeof(double)*nmax);
  if (forcez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcez,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  
  forcexmulti = (double*) malloc(sizeof(double)*obj_num);
  if (forcexmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcexmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  forceymulti = (double*) malloc(sizeof(double)*obj_num);
  if (forceymulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forceymulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  forcezmulti = (double*) malloc(sizeof(double)*obj_num);
  if (forcezmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(forcezmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torquex = (double*) malloc(sizeof(double)*nmax);
  if (torquex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }

  memset(torquex,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  torquey = (double*) malloc(sizeof(double)*nmax);
  if (torquey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquey,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;
  torquez = (double*) malloc(sizeof(double)*nmax);
  if (torquez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquez,0,sizeof(double)*nmax);
  mem_used+=nmax*sizeof(double) / 1000000L;

  torquexmulti = (double*) malloc(sizeof(double)*obj_num);
  if (torquexmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquexmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torqueymulti = (double*) malloc(sizeof(double)*obj_num);
  if (torqueymulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torqueymulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;
  torquezmulti = (double*) malloc(sizeof(double)*obj_num);
  if (torquezmulti == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(torquezmulti,0,sizeof(double)*obj_num);
  mem_used+=obj_num*sizeof(double) / 1000000L;



  macroscopicfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  macroscopicfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  macroscopicfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (macroscopicfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(macroscopicfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  localfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (localfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(localfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldx = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldx,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldy = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldy,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  incidentfieldz = (dcmplx*) malloc(sizeof(dcmplx)*nmax);
  if (incidentfieldz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(incidentfieldz,0,sizeof(dcmplx)*nmax);
  mem_used+= sizeof(dcmplx)*nmax/ 1000000L;
  polarisafield = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*3);
  if (polarisafield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(polarisafield,0,sizeof(dcmplx)*nmax*3*3);
  mem_used+= sizeof(dcmplx)*nmax*3*3/ 1000000L;
  epsilonfield = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*3);
  if (epsilonfield == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(epsilonfield,0,sizeof(dcmplx)*nmax*3*3);
  mem_used+= sizeof(dcmplx)*nmax*3*3/ 1000000L;
  kxy = (double*) malloc(sizeof(double)*nfft2d);
  if (kxy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(kxy,0,sizeof(double)*nfft2d);
  mem_used+= sizeof(double)*nfft2d/ 1000000L;
  xy = (double*) malloc(sizeof(double)*nfft2d);
  if (xy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xy,0,sizeof(double)*nfft2d);
  mem_used+= sizeof(double)*nfft2d/ 1000000L;
  eimagex = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagex == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagex,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimagey = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagey == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagey,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimagez = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimagez == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimagez,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efouriery = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efouriery == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efouriery,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincy = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincy,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  efourierincz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (efourierincz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(efourierincz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincx = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincx == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincx,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincy = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincy == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincy,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
  eimageincz = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d);
  if (eimageincz == NULL) { 
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(eimageincz,0,sizeof(dcmplx)*nfft2d*nfft2d);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d/ 1000000L;
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  FF = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FF == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FF,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  FF0 = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FF0 == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FF0,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  FFloc = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (FFloc == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFloc,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  xr = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (xr == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xr,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  xi = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3);
  if (xi == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(xi,0,sizeof(dcmplx)*nmax*3);
  mem_used+= sizeof(dcmplx)*nmax*3/ 1000000L;
  wrk = (dcmplx*) malloc(sizeof(dcmplx)*nmax*3*12);
  if (wrk == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(wrk,0,sizeof(dcmplx)*nmax*3*12);
  mem_used+= sizeof(dcmplx)*nmax*3*12/ 1000000L;

  FFTTENSORxx = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORxx == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORxx,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  FFTTENSORxy = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORxy == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORxy,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  FFTTENSORxz = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORxz == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORxz,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  FFTTENSORyy = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORyy == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORyy,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  FFTTENSORyz = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORyz == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORyz,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  FFTTENSORzz = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (FFTTENSORzz == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(FFTTENSORzz,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  
  vectx = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (vectx == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(vectx,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  vecty = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (vecty == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(vecty,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  vectz = (dcmplx*) malloc(sizeof(dcmplx)*nmax*8);
  if (vectz == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(vectz,0,sizeof(dcmplx)*nmax*8);
  mem_used+= sizeof(dcmplx)*nmax*8/ 1000000L;
  Ediffkzpos = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d*3);
  if (Ediffkzpos == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Ediffkzpos,0,sizeof(dcmplx)*nfft2d*nfft2d*3);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d*3/ 1000000L;
  
  Ediffkzneg = (dcmplx*) malloc(sizeof(dcmplx)*nfft2d*nfft2d*3);
  if (Ediffkzneg == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  
  memset(Ediffkzneg,0,sizeof(dcmplx)*nfft2d*nfft2d*3);
  mem_used+= sizeof(dcmplx)*nfft2d*nfft2d*3/ 1000000L;

  Tabdip = (int*) malloc(sizeof(int)*nmax);
  if (Tabdip == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }

  memset(Tabdip,0,sizeof(int)*nmax);
  mem_used+= sizeof(int)*nmax/ 1000000L;
  Tabmulti = (int*) malloc(sizeof(int)*nmax);
  if (Tabmulti == NULL) {
    QLOG_FATAL() << "Run::allocateVectorsMemory::Memory allocation failed";
  return -1;
  }
  memset(Tabmulti,0,sizeof(int)*nmax);
  mem_used+= sizeof(int)*nmax/ 1000000L;
  return mem_used;
}
void 
Run::cleanVectorsMemory() {
  if (xc) {free(xc); xc = NULL;}
  if (yc) {free(yc); yc = NULL;}
  if (zc) {free(zc); zc = NULL;}
  if (xcwf) {free(xcwf); xcwf = NULL;}
  if (ycwf) {free(ycwf); ycwf = NULL;}
  if (zcwf) {free(zcwf); zcwf = NULL;}
  if (xy) {free(xy); xy = NULL;}
  if (kxy) {free(kxy); kxy = NULL;}
  if (incidentfield) {free(incidentfield); incidentfield = NULL;}
  if (localfield) {free(localfield); localfield = NULL;}
  if (macroscopicfield) {free(macroscopicfield); macroscopicfield = NULL;}
  if (thetafield) {free(thetafield); thetafield = NULL;}
  if (phifield) {free(phifield); phifield = NULL;}
  if (poyntingfield) {free(poyntingfield); poyntingfield = NULL;}
  if (forcex) {free(forcex); forcex = NULL;}
  if (forcey) {free(forcey); forcey = NULL;}
  if (forcez) {free(forcez); forcez = NULL;}
  if (forcexmulti) {free(forcexmulti); forcexmulti = NULL;}
  if (forceymulti) {free(forceymulti); forceymulti = NULL;}
  if (forcezmulti) {free(forcezmulti); forcezmulti = NULL;}
  if (torquex) {free(torquex); torquex = NULL;}
  if (torquey) {free(torquey); torquey = NULL;}
  if (torquez) {free(torquez); torquez = NULL;}
  if (torquexmulti) {free(torquexmulti); torquexmulti = NULL;}
  if (torqueymulti) {free(torqueymulti); torqueymulti = NULL;}
  if (torquezmulti) {free(torquezmulti); torquezmulti = NULL;}
  if (incidentfieldx) {free(incidentfieldx); incidentfieldx = NULL;}
  if (localfieldx) {free(localfieldx); localfieldx = NULL;}
  if (macroscopicfieldx) {free(macroscopicfieldx); macroscopicfieldx = NULL;}
  if (incidentfieldy) {free(incidentfieldy); incidentfieldy = NULL;}
  if (localfieldy) {free(localfieldy); localfieldy = NULL;}
  if (macroscopicfieldy) {free(macroscopicfieldy); macroscopicfieldy = NULL;}
  if (incidentfieldz) {free(incidentfieldz); incidentfieldz = NULL;}
  if (localfieldz) {free(localfieldz); localfieldz = NULL;}
  if (macroscopicfieldz) {free(macroscopicfieldz); macroscopicfieldz = NULL;}
  if (polarisafield) {free(polarisafield); polarisafield = NULL;}
  if (epsilonfield) {free(epsilonfield); epsilonfield = NULL;}
  if (eimagex) {free(eimagex); eimagex = NULL;}
  if (eimagey) {free(eimagey); eimagey = NULL;}
  if (eimagez) {free(eimagez); eimagez = NULL;}
  if (efourierx) {free(efourierx); efourierx = NULL;}
  if (efouriery) {free(efouriery); efouriery = NULL;}
  if (efourierz) {free(efourierz); efourierz = NULL;}
  if (efourierincx) {free(efourierincx); efourierincx = NULL;}
  if (efourierincy) {free(efourierincy); efourierincy = NULL;}
  if (efourierincz) {free(efourierincz); efourierincz = NULL;}
  if (eimageincx) {free(eimageincx); eimageincx = NULL;}
  if (eimageincy) {free(eimageincy); eimageincy = NULL;}
  if (eimageincz) {free(eimageincz); eimageincz = NULL;}
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
  if (FF) {free(FF); FF = NULL;}
  if (FF0) {free(FF0); FF0 = NULL;}
  if (FFloc) {free(FFloc); FFloc = NULL;}
  if (xr) {free(xr); xr = NULL;}
  if (xi) {free(xi);  xi = NULL;}
  if (wrk) {free(wrk); wrk = NULL;}
  if (FFTTENSORxx) {free(FFTTENSORxx); FFTTENSORxx = NULL;}
  if (FFTTENSORxy) {free(FFTTENSORxy); FFTTENSORxy = NULL;}
  if (FFTTENSORxz) {free(FFTTENSORxz); FFTTENSORxz = NULL;}
  if (FFTTENSORyy) {free(FFTTENSORyy); FFTTENSORyy = NULL;}
  if (FFTTENSORyz) {free(FFTTENSORyz); FFTTENSORyz = NULL;}
  if (FFTTENSORzz) {free(FFTTENSORzz); FFTTENSORzz = NULL;}
  if (vectx) {free(vectx); vectx = NULL;}
  if (vecty) {free(vecty); vecty = NULL;}
  if (vectz) {free(vectz); vectz = NULL;}
  if (Ediffkzpos) {free(Ediffkzpos); Ediffkzpos = NULL;}
  if (Ediffkzneg) {free(Ediffkzneg); Ediffkzneg = NULL;}
  if (Tabdip) {free(Tabdip); Tabdip = NULL;}
  if (Tabmulti) {free(Tabmulti); Tabmulti = NULL;}
}
double* 
Run::getIncidentField(){
  return incidentfield;
}
double* 
Run::getLocalField(){
  return localfield;
}
double* 
Run::getMacroscopicField(){
  return macroscopicfield;
}
double* 
Run::getXc(){
  return xc;
}
double* 
Run::getYc(){
  return yc;
}
double* 
Run::getZc(){
  return zc;
}
double* 
Run::getXcWF(){
  return xcwf;
}
double* 
Run::getYcWF(){
  return ycwf;
}
double* 
Run::getZcWF(){
  return zcwf;
}
double* 
Run::getThetaField(){
  return thetafield;
}
double* 
Run::getPhiField(){
  return phifield;
}
double* 
Run::getPoyntingField(){
  return poyntingfield;
}
double* 
Run::getForceX(){
  return forcex;
}
double* 
Run::getForceY(){
  return forcey;
}
double* 
Run::getForceZ(){
  return forcez;
}
double* 
Run::getForceXMulti(){
  return forcexmulti;
}
double* 
Run::getForceYMulti(){
  return forceymulti;
}
double* 
Run::getForceZMulti(){
  return forcezmulti;
}
double* 
Run::getTorqueX(){
  return torquex;
}
double* 
Run::getTorqueY(){
  return torquey;
}
double* 
Run::getTorqueZ(){
  return torquez;
}
double* 
Run::getTorqueXMulti(){
  return torquexmulti;
}
double* 
Run::getTorqueYMulti(){
  return torqueymulti;
}
double* 
Run::getTorqueZMulti(){
  return torquezmulti;
}
dcmplx* 
Run::getIncidentFieldX(){
  return incidentfieldx;
}
dcmplx* 
Run::getLocalFieldX(){
  return localfieldx;
}
dcmplx* 
Run::getMacroscopicFieldX(){
  return macroscopicfieldx;
}
dcmplx* 
Run::getIncidentFieldY(){
  return incidentfieldy;
}
dcmplx* 
Run::getLocalFieldY(){
  return localfieldy;
}
dcmplx* 
Run::getMacroscopicFieldY(){
  return macroscopicfieldy;
}
dcmplx* 
Run::getIncidentFieldZ(){
  return incidentfieldz;
}
dcmplx* 
Run::getLocalFieldZ(){
  return localfieldz;
}
dcmplx* 
Run::getMacroscopicFieldZ(){
  return macroscopicfieldz;
}
dcmplx*
Run::getPolarisaField() {
  return (dcmplx*)polarisafield;
}
dcmplx*
Run::getEpsilonField() {
  return (dcmplx*)epsilonfield;
}
double*
Run::getXY() {
  return xy;
}
double*
Run::getKXY() {
  return kxy;
}
dcmplx*
Run::getEimageX() {
  return (dcmplx*)eimagex;
}
dcmplx*
Run::getEimageY() {
  return (dcmplx*)eimagey;
}
dcmplx*
Run::getEimageZ() {
  return (dcmplx*)eimagez;
}
dcmplx*
Run::getEfourierX() {
  return (dcmplx*)efourierx;
}
dcmplx*
Run::getEfourierY() {
  return (dcmplx*)efouriery;
}
dcmplx*
Run::getEfourierZ() {
  return (dcmplx*)efourierz;
}
dcmplx*
Run::getEfourierincX() {
  return (dcmplx*)efourierincx;
}
dcmplx*
Run::getEfourierincY() {
  return (dcmplx*)efourierincy;
}
dcmplx*
Run::getEfourierincZ() {
  return (dcmplx*)efourierincz;
}
dcmplx*
Run::getEimageincX() {
  return (dcmplx*)eimageincx;
}
dcmplx*
Run::getEimageincY() {
  return (dcmplx*)eimageincy;
}
dcmplx*
Run::getEimageincZ() {
  return (dcmplx*)eimageincz;
}
dcmplx*
Run::getFF() {
  return (dcmplx*)FF;
}
dcmplx*
Run::getFF0() {
  return (dcmplx*)FF0;
}
dcmplx*
Run::getFFloc() {
  return (dcmplx*)FFloc;
}
dcmplx*
Run::getxr() {
  return (dcmplx*)xr;
}
dcmplx*
Run::getxi() {
  return (dcmplx*)xi;
}
dcmplx*
Run::getwrk() {
  return (dcmplx*)wrk;
}
dcmplx*
Run::getFFTTENSORxx() {
  return (dcmplx*)FFTTENSORxx;
}
dcmplx*
Run::getFFTTENSORxy() {
  return (dcmplx*)FFTTENSORxy;
}
dcmplx*
Run::getFFTTENSORxz() {
  return (dcmplx*)FFTTENSORxz;
}
dcmplx*
Run::getFFTTENSORyy() {
  return (dcmplx*)FFTTENSORyy;
}
dcmplx*
Run::getFFTTENSORyz() {
  return (dcmplx*)FFTTENSORyz;
}dcmplx*
Run::getFFTTENSORzz() {
  return (dcmplx*)FFTTENSORzz;
}
dcmplx*
Run::getvectx() {
  return (dcmplx*)vectx;
}
dcmplx*
Run::getvecty() {
  return (dcmplx*)vecty;
}
dcmplx*
Run::getvectz() {
  return (dcmplx*)vectz;
}
dcmplx*
Run::getEdiffkzpos() {
  return (dcmplx*)Ediffkzpos;
}
dcmplx*
Run::getEdiffkzneg() {
  return (dcmplx*)Ediffkzneg;
}
int*
Run::getTabdip() {
  return (int*)Tabdip;
}
int*
Run::getTabmulti() {
  return (int*)Tabmulti;
}
void 
Run::setName(QString _runname){
  runname = _runname;
}
void 
Run::setObjectSubunits(int _objectsubunits){
  objectsubunits = _objectsubunits;
}
void 
Run::setMeshSubunits(int _meshsubunits){
  meshsubunits = _meshsubunits;
}
void 
Run::setNmaxpp(int _nmaxpp){
  nmaxpp = _nmaxpp;
}
void 
Run::setMeshSize(double _meshsize){
  meshsize = _meshsize;
}
void 
Run::setLambda10n(double _lambda10n){
  lambda10n = _lambda10n;
}
void 
Run::setK0(double _k0){
  k0 = _k0;
}
void 
Run::setToleranceObtained(double _toleranceobtained){
  toleranceobtained = _toleranceobtained;
}
void 
Run::setNumberofAx1(int _numberofax1){
  numberofax1 = _numberofax1;
}
void 
Run::setNumberofAx2(int _numberofax2){
  numberofax2 = _numberofax2;
}
void 
Run::setReflectivity(double _reflectivity){
  reflectivity = _reflectivity;
}
void 
Run::setTransmittivity(double _transmittivity){
  transmittivity = _transmittivity;
}
void 
Run::setAbsorptivity(double _absorptivity){
  absorptivity = _absorptivity;
}
void 
Run::setExtinctionCrossection(double _extinctioncrosssection){
  extinctioncrosssection = _extinctioncrosssection;
}
void 
Run::setAbsorbingCrossection(double _absorbingcrosssection){
  absorbingcrosssection = _absorbingcrosssection;
}
void 
Run::setScatteringCrossection(double _scatteringcrosssection){
  scatteringcrosssection = _scatteringcrosssection;
}
void 
Run::setScatteringCrossectionWithIntegration(double _scatteringcrosssectionwithintegration){
  scatteringcrosssectionwithintegration = _scatteringcrosssectionwithintegration;
}
void 
Run::setScatteringAssymetricParam(double _scatteringassymetricparam){
  scatteringassymetricparam = _scatteringassymetricparam;
}
void 
Run::setIrra(double _irra){
  irra = _irra;
}
void 
Run::setE0(dcmplx _E0){
  E0 = _E0;
}
void 
Run::setOpticalForcex(double _opticalforcex){
  opticalforcex = _opticalforcex;
}
void 
Run::setOpticalForcey(double _opticalforcey){
  opticalforcey = _opticalforcey;
}
void 
Run::setOpticalForcez(double _opticalforcez){
  opticalforcez = _opticalforcez;
}
void 
Run::setOpticalForceModulus(double _opticalforcemodulus){
  opticalforcemodulus = _opticalforcemodulus;
}
void 
Run::setOpticalTorquex(double _opticaltorquex){
  opticaltorquex = _opticaltorquex;
}
void 
Run::setOpticalTorquey(double _opticaltorquey){
  opticaltorquey = _opticaltorquey;
}
void 
Run::setOpticalTorquez(double _opticaltorquez){
  opticaltorquez = _opticaltorquez;
}
void 
Run::setOpticalTorqueModulus(double _opticaltorquemodulus){
  opticaltorquemodulus = _opticaltorquemodulus;
}

QString
Run::getName(){
  return runname;
}
int
Run::getObjectSubunits(){
  return objectsubunits;
}
int
Run::getMeshSubunits(){
  return meshsubunits;
}
int
Run::getNmaxpp(){
  return nmaxpp;
}
double
Run::getMeshSize(){
  return meshsize;
}
double
Run::getLambda10n(){
  return lambda10n;
}
double
Run::getK0(){
  return k0;
}
double
Run::getToleranceObtained(){
  return toleranceobtained;
}
int
Run::getNumberofAx1(){
  return numberofax1;
}
int
Run::getNumberofAx2(){
  return numberofax2;
}
double
Run::getReflectivity(){
  return reflectivity;
}
double
Run::getTransmittivity(){
  return transmittivity;
}
double
Run::getAbsorptivity(){
  return absorptivity;
}
double
Run::getExtinctionCrossection(){
  return extinctioncrosssection;
}
double
Run::getAbsorbingCrossection(){
  return absorbingcrosssection;
}
double
Run::getScatteringCrossection(){
  return scatteringcrosssection;
}
double
Run::getScatteringCrossectionWithIntegration(){
  return scatteringcrosssectionwithintegration;
}
double
Run::getScatteringAssymetricParam(){
  return scatteringassymetricparam;
}
double
Run::getIrra(){
  return irra;
}
dcmplx
Run::getE0(){
  return E0;
}
double
Run::getOpticalForcex(){
  return opticalforcex;
}
double
Run::getOpticalForcey(){
  return opticalforcey;
}
double
Run::getOpticalForcez(){
  return opticalforcez;
}
double
Run::getOpticalForceModulus(){
  return opticalforcemodulus;
}
double
Run::getOpticalTorquex(){
  return opticaltorquex;
}
double
Run::getOpticalTorquey(){
  return opticaltorquey;
}
double
Run::getOpticalTorquez(){
  return opticaltorquez;
}
double
Run::getOpticalTorqueModulus(){
  return opticaltorquemodulus;
}


