#ifndef CDMLIB_H
#define CDMLIB_H

extern"C" {
int cdmlib_( // input file cdm.in
	double *lambda,	char *beam, char *object, char *trope, 
	char materiau[][64], int  *nnnr, double *tolinit, char* methodeit,
        char* methode, int *integ,
        int *nlecture, char *filereread, int* nmat, char* fileh5,
	// ouput file cdm.out
	int *nlocal,int *nmacro, int *nsection,	int *nsectionsca, 
        int *nquickdiffracte, int *nrig, 
	int *nforce, int *nforced, int *ntorque, int *ntorqued,
        int *nproche, int *nlentille, int *nquicklens, int *nenergie, int *nobjet,
        // sphere.in / cube.in files (includes multiple spheres)
        double *density, double *side, double *sidex, double *sidey,  double *sidez, double *hauteur,
        int *numberobjet, double *rayon, double *xg, double *yg, double *zg,
	dcmplx *eps, dcmplx *epsani, double *lc, double *hc, int *ng,
	// ellipsoid
        double *demiaxea, double *demiaxeb, double *demiaxec, 
	double *thetaobj, double *phiobj, double *psiobj, char* namefileobj,
	// planewavecircular.in / planewavecircular.in files
	double *theta, double *phi, double *pp, double *ss, double *P0, double *W0,
        double *xgaus, double *ygaus, double *zgaus, char* namefileinc,
        double *thetam, double *phim, double *ppm, double *ssm, dcmplx *E0m, int *nbinc,
	// return info string
	char *infostr, int *nstop,
	// return scalar results
	int *subunitsobject, int *subunitsmesh, double *meshsize,
        double *lambda10n, double *k0, double *toleranceobtained,
        int *numberofax1, int *numberofax2,
        double *absorptivity, double *reflectivity, double *transmittivity,
        double *extinctioncrosssection, double *absorbingcrosssection, 
        double *crosssection, double *crosssectionpoynting, double *assymetricparam, 
	double *irra, dcmplx *E0,
        double *opticalforce, double *opticalforcedensity,
        double *opticaltorque, double *opticaltorquedensity,
        int *nxm, int *nym, int *nzm, int *nxmp, int *nymp, int *nzmp, int *nmaxpp,
        double *incidentfield, double *localfield, double *macroscopicfield,
	double *xc, double *yc, double *zc, double *xcwf, double *ycwf, double *zcwf,
	int *ntheta, int *nphi, double *thetafield, double *phifield, double *poyntingfield,
	double *forcex, double *forcey, double *forcez,
        double *forcexmulti, double *forceymulti, double *forcezmulti,
        double *torquex, double *torquey, double *torquez,
        double *torquexmulti, double *torqueymulti, double *torquezmulti,
        dcmplx *incidentfieldx, dcmplx *incidentfieldy, dcmplx *incidentfieldz,
	dcmplx *localfieldx, dcmplx *localfieldy, dcmplx *localfieldz,
        dcmplx *macroscopicfieldx, dcmplx *macroscopicfieldy, dcmplx *macroscopicfieldz,
	dcmplx *polarisa, dcmplx *epsilon, int *nfft2d, dcmplx *eimagex, dcmplx *eimagey,
        dcmplx *eimagez, dcmplx *eimageincx, dcmplx *eimageincy, dcmplx *eimageincz,
        dcmplx *efourierx, dcmplx *efouriery, dcmplx *efourierz,
	dcmplx *efourierincx, dcmplx *efourierincy, dcmplx *efourierincz,
        double *kxy, double *xy, double *numaper, double *numaperinc, double *gross, double *zlens, int *ntypemic, int *nside,
//****************************************************
//     tableaux utilises que dans cdmlib
//****************************************************
//     taille double complex (3*nxm*nym*nzm)
        dcmplx *FF, dcmplx *FF0, dcmplx *FFloc, dcmplx *xr, dcmplx *xi,
//     taille double complex (3*nxm*nym*nzm,12)
        dcmplx *wrk,
//     taille double complex (8*nxm*nym*nzm)
        dcmplx *FFTTENSORxx, dcmplx *FFTTENSORxy, dcmplx *FFTTENSORxz, 
        dcmplx *FFTTENSORyy, dcmplx *FFTTENSORyz,
        dcmplx *FFTTENSORzz, dcmplx *vectx, dcmplx *vecty, dcmplx *vectz,
//     taille double complex (nfft2d,nfft2d,3)
        dcmplx *Ediffkzpos, dcmplx *Ediffkzneg,
//     taille entier (nxm*nym*nzm)
        int *Tabdip, int *Tabmulti,
//     taille entier (nfft2d)
        int *Tabfft2	
	);
}

#endif
