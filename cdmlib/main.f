      use HDF5

      implicit none
c     integer
      integer ii,jj,kk,i,j,k,nstop
      integer  nlocal,nmacro,nsection,nsectionsca,nforce ,nforced
     $     ,ntorque,ntorqued,nsens,nproche,nlecture,nquickdiffracte,nrig
     $     ,nquad,nenergie,nmat

c     variables for the object
      integer nbsphere3,nbsphere,ndipole,IP(3),test,numberobjetmax
     $     ,numberobjet,ierror
      parameter (numberobjetmax=1)
      integer nx,ny,nz,nx2,ny2,nxy2,nz2,nxm,nym,nzm,nxmp,nymp,nzmp
     $     ,ntotal,nmaxpp
      integer subunit,nphi,ntheta
      parameter (nxm=50,nym=50,nzm=50,nphi=72,ntheta=36)
c     definition of the size for the code
      INTEGER nmax, ntotalm

c     variables for the positions
      integer ng
      double precision x0,y0,z0,x,y,z,xx0,yy0,zz0,rayon,density,side
     $     ,sidex,sidey,sidez,hauteur,xgmulti(numberobjetmax)
     $     ,ygmulti(numberobjetmax) ,zgmulti(numberobjetmax)
     $     ,rayonmulti(numberobjetmax),demiaxea ,demiaxeb,demiaxec
     $     ,thetaobj,phiobj,psiobj,lc,hc
      double precision aretecube
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: xs,ys,zs,xswf,yswf
     $     ,zswf
      DOUBLE PRECISION,DIMENSION((ntheta+1)*nphi)::thetafield,phifield
     $     ,poyntingfield
      double precision pi,lambda,lambda10n,k0,k03,epi,epr,c

c     variables for the material
      double precision eps0
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex epsmulti(numberobjetmax)
     $     ,epsanimulti(3,3,numberobjetmax)
      character(2) methode,cm
      character (64), DIMENSION(numberobjetmax) :: materiaumulti
      character(64) materiau,object,beam,namefileobj,namefileinc
     $     ,filereread
      character(3) trope

c     variables for the incident field and local field
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: incidentfield,
     $     localfield,macroscopicfield,forcex,forcey,forcez, torquex
     $     ,torquey,torquez
      double precision forcexmulti(numberobjetmax)
     $     ,forceymulti(numberobjetmax),forcezmulti(numberobjetmax)
     $     ,torquexmulti(numberobjetmax),torqueymulti(numberobjetmax)
     $     ,torquezmulti(numberobjetmax)
      double precision ss,pp,theta,phi,I0
      integer nbinc
      double precision thetam(10), phim(10), ppm(10), ssm(10)
      double complex E0m(10)
      double complex Eloc(3),Em(3),E0,uncomp,icomp,zzero
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldx
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldy
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldz
      double complex, dimension(nxm*nym*nzm) :: localfieldx
      double complex, dimension(nxm*nym*nzm) :: localfieldy
      double complex, dimension(nxm*nym*nzm) :: localfieldz
      double complex, dimension(nxm*nym*nzm) :: incidentfieldx
      double complex, dimension(nxm*nym*nzm) :: incidentfieldy
      double complex, dimension(nxm*nym*nzm) :: incidentfieldz
      double complex propaesplibre(3,3)
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc

c     Green function
      integer, dimension(nxm*nym*nzm) :: Tabdip,Tabmulti
      integer indice
      double complex, dimension(8*nxm*nym*nzm) :: FFTTENSORxx,
     $     FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,
     $     FFTTENSORzz,vectx,vecty,vectz

      double precision forcet(3),forcem,forcemie
      double precision couplet(3),couplem
      double complex Eder(3,3)
      
c     computation of the cross section
      integer iphi,itheta
      double precision MIECEXT,MIECABS,MIECSCA ,GSCA,Cext,normal(3)
     $     ,deltatheta,deltaphi,Csca,Cscai,Cabs,gasym,thetas,phis
     $     ,efficacite,efficaciteref,efficacitetrans
      double complex ctmp
      
c     variables for the iterative method
      INTEGER ldabi, nlar
      integer nnnr,ncompte
      integer NLIM,ndim,nou,maxit,nstat,nloop,STEPERR
      DOUBLE PRECISION  NORM,TOL,norm1,norm2,tolinit,tol1
      double complex ALPHA,BETA,GPETA,DZETA,R0RN

c     COMMON /ONTHEHEAP/ b,xr,xi,wrk
      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      
c     double complex wrk(*), xi(*), xr(*), b(*)
c     POINTER ( xr_p, xr ), ( b_p, b )
c     POINTER ( wrk_p, wrk ), ( xi_p, xi)

c     Poynting vector
      integer nr,nrmax,nw,nwmax
      double precision Poyntinginc

c     Info string
      character(64) infostr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nouvelle variable a passer en argument d'entree
c     power et diametre
      double precision P0,w0,xgaus,ygaus,zgaus,quatpieps0
      character(12) methodeit

c     nouvelle variable de sortie Irra     
      double precision irra
      

c     nouvelle variable
      integer nloin

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Creation des nouvelles variables
      integer na

c     variable pour avoir l'image a travers la lentille
      integer nlentille,nobjet,nfft2d,nfft2d2,nquicklens,ntypemic,nside
      parameter (nfft2d=512)
      double precision kx,ky,kz,deltakx,deltaky,numaper,deltax,gross
     $     ,numaperinc,zlens
      double precision kxy(nfft2d),xy(nfft2d)
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Eimageincx(nfft2d*nfft2d)
     $     ,Eimageincy(nfft2d *nfft2d) ,Eimageincz(nfft2d*nfft2d)
     $     ,Efourierx(nfft2d*nfft2d) ,Efouriery(nfft2d*nfft2d)
     $     ,Efourierz(nfft2d*nfft2d),Efourierincx(nfft2d*nfft2d)
     $     ,Efourierincy(nfft2d*nfft2d) ,Efourierincz(nfft2d*nfft2d)
     $     ,Ediffkzpos(nfft2d,nfft2d,3) ,Ediffkzneg(nfft2d ,nfft2d,3)

      character(LEN=100) :: h5file

c     data
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)

      nxmp=0
      nymp=0
      nzmp=0
c********************************************************
c     Defini les caracteristiques de l'onde incidente
c********************************************************
      lambda=632.8d0
      P0=1.d0               
      w0=lambda*10.d0       

c********************************************************
c     Defini le type d'onde incidente
c********************************************************
      beam='pwavelinear'
c     beam='wavelinearmulti'
c      beam='antenna'
c     beam='pwavecircular'
c     beam='gparawavecircular'
c      beam='gparawavelinear'
c      beam='gwavelinear'  
c     beam='gwavecircular'
c     beam='arbitrary' 
  
      
      if (beam(1:11).eq.'pwavelinear') then
c     theta=18.445351218553004d0
         theta=0.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
      elseif (beam(1:15).eq.'wavelinearmulti') then
         nbinc=1
         thetam(1)=0.d0
         phim(1)=0.d0
         ppm(1)=0.d0
         ssm(1)=1.d0
         E0m(1)=1.d0
      elseif (beam(1:7).eq.'antenna') then
         xgaus=-10000000000.d0
         ygaus=0.d0
         zgaus=0.d0
         theta=0.d0
         phi=0.d0

      elseif (beam(1:13).eq.'pwavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0
      elseif (beam(1:11).eq.'gwavelinear') then
         theta=0.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:13).eq.'gwavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:14).eq.'gfftwavelinear') then
         theta=0.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:16).eq.'gfftwavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:15).eq.'gparawavelinear') then
         theta=0.d0
         phi=0.d0
         pp=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:17).eq.'gparawavecircular') then
         theta=0.d0
         phi=0.d0
         ss=1.d0
         xgaus=0.d0
         ygaus=0.d0
         zgaus=0.d0
      elseif (beam(1:9).eq.'arbitrary') then
         namefileinc='incarbitrary.in'
      endif

c********************************************************
c     Defini le type d'objet
c********************************************************
      object='arbitrary'
c      object='inhomosphere'
c      object='randomsphere'
c     object='cube'
c     object='cuboid'
c      object='inhomocuboid'      
c      object='nspheres' 
c      object='ellipsoid'
c      object='cylinder' 
c      object='sphereconcentric'
c      object='arbitrary'

      write(*,*) 'coucou',object
      if (object(1:6).eq.'sphere') then
         numberobjet=1
         rayonmulti(1)=100.d0         
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
      elseif (object(1:6).eq.'inhomosphere') then
         numberobjet=1
         rayonmulti(1)=1000.d0
         lc=lambda
         hc=0.2d0
         materiaumulti(1)='xx'
         ng=2
      elseif (object(1:13).eq.'inhomocuboid1') then
         numberobjet=1
         sidex=2000.d0
         sidey=2000.d0
         sidez=2000.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         lc=lambda
         hc=0.2d0
         materiaumulti(1)='xx'
         ng=2
         write(*,*) 'tt',sidex,sidey,sidez,lc,hc,ng
      elseif (object(1:13).eq.'inhomocuboid2') then
         numberobjet=1
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         aretecube=25.d0
         lc=lambda
         hc=0.2d0
         materiaumulti(1)='xx'
         ng=2
         write(*,*) 'tt',sidex,sidey,sidez,lc,hc,ng

      elseif (object(1:13).eq.'randomsphere1') then
         numberobjet=1
         rayonmulti(1)=100.d0         
         sidex=2000.d0
         sidey=2000.d0
         sidez=2000.d0
         materiaumulti(1)='xx'
         ng=2
      elseif (object(1:13).eq.'randomsphere2') then
         numberobjet=1
         rayonmulti(1)=100.d0         
         aretecube=25.d0
         materiaumulti(1)='xx'
         ng=2
      elseif (object(1:4).eq.'cube') then
         numberobjet=1
         side=100.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=10
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:7).eq.'cuboid1') then
         numberobjet=1
         sidex=20.d0
         sidey=40.d0
         sidez=10.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=10
         phiobj=0.d0
         thetaobj=0.d0
         psiobj=0.d0
         epsmulti(1)=(2.000d0,0.5d0)
         write(*,*) 'cuboid'
      elseif (object(1:7).eq.'cuboid2') then
         numberobjet=1
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         aretecube=25.d0
         epsmulti(1)=(2.000d0,0.5d0)
         write(*,*) 'cuboid'
      elseif (object(1:8).eq.'nspheres') then
         numberobjet=2
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
c         xgmulti(2)=30.d0
c         ygmulti(2)=0.d0
c         zgmulti(2)=0.d0
c         rayonmulti(2)=15.d0
c         epsmulti(2)=(1.25d0,0.d0)
      elseif (object(1:9).eq.'ellipsoid') then
         demiaxea=10.d0
         demiaxeb=20.d0
         demiaxec=40.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         thetaobj=45.d0
         phiobj=45.d0
         psiobj=0.d0
      elseif (object(1:8).eq.'cylinder') then
         rayon=10.d0
         hauteur=20.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         thetaobj=45.d0
         phiobj=45.d0
      elseif (object(1:16).eq.'sphereconcentric') then
         numberobjet=2      
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
c         rayonmulti(2)=20.d0
c         epsmulti(2)=(1.25d0,0.d0)        
         thetaobj=0.d0
         phiobj=0.d0
      elseif (object(1:9).eq.'arbitrary') then
         numberobjet=1
         namefileobj='arbit.in'
         open(15,file=namefileobj,status='old',iostat=ierror)
         if (ierror.ne.0) then
            write(*,*) 'bad namefile for arbitrary'
            stop
         endif
         read(15,*) nx,ny,nz
         read(15,*) aretecube
         rewind(15)
         close(15)

         if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
            write(*,*) 'Size of the table too small'
            stop
         endif

         
c     definir ici l'objet arbitraire
      endif
      if (numberobjet.gt.numberobjetmax) then
         write(*,*) 'redimensionner numberobjetmax',numberobjet
     $        ,numberobjetmax
         stop
      endif
         
c********************************************************
c     Defini la nature de la permittivite
c********************************************************
      trope='iso'
c     trope='ani'
      if (trope(1:3).eq.'iso') then
         epsmulti(1)=(2.0d0,0.d0)
      else
         epsanimulti(1,1,1)=(2.d0,0.d0)
         epsanimulti(2,1,1)=0.d0
         epsanimulti(3,1,1)=0.d0
         epsanimulti(1,2,1)=0.d0
         epsanimulti(2,2,1)=(2.d0,0.d0)
         epsanimulti(3,2,1)=0.d0
         epsanimulti(1,3,1)=0.d0
         epsanimulti(2,3,1)=0.d0
         epsanimulti(3,3,1)=(2.d0,0.d0)
      endif
c********************************************************
c     Defini la polarisabilite choisie
c********************************************************
c      methode='CM'
      methode='RR'
c      methode='GB'
c      methode='LA'
c      methode='LR'
      
c********************************************************
c     Defini la methode iterative choisie
c********************************************************
      methodeit='GPBICG1'        
c     methodeit='GPBICG2'       
c      methodeit='GPBICGsafe'    
c      methodeit='GPBICGAR1'     
c     methodeit='GPBICGAR2'     
c      methodeit='QMRCLA
c      methodeit='TFQMR' 
c     methodeit='CG'       
c      methodeit='BICGSTAB' 
c      methodeit='QMRBICGSTAB1'
c     methodeit='QMRBICGSTAB2'
c      methodeit='GPBICOR'
c********************************************************
c     Defini la tolerance de la methode iterative
c********************************************************
      tolinit=1.d-4
c********************************************************
c     Defini la discretisation de l'obejet
c********************************************************
      nnnr=10

c********************************************************
c     Defini les calculs demandes
c******************************************************** 
      nproche=1 !0 calcul le champ dans l'objet, 1 dans le cube
                !contenant l'objet 2 dans la boite nxm,nym,nzm
      nlocal=0!  1: calcul le champ local
      nmacro=0 !  1: calcul le champ macro
      nsection=1! 1: calcul les sections efficaces
      nsectionsca=0 !1: calcul C_sca, Poynting et g par rayonnement des dipoles
      nquickdiffracte=0  !1: calcul C_sca, Poynting et g par avec la FFT
      nforce=0 ! 1: Calul la force optique
      nforced=0 ! 1: Calcul la densite de force
      ntorque=0! 1: Calul le couple optique
      ntorqued=0 ! 1: Calcul la densite de couple
      nrig=0 ! 1: calcul le champ par Born renormalise, 0 rigoureux
      nlecture=0 !1: relis les dipoles deja calcules
c      filereread='nom'
      nlentille=0 !1 Calcul l'objet vu a travers une lentille situee du
                  !cote des z positifs, foyer place a l'origine
      nquicklens=0!1: calcul avec FFT la vue a travers la lentille
      numaper=1.d0              ! ouverture numerique de la lentille
      numaperinc=0.0d0              ! ouverture numerique du condenseur
      ntypemic=0 ! type microscope 0 holo,1 brightfield, 2 darkfield phase
      nside=0 ! ! microscope in reflexion (1)  or transmission (0)
      nobjet=0                  ! 1: calcul juste la forme de l'objet
      nquad=0 !0 -> 5 defini le niveau d'integration du tenseur
      nenergie=0
      nmat=1 ! 0 écrit mat  file,  1 n ecrit pas, 2 hdf5 file
      h5file='ifdda.h5'
      gross=100.d0              ! grossissement
      zlens=0.d0 !position du foyer image de la lentille
      write(*,*) 'cdmlib',sidex,sidey,sidez

 
      call cdmlib(
c     input file cdm.in
     $     lambda,beam,object,trope,
     $     materiaumulti,nnnr,tolinit,methodeit,methode,nquad,nlecture
     $     ,filereread,nmat,h5file,
c     output file cdm.out
     $     nlocal,nmacro,nsection,nsectionsca,nquickdiffracte,nrig,
     $     nforce,nforced,ntorque,ntorqued,nproche,nlentille,nquicklens,
     $     nenergie,nobjet,
c     cube, sphere (includes multiple)
     $     density,side, sidex, sidey, sidez, hauteur,
     $     numberobjet, rayonmulti, xgmulti, ygmulti, zgmulti,
     $     epsmulti, epsanimulti,lc,hc,ng,
c     ellipsoid+arbitrary
     $     demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,
     $     namefileobj,
c     planewavecircular.in / planewavelinear.in files
     $     theta, phi, pp, ss, P0, w0, xgaus, ygaus, zgaus,namefileinc,
c     ondeplane multiple
     $     thetam, phim, ppm, ssm,E0m,nbinc,
c     return info stringf
     $     infostr, nstop,
c     return scalar results
     $     nbsphere, ndipole, aretecube,
     $     lambda10n, k0, tol1, ncompte, nloop,
     $     efficacite,efficaciteref,efficacitetrans,
     $     Cext, Cabs, Csca, Cscai, gasym, irra, E0,
     $     forcet, forcem,
     $     couplet, couplem,
     $     nxm, nym, nzm, nxmp, nymp, nzmp, nmaxpp,
     $     incidentfield, localfield, macroscopicfield,
     $     xs, ys, zs, xswf, yswf, zswf,
     $     ntheta, nphi, thetafield, phifield, poyntingfield,
     $     forcex,forcey,forcez,forcexmulti,forceymulti,forcezmulti,
     $     torquex,torquey,torquez,torquexmulti,torqueymulti,
     $     torquezmulti,
     $     incidentfieldx, incidentfieldy, incidentfieldz,
     $     localfieldx, localfieldy, localfieldz,
     $     macroscopicfieldx, macroscopicfieldy, macroscopicfieldz,
     $     polarisa,epsilon,
     $     nfft2d,Eimagex,Eimagey,Eimagez,Eimageincx,Eimageincy,
     $     Eimageincz,Efourierx,Efouriery,Efourierz,Efourierincx,
     $     Efourierincy,Efourierincz,kxy,xy,numaper,numaperinc,gross,
     $     zlens,ntypemic,nside,
c****************************************************
c     tableaux utilises que dans cdmlib
c****************************************************
c     taille double complex (3*nxm*nym*nzm)
     $     FF,FF0,FFloc,xr,xi,
c     taille double complex (3*nxm*nym*nzm,12)
     $     wrk,
c     taille double complex (8*nxm*nym*nzm)
     $     FFTTENSORxx, FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,
     $     FFTTENSORzz,vectx,vecty,vectz,
c     taille double complex (nfft2d,nfft2d,3)
     $     Ediffkzpos,Ediffkzneg,
c     taille entier (nxm*nym*nzm)
     $     Tabdip,Tabmulti)
c     output
      if (nstop.eq.1) then
         write(*,*) infostr
         stop
      endif
      write(*,*) 'infostr',  infostr
      pi=dacos(-1.d0)
      if (materiau.ne.'xx') then     
         write(*,*) 'Relative permittivity',epsmulti(1)
      else
         if (trope.eq.'iso') then          
            write(*,*) 'Relative permittivity',epsmulti(1)
         else 
            do i=1,3
               do j=1,3
                  write(*,*) 'Relative permittivity',epsanimulti(i,j,1)
     $                 ,i,j
               enddo
            enddo
         endif
      endif

      write(*,*) 'Object under study ',object
      write(*,*) 'Nombre Object',numberobjet
      write(*,*) 'number of subunit for the object',nbsphere
      write(*,*) 'number of subunit for the mesh ',ndipole
      write(*,*) 'mesh size',aretecube
      write(*,*) 'lambda/(10n)',lambda/10.d0/cdabs(cdsqrt(epsmulti(1)))

      write(*,*) '******* Compute the incident field *******'
      write(*,*) 'Beam used',beam
      write(*,*) 'k0=',k0     
      write(*,*) 'theta=',theta
      write(*,*) 'phi=',phi
      write(*,*) 'Irradiance',Irra
      write(*,*) 'Field',E0
      I0=cdabs(E0)**2

      write(*,*) '***** Solve the linear system *****'
      write(*,*) 'Tolerance asked for the iterative method   ',tolinit
      write(*,*) 'Tolerance obtained for the iterative method',tol1
      write(*,*) 'Number of product Ax for the iterative method'
     $     ,ncompte,nloop

      if (nsection.eq.1) then      
         rayon=rayonmulti(1)*1.d-9
         write(*,*) 'mie',epsmulti(1),rayon,lambda
         CALL CALLBHMIE(1.d0,epsmulti(1),rayon,lambda,MIECEXT,MIECABS
     $        ,MIECSCA,GSCA)
         write(*,*) 'MIECEXT',MIECEXT,MIECABS,MIECSCA,GSCA
         write(*,*) 'force',(MIECEXT-GSCA*MIECSCA)/8.d0/pi*quatpieps0

         write(*,*) 'extinction cross section',Cext
         write(*,*) 'absorbing cross section ',Cabs
         write(*,*) 'scattering cross section',Csca
         write(*,*) 'cos',gasym
      endif
      if (nsectionsca.eq.1) then
         write(*,*) 'scattering cross section with integration',Cscai
         write(*,*) 'scattering asymetric parameter',gasym
      endif

      if (nforce.eq.1) then
         write(*,*) '****** Compute the optical force *********'
         write(*,*) 'optical force x',forcet(1)
         write(*,*) 'optical force y',forcet(2)
         write(*,*) 'optical force z',forcet(3)
         
         forcemie=(MIECext-GSCA*MIECsca)/8.d0/pi*I0*quatpieps0
         write(*,*) 'modulus of the force',forcem,'Mie',forcemie
      endif
      if (ntorque*nforce.eq.1) then
         write(*,*) '********* Compute the optical torque *********'
         write(*,*) 'optical torque x',couplet(1)
         write(*,*) 'optical torque y',couplet(2)
         write(*,*) 'optical torque z',couplet(3)
         write(*,*) 'modulus of the optical torque',couplem
         write(*,*) 'couple Mie',MIECABS/8.d0/k0/pi*I0*quatpieps0
      endif
      if (numberobjet.ne.1) then
         do i=1,numberobjet
            write(*,*) 'forcex',forcexmulti(i)
            write(*,*) 'forcey',forceymulti(i)
            write(*,*) 'forcez',forcezmulti(i)
         enddo
      endif

      end
