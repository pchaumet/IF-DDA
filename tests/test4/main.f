      PROGRAM test_cdmlib
#ifdef USE_HDF5
      use HDF5
#endif

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
      character(2) polarizability
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
      parameter (nfft2d=128)
      integer tabfft2(nfft2d)
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

c     constant
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)

      nxmp=1
      nymp=2
      nzmp=3
c     DATA INPUT
      lambda=632.8d0       !wavelength
      P0=1.d0              !power    
      w0=lambda*10.d0      !waist

c********************************************************
c     Define the polarizability
c********************************************************

c     polarizability='CM'  !Clausius Mossotti
c     polarizability='RR'  !Clausius Mossotti with radiative reaction
      polarizability='LA'  !Polarizability defines by Lakthakia
c     polarizability='LR'  !Polarizability defines by Draine
c     polarizability='GB'  !Polarizability with first Mie coefficient
c     polarizability='PS'  !Polarizability for a sphere with local correction
c********************************************************
c     End polarizability
c********************************************************
      
c********************************************************
c     Define the incident wave
c********************************************************
      beam='pwavelinear'        !Linear Gaussian wave
c     beam='wavelinearmulti'    !Multiple linear plane wave
c     beam='antenna'            !Antenna define by dipolar emitter
c     beam='pwavecircular'      !Circular plane wave
c     beam='gparawavecircular'  !Linear Gaussian wave with paraxial approximation
c     beam='gparawavelinear'    !Circular Gaussian wave with paraxial approximation
c     beam='gwavelinear'        !Linear Gaussian wave rigourous
c     beam='gwavecircular'      !Circular Gaussian wave rigourous
c     beam='arbitrary'          !Arbitrary wave defines by the user
c     beam='gfftwavecircular'   !Circular Gaussian wave rigourous for the first plane and propagates with FFT
c     beam='gfftwavelinear'     !Linear Gaussian wave rigourous for the first plane and propagates with FFT
      
      if (beam(1:11).eq.'pwavelinear') then
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
c*******************************************************
c     End incident wave
c*******************************************************
      
c********************************************************
c     Define the object
c********************************************************
      object='sphere'           !Sphere
c     object='inhomosphere'     !Inhomogeneous sphere
c     object='cube'             !Cube
c     object='cuboid1'          !Cuboid with side given
c     object='cuboid2'          !Cuboid with nx,ny,nz,aretecube given
c     object='inhomocuboid1'    !Inhomogeneous cuboid with side given
c     object='inhomocuboid2'    !Inhomogeneous cuboid with nx,ny,nz,aretecube given
c     object='ellipsoid'        !Ellipsoid
c     object='nspheres'         !Multiple spheres with same radius
c     object='cylinder'         !Cylinder
c     object='concentricsphere' !Concentric sphere
c     object='randomsphere1'    !Random spheres in box with side given
c     object='randomsphere2'    !Random spheres in box with nx,ny,nz,aretecube given
c     object='arbitrary'        !Arbitrary object


      
      if (object(1:6).eq.'sphere') then
         numberobjet=1
         rayonmulti(1)=500.d0         
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         nnnr=40
         epsmulti(1)=(1.5d0,0.d0)
      elseif (object(1:6).eq.'inhomosphere') then
         numberobjet=1
         rayonmulti(1)=1000.d0
         lc=lambda
         hc=0.2d0
         materiaumulti(1)='xx'
         ng=2
         nnnr=20
         epsmulti(1)=(1.01d0,0.d0)
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
         nnnr=20
         epsmulti(1)=(1.01d0,0.d0)
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
         epsmulti(1)=(1.01d0,0.d0)
      elseif (object(1:13).eq.'randomsphere1') then
         numberobjet=1
         rayonmulti(1)=100.d0         
         density=0.2d0
         sidex=2000.d0
         sidey=2000.d0
         sidez=2000.d0
         materiaumulti(1)='xx'
         ng=2
         nnnr=20
         epsmulti(1)=(1.01d0,0.d0)
      elseif (object(1:13).eq.'randomsphere2') then
         numberobjet=1
         density=0.2d0
         rayonmulti(1)=100.d0         
         aretecube=25.d0
         materiaumulti(1)='xx'
         ng=2
         epsmulti(1)=(1.01d0,0.d0)
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
         nnnr=20
      elseif (object(1:7).eq.'cuboid2') then
         numberobjet=1
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         materiaumulti(1)='xx'
         aretecube=25.d0
         epsmulti(1)=(2.000d0,0.5d0)
      elseif (object(1:8).eq.'nspheres') then
         numberobjet=2
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
         xgmulti(2)=30.d0
         ygmulti(2)=0.d0
         zgmulti(2)=0.d0
         rayonmulti(2)=15.d0
         epsmulti(2)=(1.25d0,0.d0)
         nnnr=20
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
         nnnr=20
         epsmulti(1)=(1.01d0,0.d0)
      elseif (object(1:8).eq.'cylinder') then
         rayon=10.d0
         hauteur=20.d0
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         thetaobj=45.d0
         phiobj=45.d0
         nnnr=20
         epsmulti(1)=(1.01d0,0.d0)
      elseif (object(1:16).eq.'sphereconcentric') then
         numberobjet=2      
         rayonmulti(1)=10.d0
         epsmulti(1)=(2.25d0,0.d0)
         xgmulti(1)=0.d0
         ygmulti(1)=0.d0
         zgmulti(1)=0.d0
         rayonmulti(2)=20.d0
         epsmulti(2)=(1.25d0,0.d0)        
         thetaobj=0.d0
         phiobj=0.d0
         nnnr=20
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
      endif
      if (numberobjet.gt.numberobjetmax) then
         write(*,*) 'redimensionner numberobjetmax',numberobjet
     $        ,numberobjetmax
         stop
      endif
c*******************************************************
c     End object
c*******************************************************
      
c********************************************************
c     Define if object is isotropic or not
c********************************************************
      trope='iso'
c     trope='ani'
      if (trope(1:3).eq.'ani') then
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

      
c*******************************************************
c     Define the iterative method used
c*******************************************************
      methodeit='GPBICG1'      
c     methodeit='GPBICG2'
c     methodeit='GPBICGplus'
c     methodeit='GPBICGsafe'          
c     methodeit='GPBICGAR'    
c     methodeit='GPBICGAR2'
c     methodeit='BICGSTARPLUS'
c     methodeit='GPBICOR'
c     methodeit='CORS'
c     methodeit='QMRCLA'          
c     methodeit='TFQMR'       
c     methodeit='CG'           
c     methodeit='BICGSTAB'    
c     methodeit='QMRBICGSTAB1' 
c     methodeit='QMRBICGSTAB2' 


c     Define the tolerance of the iterative method
      tolinit=1.d-4
c*******************************************************
c     End iterative method used
c*******************************************************


c********************************************************
c     define all the options
c********************************************************
      nobjet=0                  ! 1 compute only the position of the
                                ! dipole, all the other options are
                                ! disabled.
      
c     nproche adjust the size of near field domain. Near field (0)
c     inside the object, (1) inside a cuboid which contains the object,
c     (2) inside the boxnx+2*nxmp,ny+2*nymp,nz+2*nzmp
      nproche=1
      
      nxmp=0                    ! if nproche=2 used then the addsize along x : nx+2*nxmp
      nymp=0                    ! if nproche=2 used then the addsize along y : ny+2*nymp
      nzmp=0                    ! if nproche=2 used then the addsize along z : nz+2*nzmp
      nlocal=0                  ! 0 do not compute the local field, 1 compute the local field
      nmacro=0                  ! 0 do not compute the macroscopic field, 1 compute the macroscopic field

c     1 reread or create a file which contains the local field at each
c     position. Avoid to compute again the local field if the
c     configuration is the same, i.e. keep the same local field.
      nlecture=0               
      filereread='toto'         ! name fo the file if reread the local field.
      
c     nrig adjust the ways used to compute the near field. (0) compute
c     rigorously the near field in solving the near field equation, (1)
c     use renormalized Born approximation, (2) use Born approximation,
c     (3) use Born series at order 1, (4) renormalized Rytov, (5) Rytov,
c     (6) BPM, (7) renormalized BPM.
      nrig=0

      nforce=0                  ! (0) Do not compute or (1) compute the optical force.
      nforced=0                 ! (0) Do not compute or (1) compute the density of optical force.
      ntorque=0                 ! (0) Do not compute or (1) compute the optical torque.
      ntorqued=0                ! (0) Do not compute or (1) compute the density of optical torque.
      
      nsection=0                ! 0 do not compute the cross section, 1 compute the cross section. 
      nsectionsca=0             !1: calcul C_sca, Poynting and g with radiating dipole.
      nquickdiffracte=0         ! 0 compute far field classically, 1 compute far field with FFT.      
      nside=1                   ! compute microscope in reflexion (1), or transmission (0).
      nlentille=1               ! Compute microscopy.
      nquicklens=1              ! Compute microscopy with FFT (1) or wihtout FFT (0).      
      numaper= 0.9d0            ! Numerical aperture for the microscope.
      zlens=500.d0              ! Position of lens.
      ntypemic=2                ! Type of microscope: O Holographic, 1 Bright field, 2 Dark field
      gross=100.d0              ! Manyfing factor for the microscope
      numaperinc=0.8d0          ! Numerical aperture for the condenser lens.

      nenergie=0                ! 0 Do not compute energy, 1 compute energy conservation.

      nmat=0                    ! 1 Do not save the data, 0 save the data in mat file, 2 save the data in one hdf5 file.
      h5file='ifdda.h5'         ! name of the hdf5 file

      nquad=0                   !0 -> 5 define the level of integration of the Green tensor.

 
c*******************************************************
c     End options
c*******************************************************
      

c     compute size when meshsize is given
      if (object(1:13).eq.'inhomocuboid2'.or.object(1:7).eq.'cuboid2')
     $     then
         nx=nxm-2*nxmp
         ny=nym-2*nymp
         nz=nzm-2*nzmp
      else
         if (nx+2*nxmp.gt.nxm) then
            write(*,*) 'pb with size: increase nxm'
            stop
         endif
         if (ny+2*nymp.gt.nym) then
            write(*,*) 'pb with size: increase nym'
            stop
         endif
         if (nz+2*nzmp.gt.nzm) then
            write(*,*) 'pb with size: increase nzm'
            stop
         endif
      endif

      
      call cdmlib(
c     input file cdm.in
     $     lambda,beam,object,trope,
     $     materiaumulti,nnnr,tolinit,methodeit,polarizability,nquad
     $     ,nlecture,filereread,nmat,h5file,
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
     $     Tabdip,Tabmulti,tabfft2)
c     output
  
      if (nstop.eq.0) then
         write(*,*) '***********************************************'
         write(*,*) 'Computation finished without problem:'
         write(*,*) '***********************************************'
      else
         write(*,*) '***********************************************'
         write(*,*) 'Programm finished with problem:'
         write(*,*) infostr
         write(*,*) '***********************************************'
      endif
      

      end
