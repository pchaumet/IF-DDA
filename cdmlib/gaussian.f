c     Authors: P. C. Chaumet 
c     Date: 03/10/2009

c     Purpose: These routines compute a Gaussian beam in the case of
c     highly nonparaxial field ($w_0<\lambda$). The computation is done
c     in the vectorial case for a polarized Gaussian beam and is valid
c     whatever the position: far away or close to/in the waist. It is
c     also computed from the the power of the Gaussian beam the
c     intensity and the irradiance of the field at the origin. Notice
c     that the origin of the beam is taken at the center of the
c     waist. The profile of the electric field at the origin is
c     E=E0e^{-rho^2/2/w_0^2}
c     
c     Purpose (french): Cet ensemble de routines calcul le champ
c     electric d'un faisceau Gaussian dans le cas general d'un faisceau
c     pouvant etre tres fortement non paraxial ($w_0<\lambda$). Le
c     calcul est fait dans le cas vectoriel a trois dimensions pour un
c     faisceau polarise.  Le calcul du faisceau gaussian est base sur un
c     developpement du champ en spectre d'onde plane. La position
c     d'observation peut etre a priori n'importe ou loin comme
c     proche/dans le wasit; teste jusqu'a x longeur d'onde: En champ
c     lontain des approximations analytiques se trouvent dans la
c     litterature. Il existe aussi une routine qui sachant la puissance
c     en watt du faisceau Gaussian donne l'intensite du champ a
c     l'origine ainsi que l'irradiance a cette meme position. A noter
c     que quand je parle d'origine il s'agit toujours du centre du waist
c     du faisceau Gaussian. Le profil du champ electrique dans le waist
c     s'ecrit: E=E0e^{-rho^2/2/w_0^2}


c     Reference: if you use this routine in your research, please
c     reference, as appropriate: P. C. Chaumet, Fully vectorial highly
c     non paraxial beam close to the waist. J. Opt. Soc. Am. A 23, 3197
c     (2006) and P. C. Chaumet, B. Pouligny, R. Dimova and N. Sojic,
c     Optical tweezers in interaction with an apertureless probe
c     J. Appl. Phys. 102, 024915 (2007).

c     license: GNU GPL

c*******************************************************************
c*******************************************************************

c     Computation of the intensity and irradiance of the beam at the
c     origin from the power of the Gaussian beam.

c     P: power (Watt)

c     irra: Irradiance (Watt/m^2)

c     w0: waist (m)

c     k0: wavenumber (m^{-1})

c     Igauss: Intensity of the field. Then the field is E0=dsqrt(Igauss)
c     (V/m)

c*******************************************************************
c*******************************************************************

      subroutine gaussianpuissance(P,irra,w0,k0,E0)
      implicit none
      double precision P,w0,k0,Igauss,c,pi,u,v,const,const0,const1,eps0
      double complex icomp,zarg,uncomp,zf,E0
      logical flag
      double precision Ipower,ww0,irra

      external Ipower
      common/comgaussian0/ww0
      ww0=w0

      c=299792458.d0
      pi=dacos(-1.d0)
      eps0=1.d0/(c*c*4.d0*pi*1.d-7)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      
      const0=pi/4.d0*w0*w0*c*eps0
      const1=k0*w0

      zarg=const1*uncomp
      call WOFZ (zarg, U, V, FLAG)
      zf=(u*uncomp+v*icomp)
c      if (dimag(zf).ne.0.d0) write(*,*) 'puissance imag',zarg,zf
      
      const=const0*(1.d0+(const1*const1-1.d0)/const1*dsqrt(pi)/2.d0
     $     *dimag(zf))

      Igauss=P/const

c     Compute irradiance
      
      zarg=const1*uncomp/2.d0
      call WOFZ (zarg, U, V, FLAG)
      zf=(u*uncomp+v*icomp)
c      if (dimag(zf).ne.0.d0) write(*,*) 'irradiance imag',zarg,zf     
      E0=uncomp*dsqrt(Igauss)


      zarg=const1*uncomp/dsqrt(2.d0)
      call WOFZ (zarg, U, V, FLAG)
      zf=(u*uncomp+v*icomp)
      irra=c*eps0/2.d0*(1.d0/2.d0-dsqrt(2.d0*pi)*(const1*const1
     $     -1.d0)/const1/4.d0*dreal(icomp*(cdexp(-zarg*zarg)
     $     -dconjg(zf))))*Igauss

      end

c*****************************************************************
c*****************************************************************

c     subroutine to call to compute the Gaussian beam.
c     The origin is defined at the center of the waist

c     w0: waist (m)

c     k0: wavenumber (m^{-1})

c     E0, magnitude of the field at the origin, one can uses the routine
c     gaussianpuissance to know this field if only the power of the
c     Gaussian beam is known. (V/m)

c     x,y,z location where the electric field is wanted (m).

c     x0,y0,z0 location of the center of the waist (origin) (m).

c     Ex,Ey,Ez electric field at the position x,y,z (V/m)

c     vu(3) unit vector which gives the direction of the Gaussian
c     beam. It should be perpendicular to the field at the origin.

c     tol: tolerance used when the integration is performed

c     nloin: in the angular spectrum representation one can uses only
c     the propagating waves (0), the evanescent waves (1) or both (2).

c*****************************************************************
c*****************************************************************

      subroutine gaussianchamp(x,y,z,x0,y0,z0,thetat,phit,w0,k0,ss,pp,E0
     $     ,Ex,Ey,Ez,tol,nloin,nstop,infostr)
      implicit none
      integer nloin,nstop
      double precision x,y,z,xx,yy,zz,x0,y0,z0,a,k0,w0,sinthe,costhe
     $     ,sin2the,cos2the,tol,rpos(3),vu(3),Enorm,theta,phi,vn,eps
     $     ,mat(3,3),thetat,phit,ss,pp,s,p,pi,Em
      double complex E0x,E0y,E0z,E00x,E00y,E00z,Ex,Ey,Ez,Epx,Epy,Epz
     $     ,Ixxp,Ixxe,Izzp,Izze,Ixx,Izz,icomp,ctmp,E0,E0s,E0p
      character(64) infostr

      call gaussianparalinear(x,y,z,x0,y0,z0,thetat,phit,w0,k0,ss,pp,E0
     $     ,Ex,Ey,Ez,nstop,infostr)
      Em=dsqrt(cdabs(Ex)+cdabs(Ey)+cdabs(Ez))
c      write(*,*) 'gauss',x,y,z,x0,y0,z0,Em,cdabs(E0)
      if (Em.le.cdabs(E0)*1.d-7) then
         Ex=0.d0
         Ey=0.d0
         Ez=0.d0
         return
      endif

      
c     precision machine
      eps=1.d-8
      pi=dacos(-1.d0)

      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

      if (pp.lt.0.d0.or.pp.gt.1.d0) then
         nstop=-1
         infostr='problem in plance incident wave'
         return
      endif

      p=dsqrt(pp)
      s=dsqrt(1.d0-pp)

      E0s=E0*s
      E0p=E0*p
      vu(1)=dsin(theta)*dcos(phi)
      vu(2)=dsin(theta)*dsin(phi)
      vu(3)=dcos(theta)

      E0y=E0s*dcos(phi)+E0p*dcos(theta)*dsin(phi)
      E0x=-E0s*dsin(phi)+E0p*dcos(theta)*dcos(phi)      
      E0z=-E0p*dsin(theta)

c     check if the unit vector is perpendicular to the incident field.
      Enorm=dsqrt(cdabs(E0x)*cdabs(E0x)+cdabs(E0y)*cdabs(E0y)+cdabs(E0z)
     $     *cdabs(E0z))
      ctmp=E0x*vu(1)+E0y*vu(2)+E0z*vu(3)
      if (cdabs(ctmp)/Enorm.ge.eps) then
         nstop=-1
         infostr='unit vector not perpendicular to the incident field'
         return
      endif
c     check if vv is a unit vector
      vn=vu(1)*vu(1)+vu(2)*vu(2)+vu(3)*vu(3)
      if (vn.le.1.d0-eps.or.vn.ge.1.d0+eps) then
         nstop=-1
         infostr='vu is not a unit vector'
         return
      endif

c     translation
      rpos(1)=x-x0
      rpos(2)=y-y0
      rpos(3)=z-z0
c     compute the matrix of rotation from the unit vector
      vn=vu(1)*vu(1)+vu(2)*vu(2)
      if (vn.le.eps) then
         mat(1,1)=1.d0
         mat(1,2)=0.d0
         mat(1,3)=0.d0
         mat(2,1)=0.d0
         mat(2,2)=1.d0
         mat(2,3)=0.d0
         mat(3,1)=0.d0
         mat(3,2)=0.d0
         mat(3,3)=1.d0
      else
         mat(1,1)=vu(3)+(1.d0-vu(3))*vu(2)*vu(2)/vn
         mat(1,2)=-(1.d0-vu(3))*vu(1)*vu(2)/vn
         mat(1,3)=-vu(1)
         mat(2,1)=mat(1,2)
         mat(2,2)=vu(3)+(1.d0-vu(3))*vu(1)*vu(1)/vn
         mat(2,3)=-vu(2)
         mat(3,1)=vu(1)
         mat(3,2)=vu(2)
         mat(3,3)=vu(3)
      endif

c     compute the new incident field and position
      xx=mat(1,1)*rpos(1)+mat(1,2)*rpos(2)+mat(1,3)*rpos(3)
      yy=mat(2,1)*rpos(1)+mat(2,2)*rpos(2)+mat(2,3)*rpos(3)
      zz=mat(3,1)*rpos(1)+mat(3,2)*rpos(2)+mat(3,3)*rpos(3)
      
c      write(*,*) 'E',E0x,E0y,E0z

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z

c      write(*,*) 'E',E00x,E00y,E00z

      icomp=(0.d0,1.d0)
      a=dsqrt(xx*xx+yy*yy)
c      write(*,*) a,zz,w0,k0,tol,nloin
      call gaussian(a,zz,w0,k0,Ixxp,Ixxe,Izzp,Izze,tol,nloin)
      
      Ixx=Ixxp+Ixxe
      Izz=Izzp+Izze

      sinthe=xx/a
      costhe=yy/a
      if (a.eq.0.d0) then
         sinthe=0.d0
         costhe=0.d0
      endif
      sin2the=2.d0*sinthe*costhe
      cos2the=costhe*costhe-sinthe*sinthe
      
      Epx=Ixx*E00x
      Epy=Ixx*E00y
      Epz=-icomp*Izz*(E00x*sinthe+E00y*costhe)

c     compute the final field with the inverse rotation matrix
      Ex=mat(1,1)*Epx+mat(2,1)*Epy+mat(3,1)*Epz
      Ey=mat(1,2)*Epx+mat(2,2)*Epy+mat(3,2)*Epz
      Ez=mat(1,3)*Epx+mat(2,3)*Epy+mat(3,3)*Epz

      end
c*****************************************************************
c*****************************************************************
c*****************************************************************
      subroutine gaussianchampcirc(x,y,z,x0,y0,z0,thetat,phit,w0,k0,s,E0
     $     ,Ex,Ey,Ez,tol,nloin,nstop,infostr)
      implicit none
      integer nloin,nstop
      double precision x,y,z,xx,yy,zz,x0,y0,z0,a,k0,w0,sinthe,costhe
     $     ,sin2the,cos2the,tol,rpos(3),vu(3),Enorm,theta,phi,vn,eps
     $     ,mat(3,3),thetat,phit,s,pi,us,Em
      double complex E0x,E0y,E0z,E00x,E00y,E00z,Ex,Ey,Ez,Epx,Epy,Epz
     $     ,Ixxp,Ixxe,Izzp,Izze,Ixx,Izz,icomp,ctmp,E0,E0s,E0p
      character(64) infostr
      
      call gaussianparacirc(x,y,z,x0,y0,z0,thetat,phit,w0,k0,s,E0 ,Ex,Ey
     $     ,Ez,nstop,infostr)
      Em=dsqrt(cdabs(Ex)+cdabs(Ey)+cdabs(Ez))
      if (Em.le.cdabs(E0)*1.d-7) then
         Ex=0.d0
         Ey=0.d0
         Ez=0.d0
         return
      endif
      
c     precision machine
      eps=1.d-8
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      us=1.d0/dsqrt(2.d0)
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

      E0s=E0/dsqrt(2.d0)
      E0p=E0*icomp*dsign(us,s)

      vu(1)=dsin(theta)*dcos(phi)
      vu(2)=dsin(theta)*dsin(phi)
      vu(3)=dcos(theta)

      E0y=E0s*dcos(phi)+E0p*dcos(theta)*dsin(phi)
      E0x=-E0s*dsin(phi)+E0p*dcos(theta)*dcos(phi)      
      E0z=-E0p*dsin(theta)

c     check if the unit vector is perpendicular to the incident field.
      Enorm=dsqrt(cdabs(E0x)*cdabs(E0x)+cdabs(E0y)*cdabs(E0y)+cdabs(E0z)
     $     *cdabs(E0z))
      ctmp=E0x*vu(1)+E0y*vu(2)+E0z*vu(3)
      if (cdabs(ctmp)/Enorm.ge.eps) then
         write(*,*)
     $        'unit vector not perpendicular to the incident field'
         stop
      endif
c     check if vv is a unit vector
      vn=vu(1)*vu(1)+vu(2)*vu(2)+vu(3)*vu(3)
      if (vn.le.1.d0-eps.or.vn.ge.1.d0+eps) then
         write(*,*) 'vu is not a unit vector'
         stop
      endif

c     translation
      rpos(1)=x-x0
      rpos(2)=y-y0
      rpos(3)=z-z0
c     compute the matrix of rotation from the unit vector
      vn=vu(1)*vu(1)+vu(2)*vu(2)
      if (vn.le.eps) then
         mat(1,1)=1.d0
         mat(1,2)=0.d0
         mat(1,3)=0.d0
         mat(2,1)=0.d0
         mat(2,2)=1.d0
         mat(2,3)=0.d0
         mat(3,1)=0.d0
         mat(3,2)=0.d0
         mat(3,3)=1.d0
      else
         mat(1,1)=vu(3)+(1.d0-vu(3))*vu(2)*vu(2)/vn
         mat(1,2)=-(1.d0-vu(3))*vu(1)*vu(2)/vn
         mat(1,3)=-vu(1)
         mat(2,1)=mat(1,2)
         mat(2,2)=vu(3)+(1.d0-vu(3))*vu(1)*vu(1)/vn
         mat(2,3)=-vu(2)
         mat(3,1)=vu(1)
         mat(3,2)=vu(2)
         mat(3,3)=vu(3)
      endif

c     compute the new incident field and position
      xx=mat(1,1)*rpos(1)+mat(1,2)*rpos(2)+mat(1,3)*rpos(3)
      yy=mat(2,1)*rpos(1)+mat(2,2)*rpos(2)+mat(2,3)*rpos(3)
      zz=mat(3,1)*rpos(1)+mat(3,2)*rpos(2)+mat(3,3)*rpos(3)
      
c      write(*,*) 'E',E0x,E0y,E0z

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z

c      write(*,*) 'E',E00x,E00y,E00z

      a=dsqrt(xx*xx+yy*yy)
c      write(*,*) a,zz,w0,k0,tol,nloin
      call gaussian(a,zz,w0,k0,Ixxp,Ixxe,Izzp,Izze,tol,nloin)
      
      Ixx=Ixxp+Ixxe
      Izz=Izzp+Izze

      sinthe=xx/a
      costhe=yy/a
      if (a.eq.0.d0) then
         sinthe=0.d0
         costhe=0.d0
      endif
      sin2the=2.d0*sinthe*costhe
      cos2the=costhe*costhe-sinthe*sinthe
      
      Epx=Ixx*E00x
      Epy=Ixx*E00y
      Epz=-icomp*Izz*(E00x*sinthe+E00y*costhe)

c     compute the final field with the inverse rotation matrix
      Ex=mat(1,1)*Epx+mat(2,1)*Epy+mat(3,1)*Epz
      Ey=mat(1,2)*Epx+mat(2,2)*Epy+mat(3,2)*Epz
      Ez=mat(1,3)*Epx+mat(2,3)*Epy+mat(3,3)*Epz

      end
c*****************************************************************
c*****************************************************************
c*****************************************************************
      subroutine gaussian(a,zzz,w0,k0,Ixxp,Ixxe,Izzp,Izze,tol,nloin)
      implicit none
      integer nloin,ishanks,lordref,nlda,n
      parameter (nlda=2)
      double precision a,z,k0,w0,aa,zz,kk0,ww0,zzz
      double precision erreur(nlda)
      double complex Ixxp,Ixxe,Izzp,Izze,result(nlda)

      double precision epsabs,epsrel,borneinf,bornesup,tol
      external Ipropagaussian,Ievagaussian

      common/comgaussian/aa,zz,kk0,ww0
!$OMP THREADPRIVATE(/comgaussian/)
      if (nloin.gt.2.or.nloin.lt.0) then
         write(*,*) 'not right value of nloin',nloin
         stop
      endif

      if (zzz.lt.0.d0) then
         z=-zzz
      else
         z=zzz
      endif   

c     reinitialise a cause du common
      aa=a
      zz=z
      kk0=k0*k0
      ww0=w0*w0
      epsrel=tol
      epsabs=0.d0
      ishanks=0
      n=2

      if (nloin.eq.0.or.nloin.eq.1) then
c     PREMIERE INTEGRALE DE 0 A K0        
         borneinf=max(0.d0,kk0-100.d0/ww0)
         borneinf=dsqrt(borneinf)          
c         write(*,*) 'borneinf',borneinf,k0
         bornesup=k0
        
         call intgausskronrodpattersonmulti(borneinf,bornesup,ishanks
     $        ,result,erreur,epsrel,epsabs,lordref,Ipropagaussian,nlda
     $        ,n)
      endif           
      Ixxp=result(1)
      Izzp=result(2)

      if (nloin.eq.0.or.nloin.eq.2) then
c     DEUXIEME INTEGRALE DE 0 A L INFINI        
         borneinf=0.d0
         bornesup=min(15.d0/(z+1.d-10),10.d0/w0)
        
         call intgausskronrodpattersonmulti(borneinf,bornesup,ishanks
     $        ,result,erreur,epsrel,epsabs,lordref,Ievagaussian
     $        ,nlda,n) 
         Ixxe=result(1)
         Izze=-(0.d0,1.d0)*result(2)
      endif

      if (zzz.lt.0.d0) then
         Ixxp=dconjg(Ixxp)
         Izzp=dconjg(Izzp)
         Ixxe=dconjg(Ixxe)
         Izze=dconjg(Izze)
      endif

      end
c*******************************************************************
c*******************************************************************
c*******************************************************************
      subroutine Ipropagaussian(kz,n,nlda,Imulti)
      implicit none
      integer n,nlda
      double precision a,z,k0,w0,kz,kpara,dbesj0,dbesj1
      double complex cons,Imulti(nlda)

      common/comgaussian/a,z,k0,w0
!$OMP THREADPRIVATE(/comgaussian/)     
      kpara=dsqrt(k0-kz*kz)
      cons=w0*dexp(-kpara*kpara*w0/2.d0)*cdexp((0.d0,1.d0)*kz*z)    
      Imulti(1)=cons*kz*dbesj0(a*kpara)
      Imulti(2)=cons*kpara*dbesj1(a*kpara)

      end
c*******************************************************************
c*******************************************************************
c*******************************************************************
      subroutine Ievagaussian(kz,n,nlda,Imulti)
      implicit none
      integer n,nlda
      double precision a,z,k0,w0,kz,kpara,dbesj0,dbesj1,cons
      double complex Imulti(nlda)

      common/comgaussian/a,z,k0,w0
!$OMP THREADPRIVATE(/comgaussian/)           
      kpara=dsqrt(k0+kz*kz)
      cons=w0*dexp(-kpara*kpara*w0/2.d0)*dexp(-kz*z)    
      Imulti(1)=cons*kz*dbesj0(a*kpara)
      Imulti(2)=cons*kpara*dbesj1(a*kpara)

      end
