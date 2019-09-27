c     programme pour faire la derivee d un faisceau gaussian
c     propagation suivant z positif
c*****************************************************************
      subroutine gaussianchampd(x,y,z,x0,y0,z0,thetat,phit,w0,k0,ss,pp
     $     ,E0,Ed,tol,nloin)
      implicit none 
      integer nloin,j
      double precision x,y,z,x0,y0,z0,a,k0,w0,sinthe,costhe,sin2the
     $     ,cos2the,xx,yy,zz,tol,rpos(3),vu(3),mat(3,3),Enorm,vn,eps,ss
     $     ,pp,s,p ,thetat,phit,theta,phi,pi
      double complex E0x,E0y,E0z,E00x,E00y,E00z,Ed(3,3),Edp(3,3),icomp
     $     ,E0s,E0p,E0
      double complex Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p,Izdx0e,Izdx2p
     $     ,Izdx2e,Ixdx,Ixdz,Izdx0,Izdx2,ctmp

c     precision machine
      eps=1.d-12
      pi=dacos(-1.d0)

      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

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

      icomp=(0.d0,1.d0)
      a=dsqrt(xx*xx+yy*yy)

      call gaussiand(a,zz,w0,k0,Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p
     $     ,Izdx0e,Izdx2p,Izdx2e,tol,nloin)

      Ixdx=Ixdxp+Ixdxe
      Ixdz=Ixdzp+Ixdze
      Izdx0=Izdx0p+Izdx0e
      Izdx2=Izdx2p+Izdx2e

      sinthe=x/a
      costhe=y/a
      if (a.eq.0.d0) then
         sinthe=0.d0
         costhe=0.d0
      endif
      sin2the=2.d0*sinthe*costhe
      cos2the=costhe*costhe-sinthe*sinthe

c     composante x
      Edp(1,1)=-E0x*sinthe*Ixdx
      Edp(1,2)=-E0x*costhe*Ixdx
      Edp(1,3)=icomp*E0x*Ixdz
      
c     composante y
      Edp(2,1)=-E0y*sinthe*Ixdx
      Edp(2,2)=-E0y*costhe*Ixdx
      Edp(2,3)=icomp*E0y*Ixdz
      
c     composante z
      Edp(3,1)=icomp*E0x*(-cos2the*Izdx2-Izdx0)+icomp*E0y*sin2the*Izdx2
      Edp(3,2)=icomp*E0y*(cos2the*Izdx2-Izdx0)+icomp*E0x*sin2the*Izdx2
      Edp(3,3)=-Ed(1,1)-Ed(2,2)
      
c     compute the final field with the inverse rotation matrix
      do j=1,3
         Ed(1,j)=mat(1,1)*Edp(1,j)+mat(2,1)*Edp(2,j)+mat(3,1)*Edp(3,j)
         Ed(2,j)=mat(1,2)*Edp(1,j)+mat(2,2)*Edp(2,j)+mat(3,2)*Edp(3,j)
         Ed(3,j)=mat(1,3)*Edp(1,j)+mat(2,3)*Edp(2,j)+mat(3,3)*Edp(3,j)
      enddo


      end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************
      subroutine gaussianchampdcirc(x,y,z,x0,y0,z0,thetat,phit,w0,k0,s
     $     ,E0,Ed,tol,nloin)
      implicit none 
      integer nloin,j
      double precision x,y,z,x0,y0,z0,a,k0,w0,sinthe,costhe,sin2the
     $     ,cos2the,xx,yy,zz,tol,rpos(3),vu(3),mat(3,3),Enorm,vn,eps,s
     $     ,thetat,phit,theta,phi,pi,us
      double complex E0x,E0y,E0z,E00x,E00y,E00z,Ed(3,3),Edp(3,3),icomp
     $     ,E0s,E0p,E0
      double complex Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p,Izdx0e,Izdx2p
     $     ,Izdx2e,Ixdx,Ixdz,Izdx0,Izdx2,ctmp

c     precision machine
      eps=1.d-12
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

      call gaussiand(a,zz,w0,k0,Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p
     $     ,Izdx0e,Izdx2p,Izdx2e,tol,nloin)

      Ixdx=Ixdxp+Ixdxe
      Ixdz=Ixdzp+Ixdze
      Izdx0=Izdx0p+Izdx0e
      Izdx2=Izdx2p+Izdx2e

      sinthe=x/a
      costhe=y/a
      if (a.eq.0.d0) then
         sinthe=0.d0
         costhe=0.d0
      endif
      sin2the=2.d0*sinthe*costhe
      cos2the=costhe*costhe-sinthe*sinthe

c     composante x
      Edp(1,1)=-E0x*sinthe*Ixdx
      Edp(1,2)=-E0x*costhe*Ixdx
      Edp(1,3)=icomp*E0x*Ixdz
      
c     composante y
      Edp(2,1)=-E0y*sinthe*Ixdx
      Edp(2,2)=-E0y*costhe*Ixdx
      Edp(2,3)=icomp*E0y*Ixdz
      
c     composante z
      Edp(3,1)=icomp*E0x*(-cos2the*Izdx2-Izdx0)+icomp*E0y*sin2the*Izdx2
      Edp(3,2)=icomp*E0y*(cos2the*Izdx2-Izdx0)+icomp*E0x*sin2the*Izdx2
      Edp(3,3)=-Ed(1,1)-Ed(2,2)
      
c     compute the final field with the inverse rotation matrix
      do j=1,3
         Ed(1,j)=mat(1,1)*Edp(1,j)+mat(2,1)*Edp(2,j)+mat(3,1)*Edp(3,j)
         Ed(2,j)=mat(1,2)*Edp(1,j)+mat(2,2)*Edp(2,j)+mat(3,2)*Edp(3,j)
         Ed(3,j)=mat(1,3)*Edp(1,j)+mat(2,3)*Edp(2,j)+mat(3,3)*Edp(3,j)
      enddo


      end

c*****************************************************************
c*****************************************************************
c*****************************************************************
      subroutine gaussiand(a,zz,w0,k0,Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p
     $     ,Izdx0e,Izdx2p,Izdx2e,tol,nloin)
      implicit none
      integer nloin
      double precision a,z,w0,k0,zz,tol
      double complex Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p,Izdx0e,Izdx2p
     $     ,Izdx2e
      double precision expe,spi,w0k0,w02,k02,z2,zp,u,v,s2
      double complex expz,icomp,uncomp,zzz,zzzf,Itmp
      logical FLAG

      if (nloin.gt.2.or.nloin.lt.0) then
         write(*,*) 'nloin has not the right value',nloin
         stop
      endif

      if (zz.lt.0.d0) then
         z=-zz
      else
         z=zz
      endif

      if (a.eq.0.d0) then
         Ixdxp=0.d0
         Ixdxe=0.d0
         Izdx2p=0.d0
         Izdx2e=0.d0

         uncomp=(1.d0,0.d0)
         icomp=(0.d0,1.d0)
         expz=cdexp(icomp*k0*z)
         expe=dexp(-k0*k0*w0*w0/2.d0)
         zp=z/w0/dsqrt(2.d0)
         spi=dsqrt(dacos(-1.d0))
         s2=dsqrt(2.d0)
         w0k0=w0*k0/s2
         w02=w0*w0
         k02=k0*k0
         z2=z*z

         zzz=icomp*zp
         call WOFZ (zzz,u,v,FLAG)
         zzzf=u*uncomp+v *icomp
         Ixdze=icomp/w02*expe*(-z+spi*(z2+w02)/s2/w0*zzzf)
         Izdx0e=(-Ixdze-icomp*w0*k02*spi/s2*expe*zzzf)/2.d0

         zzz=icomp*zp+w0k0
         call WOFZ (zzz,u,v,FLAG)
         zzzf=u*uncomp+v *icomp

         Itmp=expz*(k0-icomp*z/w02+icomp/w0*spi/s2*(1.d0+z2/w02)*zzzf)
         Ixdzp=Itmp-Ixdze
         Izdx0p=(-Itmp-icomp*spi/s2*expz*w0*k02*zzzf)/2.d0-Izdx0e

      else
         call gaussiandcalcul(a,z,k0,w0,tol,nloin,Ixdxp,Izdx2p,Ixdzp
     $        ,Izdx0p,Ixdxe,Izdx2e,Ixdze ,Izdx0e)       
      endif

      if (zz.lt.0.d0) then
         Ixdzp=dconjg(Ixdzp)
         Izdx0p=dconjg(Izdx0p)
         Ixdxp=dconjg(Ixdxp)
         Izdx2p=dconjg(Izdx2p)
         Ixdze=dconjg(Ixdze)
         Izdx0e=dconjg(Izdx0e)
         Ixdxe=dconjg(Ixdxe)
         Izdx2e=dconjg(Izdx2e)
      endif

      end
c*****************************************************************
c*****************************************************************
c*****************************************************************
      subroutine gaussiandcalcul(a,z,k0,w0,tol,nloin,Ixdxp,Izdx2p,Ixdzp
     $     ,Izdx0p,Ixdxe,Izdx2e,Ixdze,Izdx0e)
      implicit none
      integer nloin
      double precision a,z,k0,w0,aa,zz,kk0,ww0,tol     

c     declaration
      integer nlda,n,lordref,ishanks
      parameter (nlda=4)
      double precision epsabs,epsrel,borneinf,bornesup,erreur(nlda)
     $     ,intmax
      double complex icomp,Ixdxp,Ixdxe,Ixdzp,Ixdze,Izdx0p,Izdx0e,Izdx2p
     $     ,Izdx2e,result(nlda)

      external Idgaussianpropa,Idgaussianeva

      common/comgaussian/aa,zz,kk0,ww0
!$OMP THREADPRIVATE(/comgaussian/)   

      icomp=(0.d0,1.d0)

c     reinitialise a cause du common
      aa=a
      zz=z
      kk0=k0*k0
      ww0=w0*w0
      ishanks=0
      n=4
      epsrel=tol
      epsabs=0.d0

      if (nloin.eq.0.or.nloin.eq.1) then
c     PREMIERE INTEGRALE DE 0 A K0        
         borneinf=max(0.d0,kk0-100.d0/ww0)
         borneinf=dsqrt(borneinf)
         bornesup=k0

         call intgausskronrodpattersonmulti(borneinf,bornesup,ishanks
     $        ,result,erreur,epsrel,epsabs,lordref,Idgaussianpropa,nlda
     $        ,n)

         Ixdxp=result(1)
         Ixdzp=result(2)
         Izdx0p=result(3)/2.d0
         Izdx2p=result(4)/2.d0
         intmax=max(cdabs(Ixdxp),cdabs(Ixdzp),cdabs(Izdx0p)
     $        ,cdabs(Izdx2p))
      endif

      if (nloin.eq.0.or.nloin.eq.2) then
c     DEUXIEME INTEGRALE DE 0 A L INFINI
         borneinf=0.d0
c     bornesup soit avec z soit avec w0
         bornesup=min(300.d0/(z+1.d-10),10.d0/w0)
         epsabs=epsrel*intmax
         call intgausskronrodpattersonmulti(borneinf,bornesup,ishanks
     $        ,result,erreur,epsrel,epsabs,lordref,Idgaussianeva,nlda
     $        ,n)

         Ixdxe=-icomp*result(1)
         Ixdze=-icomp*result(2)
         Izdx0e=-icomp*result(3)/2.d0
         Izdx2e=-icomp*result(4)/2.d0

      endif

      end
c*******************************************************************
c*******************************************************************
c*******************************************************************
      subroutine Idgaussianpropa(kz,n,nlda,Imulti)
      implicit none
      integer n,nlda,nn,nz
      double precision a,z,k0,w0,kz,kpara
      double complex cons,Imulti(nlda)
      double precision y(3),alp,xx

      common/comgaussian/a,z,k0,w0
!$OMP THREADPRIVATE(/comgaussian/)   
      kpara=dsqrt(k0-kz*kz)
      nn=3
      alp=0.d0
      xx=a*kpara
      call DBESJ(XX,ALP,NN,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB in gaussiand with BESSEL'
         stop
      endif
    
      cons=w0*dexp(-kpara*kpara*w0/2.d0)*cdexp((0.d0,1.d0)*kz*z)
    
      Imulti(1)=cons*kpara*kz*Y(2)
      Imulti(2)=cons*kz*kz*Y(1)
      Imulti(3)=cons*kpara*kpara*Y(1)
      Imulti(4)=cons*kpara*kpara*Y(3)
      end
c*******************************************************************
c*******************************************************************
c*******************************************************************
      subroutine Idgaussianeva(kz,n,nlda,Imulti)
      implicit none
      integer n,nlda,nn,nz
      double precision a,z,k0,w0,kz,kpara,cons
      double complex Imulti(nlda)
      double precision y(3),alp,xx

      common/comgaussian/a,z,k0,w0
!$OMP THREADPRIVATE(/comgaussian/)   
      kpara=dsqrt(k0+kz*kz)
      nn=3
      alp=0.d0
      xx=a*kpara
      call DBESJ(XX,ALP,NN,Y,NZ)
      if (NZ.ne.0) then
         write(*,*) 'PB in gaussiand with BESSEL'
         stop
      endif
    
      cons=w0*dexp(-kpara*kpara*w0/2.d0)*dexp(-kz*z)
    
      Imulti(1)=cons*kpara*kz*Y(2)*(0.d0,1.d0)
      Imulti(2)=-cons*kz*kz*Y(1)
      Imulti(3)=cons*kpara*kpara*Y(1)
      Imulti(4)=cons*kpara*kpara*Y(3)
      end
