c$$$c     main
c$$$      implicit none
c$$$      integer nloin,nstop
c$$$      double precision x,y,z,x0,y0,z0,theta,phi,w0,pi,lambda,k0,ss,pp,P0
c$$$     $     ,irra,tol
c$$$      double complex E0,Ex,Ey,Ez
c$$$      character(64) infostr
c$$$c     donnees
c$$$      
c$$$     
c$$$      theta=45.d0
c$$$      phi=45.d0
c$$$      P0=1.d0
c$$$      pi=dacos(-1.d0)
c$$$      lambda=632.8d-9
c$$$      w0=600.d-9*100.d0
c$$$      k0=2.d0*pi/lambda
c$$$      write(*,*) 'k0w0',k0*w0
c$$$      pp=0.d0    
c$$$      P0=1.d0
c$$$      x=lambda/4.d0
c$$$      y=lambda
c$$$      z=lambda/2.d0
c$$$      x0=0.d0
c$$$      y0=0.d0
c$$$      z0=0.d0
c$$$
c$$$      call gaussianpuissance(P0,irra,w0,k0,E0)
c$$$      write(*,*) E0,irra
c$$$      call gaussianpuissancepara(w0,P0,irra,E0)
c$$$      write(*,*) E0,irra
c$$$      call irradiance(P0,w0,E0,irra)
c$$$      write(*,*) E0,irra
c$$$
c$$$      call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi,Ex,Ey,Ez,nstop
c$$$     $     ,infostr)
c$$$      write(*,*) 'E1',Ex,Ey,Ez
c$$$
c$$$      call gaussianparalinear(x,y,z,x0,y0,z0,theta,phi,w0,k0,ss,pp,E0
c$$$     $     ,Ex,Ey,Ez,nstop,infostr)
c$$$      write(*,*) 'E2',Ex,Ey,Ez
c$$$
c$$$      nloin=0
c$$$      tol=1.d-6
c$$$      call gaussianchamp(x,y,z,x0,y0,z0,theta,phi,w0,k0,ss,pp,E0,Ex,Ey
c$$$     $     ,Ez,tol,nloin,nstop,infostr)
c$$$      write(*,*) 'E3',Ex,Ey,Ez
c$$$      end
c**********************************************************
c**********************************************************
c**********************************************************
      subroutine gaussianpuissancepara(w0,P0,irra,E0)
      implicit none
      double precision w0,P0,irra,c,pi,eps0
      double complex E0

c     2*pi*w0^2
      c=299792458.d0
      pi=dacos(-1.d0)
      eps0=1.d0/(c*c*4.d0*pi*1.d-7)

      E0=dsqrt(2.d0*P0/pi/c/eps0/w0/w0)*(1.d0,0.d0)
      irra=eps0*c*cdabs(E0)**2.d0/2.d0

      end
c************************************************************
c************************************************************
c************************************************************
      subroutine gaussianparalinear(x,y,z,x0,y0,z0,thetat,phit,w0,k0,ss
     $     ,pp,E0,Ex,Ey,Ez,nstop,infostr)
      implicit none
      integer nstop
      double precision x,y,z,x0,y0,z0,theta,phi,thetat,phit,w0,k0,ss
     $     ,pp,eps,pi,s,p,vu(3),rpos(3),mat(3,3),xx,yy,zz,Enorm,vn
      double complex E0,E0s,E0p,Ex,Ey,Ez,ctmp,E00x,E00y,E00z,E0x,E0y
     $     ,E0z,Epx,Epy,Epz
      character(64) infostr


c     precision machine
      eps=1.d-12
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

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z
      if (cdabs(E00z).ge.eps*cdabs(E0)) then
         nstop=-1
         infostr='probleme in paraxial gaussian beam'
         return
      endif
c      write(*,*) 'E00',E00x,E00y,E00z
      call gaussianpara(xx,yy,zz,w0,k0,E00x,E00y,Epx,Epy,Epz)

c     compute the final field with the inverse rotation matrix
      Ex=mat(1,1)*Epx+mat(2,1)*Epy
      Ey=mat(1,2)*Epx+mat(2,2)*Epy
      Ez=mat(1,3)*Epx+mat(2,3)*Epy

      end
c************************************************************
c************************************************************
c************************************************************
      subroutine gaussianparacirc(x,y,z,x0,y0,z0,thetat,phit,w0,k0,s ,E0
     $     ,Ex,Ey,Ez,nstop,infostr)
      implicit none
      integer nstop
      double precision x,y,z,x0,y0,z0,theta,phi,thetat,phit,w0,k0,ss
     $     ,eps,pi,s,p,vu(3),rpos(3),mat(3,3),xx,yy,zz,Enorm,vn,us
      double complex E0,E0s,E0p,Ex,Ey,Ez,icomp,ctmp,E00x,E00y,E00z,E0x
     $     ,E0y,E0z,Epx,Epy,Epz
      character(64) infostr

c     precision machine
      eps=1.d-12
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0
      us=1.d0/dsqrt(2.d0)
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

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z
      if (cdabs(E00z).ge.eps*cdabs(E0)) then
         nstop=-1
         infostr='probleme in paraxial gaussian beam'
         return
      endif

      call gaussianpara(xx,yy,zz,w0,k0,E00x,E00y,Epx,Epy,Epz)
     

c     compute the final field with the inverse rotation matrix
      Ex=mat(1,1)*Epx+mat(2,1)*Epy
      Ey=mat(1,2)*Epx+mat(2,2)*Epy
      Ez=mat(1,3)*Epx+mat(2,3)*Epy

      end
c*******************************************************
c     faisceau gaussian dans l'approximation paraxiale
      subroutine gaussianpara(x,y,z,w0,k0,E0x,E0y,Ex,Ey,Ez)
      implicit none
      double precision x,y,z,w0,w02,k0
      double complex E0x,E0y,Ex,Ey,Ez,icomp,const
      double precision w,w2,z0,eta,rinv,r2

      icomp=(0.d0,1.d0)

      r2=x*x+y*y
      w02=w0*w0
      z0=k0*w02

      w2=2.d0*w02*(1.d0+z*z/z0/z0)
      w=dsqrt(w2)

      eta=datan2(z,z0)
      rinv=z/(z*z+z0*z0)
      
      const=dsqrt(2.d0)*w0/w*dexp(-r2/w2)*cdexp(icomp*k0*r2*rinv/2.d0)
     $     *cdexp(icomp*(k0*z+eta))

      Ex=E0x*const
      Ey=E0y*const
      Ez=0.d0
      
      end
c*****************************************************************
c*****************************************************************
c*****************************************************************
c************************************************************
c************************************************************
c************************************************************
      subroutine gaussianparalineard(x,y,z,x0,y0,z0,thetat,phit,w0,k0,ss
     $     ,pp,E0,Ed,nstop,infostr)
      implicit none
      integer nstop,j
      double precision x,y,z,x0,y0,z0,theta,phi,thetat,phit,w0,k0,ss
     $     ,pp,eps,pi,s,p,vu(3),rpos(3),mat(3,3),xx,yy,zz,Enorm,vn
      double complex E0,E0s,E0p,Ed(3,3),Edp(3,3),ctmp,E00x,E00y
     $     ,E00z,E0x,E0y,E0z
      character(64) infostr

c     precision machine
      eps=1.d-12
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

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z
      if (cdabs(E00z).ge.eps*cdabs(E0)) then
         nstop=-1
         infostr='probleme in paraxial gaussian beam'
         return
      endif

      call gaussianparad(xx,yy,zz,w0,k0,E00x,E00y,Edp)
     
c     compute the final field with the inverse rotation matrix
      do j=1,3
         Ed(1,j)=mat(1,1)*Edp(1,j)+mat(2,1)*Edp(2,j)+mat(3,1)*Edp(3,j)
         Ed(2,j)=mat(1,2)*Edp(1,j)+mat(2,2)*Edp(2,j)+mat(3,2)*Edp(3,j)
         Ed(3,j)=mat(1,3)*Edp(1,j)+mat(2,3)*Edp(2,j)+mat(3,3)*Edp(3,j)
      enddo


      end
c************************************************************
c************************************************************
c************************************************************
      subroutine gaussianparacircd(x,y,z,x0,y0,z0,thetat,phit,w0,k0,s
     $     ,E0,Ed,nstop,infostr)
      implicit none
      integer nstop,j
      double precision x,y,z,x0,y0,z0,theta,phi,thetat,phit,w0,k0,ss
     $     ,eps,pi,s,p,vu(3),rpos(3),mat(3,3),xx,yy,zz,Enorm,vn,us
      double complex E0,E0s,E0p,Ed(3,3),Edp(3,3),icomp,ctmp,E00x,E00y
     $     ,E00z,E0x,E0y,E0z
      character(64) infostr

c     precision machine
      eps=1.d-12
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0
      us=1.d0/dsqrt(2.d0)

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

      E00x=mat(1,1)*E0x+mat(1,2)*E0y+mat(1,3)*E0z
      E00y=mat(2,1)*E0x+mat(2,2)*E0y+mat(2,3)*E0z
      E00z=mat(3,1)*E0x+mat(3,2)*E0y+mat(3,3)*E0z
      if (cdabs(E00z).ge.eps*cdabs(E0)) then
         nstop=-1
         infostr='probleme in paraxial gaussian beam'
         return
      endif

      call gaussianparad(xx,yy,zz,w0,k0,E00x,E00y,Edp)
     
c     compute the final field with the inverse rotation matrix
      do j=1,3
         Ed(1,j)=mat(1,1)*Edp(1,j)+mat(2,1)*Edp(2,j)+mat(3,1)*Edp(3,j)
         Ed(2,j)=mat(1,2)*Edp(1,j)+mat(2,2)*Edp(2,j)+mat(3,2)*Edp(3,j)
         Ed(3,j)=mat(1,3)*Edp(1,j)+mat(2,3)*Edp(2,j)+mat(3,3)*Edp(3,j)
      enddo



      end
c     calcul de la derivee du faisceau gaussian dans l'approximation
c     paraxiale.
      subroutine gaussianparad(x,y,z,w0,k0,E0x,E0y,Ed)
      implicit none
      double precision x,y,z,w0,w02,k0
      double complex E0x,E0y,Ed(3,3),icomp,const
      double precision w,w2,z0,eta,rinv,r2

      icomp=(0.d0,1.d0)

      r2=x*x+y*y
      w02=w0*w0
      z0=k0*w02

      w2=2.d0*w02*(1.d0+z*z/z0/z0)
      w=dsqrt(w2)

      eta=datan2(z,z0)
      rinv=z/(z*z+z0*z0)
      
      const=dsqrt(2.d0)*w0/w*dexp(-r2/w2)*cdexp(icomp*k0*r2*rinv/2.d0)
     $     *cdexp(icomp*(k0*z+eta))

      Ed(1,1)=E0x*const*(-1.d0/w2+icomp*k0*rinv/2.d0)*2.d0*x
      Ed(1,2)=E0x*const*(-1.d0/w2+icomp*k0*rinv/2.d0)*2.d0*y
      Ed(1,3)=E0x*const*(icomp*k0+(icomp*z0-z)/(z*z+z0*z0)+(z*k0*k0*w02
     $     +icomp*k0*(z0*z0-z*z)/2.d0)*r2/(z*z+z0*z0)/(z*z+z0*z0))

      Ed(2,1)=E0y*const*(-1.d0/w2+icomp*k0*rinv/2.d0)*2.d0*x
      Ed(2,2)=E0y*const*(-1.d0/w2+icomp*k0*rinv/2.d0)*2.d0*y
      Ed(2,3)=E0y*const*(icomp*k0+(icomp*z0-z)/(z*z+z0*z0)+(z*k0*k0*w02
     $     +icomp*k0*(z0*z0-z*z)/2.d0)*r2/(z*z+z0*z0)/(z*z+z0*z0))

      Ed(3,1)=0.d0
      Ed(3,2)=0.d0
      Ed(3,3)=0.d0

      end
