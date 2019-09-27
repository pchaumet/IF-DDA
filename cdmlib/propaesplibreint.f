c     Integration of the Green Tensor
      subroutine propaespacelibreint(xa,ya,za,x0a,y0a,z0a,k0a
     $     ,arretecube,lim,RELREQ,propaint)
      implicit none
      integer i,j
c     definition of the position of the dipole, observation, wavenumber
c     ,wavelength, spacing lattice
      double precision xa,ya,za,x0a,y0a,z0a,k0a,pi,arretecubem,lim
      double precision x,y,z,x0,y0,z0,arretecube,k0,xx0,yy0,zz0,lambda

c     Green tensor and complex number i,1 and constant
      double complex propaint(3,3),icomp,uncomp,const1,const2

c     Variables needs for the integration
      integer  KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
      parameter (nw=4000000,ndim=3,nf=12)
      double precision A(NDIM), B(NDIM), WRKSTR(NW) 
      double precision  ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
      
      double precision Id(3,3),Rab,Rtenseur(3,3),Rvect(3)

      external fonctionigtadda

      common/k0xyz/k0,x,y,z,xx0,yy0,zz0

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)

      x=xa
      y=ya
      z=za
      x0=x0a
      y0=y0a
      z0=z0a
      k0=k0a
      lambda=2.d0*pi/k0
      arretecubem=arretecube*0.1d0
      Rab=dsqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))

c     test if the dipole and the observation are confounded
      if (Rab.le.arretecubem) then
c     We use here the integration given by Lakthakia
         do i=1,3
            do j=1,3
               propaint(i,j)=0.d0
            enddo
         enddo
      elseif (Rab.ge.lim*lambda) then
c     As the observation is far away the dipole, we do not perform the
c     integration
         Rvect(1)=(x-x0)
         Rvect(2)=(y-y0)
         Rvect(3)=(z-z0)
         do i=1,3
            do j=1,3
               Id(i,j)=0.d0
               if (i.eq.j) Id(i,i)=1.d0
               Rtenseur(i,j)=Rvect(i)*Rvect(j)
            enddo
         enddo
         do i=1,3
            do j=1,3
               Rtenseur(i,j)=Rtenseur(i,j)/(Rab*Rab)
            enddo
         enddo         
         const1=(Rab*uncomp)**(-3.d0)-icomp*k0*(Rab**(-2.d0))
         const2=k0*k0/Rab*uncomp
         do i=1,3
            do j=1,3
               propaint(i,j)=((3.d0*Rtenseur(i,j)-Id(i,j))*const1+(Id(i
     $              ,j)-Rtenseur(i,j))*const2)* cdexp(icomp*k0*Rab)
            enddo
         enddo
         
      else
c     We perform the integration of the tensor
c     defintion for the integration
         MINCLS = 1000
         MAXCLS = 1000000
         KEY = 0
         ABSREQ = 0.d0
         
         A(1)=x0-arretecube/2.d0
         A(2)=y0-arretecube/2.d0
         A(3)=z0-arretecube/2.d0
         B(1)=x0+arretecube/2.d0
         B(2)=y0+arretecube/2.d0
         B(3)=z0+arretecube/2.d0
         
         xx0=1.d0
         yy0=1.d0
         zz0=1.d0
         if (dabs((z-z0)).le.arretecubem) then
            zz0=0.d0
         endif
         if (dabs((x-x0)).le.arretecubem) then
            xx0=0.d0
         endif
         if (dabs((y-y0)).le.arretecubem) then
            yy0=0.d0
         endif
         call  DCUHRE(NDIM,NF,A,B, MINCLS, MAXCLS, fonctionigtadda,
     $        ABSREQ,RELREQ,KEY,NW,0,finest,ABSEST,NEVAL,IFAIL, WRKSTR) 
         do N = 1,NF
            FINEST(N)=FINEST(N)/arretecube/arretecube/arretecube
         enddo

         if (ifail.ne.0) then
            write(*,*) 'IFAIL in IGT routine',IFAIL
            stop
         endif

         
         propaint(1,1)=FINEST(1)+icomp*FINEST(7)
         propaint(1,2)=FINEST(2)+icomp*FINEST(8)
         propaint(1,3)=FINEST(3)+icomp*FINEST(9)
         propaint(2,2)=FINEST(4)+icomp*FINEST(10)
         propaint(2,3)=FINEST(5)+icomp*FINEST(11)
         propaint(3,3)=FINEST(6)+icomp*FINEST(12)
         propaint(2,1)=propaint(1,2)
         propaint(3,1)=propaint(1,3)
         propaint(3,2)=propaint(2,3)
      endif

      end
c*************************************************************
      subroutine fonctionigtadda(ndim,zz,nfun,f)
      implicit none
      integer ndim,nfun
      double precision zz(ndim),f(nfun)
      
      integer i,j
      double precision x,y,z,x0,y0,z0,k0,Id(3,3),Rab,Rtenseur(3,3)
     $     ,Rvect(3),xx0,yy0,zz0
      double complex propaesplibre(3,3),const1,const2
      common/k0xyz/k0,x,y,z,xx0,yy0,zz0
      x0=zz(1)
      y0=zz(2)
      z0=zz(3)

      Rab=0.d0
      Rvect(1)=(x-x0)
      Rvect(2)=(y-y0)
      Rvect(3)=(z-z0)

      do i=1,3
         do j=1,3
            Id(i,j)=0.d0
            if (i.eq.j) Id(i,i)=1.d0
            Rtenseur(i,j)=Rvect(i)*Rvect(j)
         enddo
         Rab=Rab+Rvect(i)*Rvect(i)
      enddo
      Rab=dsqrt(Rab)

c     normalise pour avoir le vecteur unitaire
      do i=1,3
         do j=1,3
            Rtenseur(i,j)=Rtenseur(i,j)/(Rab*Rab)
         enddo
      enddo
    
      const1=(Rab*(1.d0,0.d0))**(-3.d0)-(0.d0,1.d0)*k0*(Rab**(-2.d0))
      const2=k0*k0/Rab*(1.d0,0.d0)
      do i=1,3
         do j=1,3
            propaesplibre(i,j)=((3.d0*Rtenseur(i,j)-Id(i,j))*const1+
     *           (Id(i,j)-Rtenseur(i,j))*const2)*
     *           cdexp((0.d0,1.d0)*k0*Rab)
         enddo
      enddo

      f(1)=dreal(propaesplibre(1,1))
      f(2)=dreal(propaesplibre(1,2))*xx0*yy0
      f(3)=dreal(propaesplibre(1,3))*xx0*zz0
      f(4)=dreal(propaesplibre(2,2))
      f(5)=dreal(propaesplibre(2,3))*yy0*zz0
      f(6)=dreal(propaesplibre(3,3))

      f(7)=dimag(propaesplibre(1,1))
      f(8)=dimag(propaesplibre(1,2))*xx0*yy0
      f(9)=dimag(propaesplibre(1,3))*xx0*zz0
      f(10)=dimag(propaesplibre(2,2))
      f(11)=dimag(propaesplibre(2,3))*yy0*zz0
      f(12)=dimag(propaesplibre(3,3))
      end
