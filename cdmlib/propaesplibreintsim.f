c     Integration of the Green Tensor
      subroutine propaespacelibreintsim(x,y,z,x0,y0,z0,k0,arretecube
     $     ,lim,RELREQ,propaint)
      implicit none
      integer i,j,k,nstep
c     definition of the position of the dipole, observation, wavenumber
c     ,wavelength, spacing lattice
      double precision xa,ya,za,pi,arretecubem,lim
      double precision x,y,z,x0,y0,z0,arretecube,k0,xx0,yy0,zz0,lambda

c     Green tensor and complex number i,1 and constant
      double complex propaint(3,3),icomp,uncomp,prop(3,3)  
      double precision Rab,tmp,RELREQ,A(3),w(0:128),ww,pas
      

      
    
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)

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
         call propa_espace_libre(x,y,z,x0,y0,z0,k0,propaint)
         
      else
c     We perform the integration of the tensor
c     defintion for the integration
         nstep=128
         call propa_espace_libre(x,y,z,x0,y0,z0,k0,prop)
c         write(*,*) 'prop',prop
         A(1)=x-arretecube/2.d0
         A(2)=y-arretecube/2.d0
         A(3)=z-arretecube/2.d0
         
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
         w(0)=1.d0
         w(nstep)=1.d0
         do i=1,nstep/2-1
            w(2*i)=2.d0
         enddo
         do i=1,nstep/2
            w(2*i-1)=4.d0
         enddo
         pas=arretecube/dble(nstep)

         propaint=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,xa,ya,za,ww,prop)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:propaint)
         do i=0,nstep
            do j=0,nstep
               do k=0,nstep
                  xa=A(1)+pas*dble(i)
                  ya=A(2)+pas*dble(j)
                  za=A(3)+pas*dble(k)
                  ww=w(i)*w(j)*w(k)
                  call propa_espace_libre(xa,ya,za,x0,y0,z0,k0,prop)
                  propaint(1,1)=propaint(1,1)+prop(1,1)*ww
                  propaint(1,2)=propaint(1,2)+prop(1,2)*ww
                  propaint(1,3)=propaint(1,3)+prop(1,3)*ww
                  propaint(2,2)=propaint(2,2)+prop(2,2)*ww
                  propaint(2,3)=propaint(2,3)+prop(2,3)*ww
                  propaint(3,3)=propaint(3,3)+prop(3,3)*ww
               enddo
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         tmp=pas/3.d0*pas/3.d0*pas/3.d0/arretecube/arretecube/arretecube

         propaint(1,1)=propaint(1,1)*tmp
         propaint(1,2)=propaint(1,2)*tmp*xx0*yy0
         propaint(1,3)=propaint(1,3)*tmp*xx0*zz0
         propaint(2,2)=propaint(2,2)*tmp
         propaint(2,3)=propaint(2,3)*tmp*yy0*zz0
         propaint(3,3)=propaint(3,3)*tmp
         propaint(2,1)=propaint(1,2)
         propaint(3,1)=propaint(1,3)
         propaint(3,2)=propaint(2,3)
c         write(*,*) 'prop',propaint
      endif

      end
