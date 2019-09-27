      subroutine dipoleincder(xdip,ydip,zdip,thetat,phit,x,y,z,aretecube
     $     ,k0,E0,Eder,test,nstop,infostr)
      implicit none
      integer nstop,test
      double precision x,y,z,k0,pi,xdip,ydip,zdip,theta,phi,thetat,phit
     $     ,aretecube,dist,vol
      double complex E0,Eder(3,3),icomp,uncomp,Stenseurd(3,3,3),p(3)

      character(64) infostr
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      vol=1.d0/(4.d0*pi)
      
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

c     calcul orientation et amplitude du dipole
      
      p(1)=E0*dsin(theta)*dcos(phi)
      p(2)=E0*dsin(theta)*dsin(phi)
      p(3)=E0*dcos(theta)
      
      dist=dsqrt((xdip-x)*(xdip-x)+(ydip-y)*(ydip-y) +(zdip-z)*(zdip-z))

      if (dist.le.aretecube/100.d0) then
         Eder=0.d0
      else
       
         call propesplibdermat(x,y,z,xdip,ydip,zdip,k0,test,Stenseurd)
         Eder(1,test)=(Stenseurd(1,1,test)*p(1)+Stenseurd(1,2,test)
     $        *p(2)+Stenseurd(1,3,test)*p(3))*vol
         Eder(2,test)=(Stenseurd(2,1,test)*p(1)+Stenseurd(2,2,test)
     $        *p(2)+Stenseurd(2,3,test)*p(3))*vol
         Eder(3,test)=(Stenseurd(3,1,test)*p(1)+Stenseurd(3,2,test)
     $        *p(2)+Stenseurd(3,3,test)*p(3))*vol

      endif


      end
