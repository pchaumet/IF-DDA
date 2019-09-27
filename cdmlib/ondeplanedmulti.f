      subroutine ondeplanedmulti(x,y,z,k0,E0m,ssm,ppm,thetam,phim,nbinc
     $     ,test,Edert,nstop,infostr)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
      implicit none
      integer nstop,nbinc,i,test
      double precision x,y,z,k0,ssm(10),ppm(10),ss,pp,thetam(10)
     $     ,phim(10),theta,phi
      double complex E0m(10),E0,Eder(3,3),Edert(3 ,3)

      character(64) infostr

      if (nbinc.gt.10) then
         infostr='too many plane wave'
         nstop=1
         return
      endif
      Edert=0.d0
      
      do i=1,nbinc

         theta=thetam(i)
         phi=phim(i)
         ss=ssm(i)
         pp=ppm(i)
         E0=E0m(i)
         call ondeplaned(x,y,z,k0,E0,ss,pp,theta,phi,test,Eder)
         if (test.ne.4) then
            Edert(1,test)=Edert(1,test)+Eder(1,test)
            Edert(2,test)=Edert(2,test)+Eder(2,test)
            Edert(3,test)=Edert(3,test)+Eder(3,test)
         elseif(test.eq.4) then
            Edert(1,1)=Edert(1,1)+Eder(1,1)
            Edert(2,1)=Edert(2,1)+Eder(2,1)
            Edert(3,1)=Edert(3,1)+Eder(3,1)
            Edert(1,2)=Edert(1,2)+Eder(1,2)
            Edert(2,2)=Edert(2,2)+Eder(2,2)
            Edert(3,2)=Edert(3,2)+Eder(3,2)
            Edert(1,3)=Edert(1,3)+Eder(1,3)
            Edert(2,3)=Edert(2,3)+Eder(2,3)
            Edert(3,3)=Edert(3,3)+Eder(3,3)
         else
            infostr='wrong value for test'
            nstop=1
            return  
         endif
            
      enddo
      
      end
