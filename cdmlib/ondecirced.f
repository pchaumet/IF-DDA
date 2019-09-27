      subroutine ondecirced(x,y,z,k0,E0,s,thetat,phit,test,Eder)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
c     si test =1 ==>x,test =2 ==>y,test =3 ==>z,test =4 ==>tous
      implicit none
      integer test,ii,jj
      double precision x,y,z,k0,s,theta,phi,pi,thetat,phit,k0x,k0y
     $     ,k0z,us
      double complex E0,icomp,uncomp,E0s,E0p,E0x,E0y,E0z,exparg
      double complex Eder(3,3),expargi,Eder1(3,3),Eder2(3,3)

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      us=1.d0/dsqrt(2.d0)
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

      E0s=E0/dsqrt(2.d0)
      E0p=E0*icomp*dsign(us,s)
      k0x=k0*dsin(theta)*dcos(phi)
      k0y=k0*dsin(theta)*dsin(phi)
      k0z=k0*dcos(theta)

      exparg=cdexp(icomp*(x*k0x+y*k0y+z*k0z))
      expargi=exparg*icomp


c     polarisation s
      E0y=E0s*dcos(phi)
      E0x=-E0s*dsin(phi)
      E0z=0        



      if (test.eq.1) then
         
         Eder1(1,1)=E0x*expargi*k0x
         Eder1(2,1)=E0y*expargi*k0x
         Eder1(3,1)=E0z*expargi*k0x


      elseif(test.eq.2) then

         Eder1(1,2)=E0x*expargi*k0y
         Eder1(2,2)=E0y*expargi*k0y
         Eder1(3,2)=E0z*expargi*k0y


      elseif(test.eq.3) then

         Eder1(1,3)=E0x*expargi*k0z
         Eder1(2,3)=E0y*expargi*k0z
         Eder1(3,3)=E0z*expargi*k0z        


      elseif(test.eq.4) then

         Eder1(1,1)=E0x*expargi*k0x
         Eder1(2,1)=E0y*expargi*k0x
         Eder1(3,1)=E0z*expargi*k0x
         Eder1(1,2)=E0x*expargi*k0y
         Eder1(2,2)=E0y*expargi*k0y
         Eder1(3,2)=E0z*expargi*k0y
         Eder1(1,3)=E0x*expargi*k0z
         Eder1(2,3)=E0y*expargi*k0z
         Eder1(3,3)=E0z*expargi*k0z   


      endif

c     polarisation p
      E0y=E0p*dcos(theta)*dsin(phi)
      E0x=E0p*dcos(theta)*dcos(phi)
      E0z=-E0p*dsin(theta)

      if (test.eq.1) then
         
         Eder2(1,1)=E0x*expargi*k0x
         Eder2(2,1)=E0y*expargi*k0x
         Eder2(3,1)=E0z*expargi*k0x

      elseif(test.eq.2) then

         Eder2(1,2)=E0x*expargi*k0y
         Eder2(2,2)=E0y*expargi*k0y
         Eder2(3,2)=E0z*expargi*k0y


      elseif(test.eq.3) then

         Eder2(1,3)=E0x*expargi*k0z
         Eder2(2,3)=E0y*expargi*k0z
         Eder2(3,3)=E0z*expargi*k0z        

      elseif(test.eq.4) then

         Eder2(1,1)=E0x*expargi*k0x
         Eder2(2,1)=E0y*expargi*k0x
         Eder2(3,1)=E0z*expargi*k0x
         Eder2(1,2)=E0x*expargi*k0y
         Eder2(2,2)=E0y*expargi*k0y
         Eder2(3,2)=E0z*expargi*k0y
         Eder2(1,3)=E0x*expargi*k0z
         Eder2(2,3)=E0y*expargi*k0z
         Eder2(3,3)=E0z*expargi*k0z   


      endif

c     somme des deux polas
      do ii=1,3
         do jj=1,3
            Eder(ii,jj)=Eder1(ii,jj)+Eder2(ii,jj)
         enddo
      enddo
      return
      end
