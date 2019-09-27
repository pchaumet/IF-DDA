      subroutine ondecirce(x,y,z,k0,E0,s,thetat,phit,Ex,Ey,Ez)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
      implicit none
      double precision x,y,z,k0,s,theta,phi,pi,thetat,phit,k0x,k0y,k0z
     $     ,us
      double complex E0,Ex,Ey,Ez,icomp,uncomp,E0s,E0p,E0x,E0y ,E0z
     $     ,exparg,Ex1,Ey1,Ez1,Ex2,Ey2,Ez2

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

c     polarisation s
      E0y=E0s*dcos(phi)
      E0x=-E0s*dsin(phi)
      E0z=0.d0
        
      Ex1=E0x*exparg
      Ey1=E0y*exparg
      Ez1=E0z*exparg

c     polarisation p
      E0y=E0p*dcos(theta)*dsin(phi)
      E0x=E0p*dcos(theta)*dcos(phi)
      E0z=-E0p*dsin(theta)

      Ex2=E0x*exparg
      Ey2=E0y*exparg
      Ez2=E0z*exparg

      Ex=Ex1+Ex2
      Ey=Ey1+Ey2
      Ez=Ez1+Ez2

      end
 
