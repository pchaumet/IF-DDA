      subroutine ondeplane(x,y,z,k0,E0,ss,pp,thetat,phit,Ex,Ey,Ez,nstop
     $     ,infostr)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
      implicit none
      integer nstop
      double precision x,y,z,k0,ss,pp,s,p,theta,phi,pi,thetat,phit,k0x
     $     ,k0y,k0z
      double complex E0,Ex,Ey,Ez,icomp,uncomp,E0s,E0p,E0x,E0y,E0z,exparg

      character(64) infostr

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)

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

      E0s=E0*s*dsign(1.d0,theta)
      E0p=E0*p*dsign(1.d0,theta)
      k0x=k0*dsin(theta)*dcos(phi)
      k0y=k0*dsin(theta)*dsin(phi)
      k0z=k0*dcos(theta)

      E0y=E0s*dcos(phi)+E0p*dcos(theta)*dsin(phi)
      E0x=-E0s*dsin(phi)+E0p*dcos(theta)*dcos(phi)     
      E0z=-E0p*dsin(theta)

      exparg=cdexp(icomp*(x*k0x+y*k0y+z*k0z))

      Ex=E0x*exparg
      Ey=E0y*exparg
      Ez=E0z*exparg
      end
