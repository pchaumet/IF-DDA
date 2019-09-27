      subroutine ondecircekxky(x,y,z,k0,E0,s,k0x,k0y,Ex,Ey,Ez)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
      implicit none
      double precision x,y,z,k0,s,theta,phi,pi,thetat,phit,k0x,k0y,k0z
     $     ,us,nsx,nsy,npx,npy,npz,np
      double complex E0,Ex,Ey,Ez,icomp,uncomp,E0s,E0p,E0x,E0y ,E0z
     $     ,exparg,Ex1,Ey1,Ez1,Ex2,Ey2,Ez2

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      us=1.d0/dsqrt(2.d0)

      E0s=E0*us
      E0p=E0*icomp*dsign(us,s)

      k0z=dsqrt(k0*k0-k0x*k0x-k0y*k0y)


      if (k0z.ge.0.99999d0*k0) then

         E0x=E0p
         E0y=E0s
         E0z=0.d0

      else
         nsx=k0y/dsqrt(k0x*k0x+k0y*k0y)
         nsy=-k0x/dsqrt(k0x*k0x+k0y*k0y)

         npx=-k0z*nsy
         npy=k0z*nsx
         npz=k0x*nsy-k0y*nsx

         np=dsqrt(npx*npx+npy*npy+npz*npz)
         npx=npx/np
         npy=npy/np
         npz=npz/np
         
         Ex=E0s*nsx+E0p*npx
         Ey=E0s*nsy+E0p*npy
         Ez=E0p*npz

      endif




      end
 
