      subroutine ondeplanekxky(x,y,z,k0,E0,ss,pp,k0x,k0y,Ex,Ey,Ez,nstop
     $     ,infostr)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     ns=kinc vect z et np=kinc vect ns base (kinc,ns,np)
c     par defaut l'onde se propage selon z positif
      implicit none
      integer nstop
      double precision x,y,z,k0,ss,pp,s,p,k0x ,k0y,k0z,nsx,nsy,npx,npy
     $     ,npz,np
      double complex E0,Ex,Ey,Ez,icomp,uncomp,E0s,E0p,E0x,E0y,E0z,exparg

      character(64) infostr

      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)

      if (pp.lt.0.d0.or.pp.gt.1.d0) then
         nstop=-1
         infostr='problem in plance incident wave'
         return
      endif

      p=dsqrt(pp)
      s=dsqrt(1.d0-pp)

      E0s=E0*s
      E0p=E0*p

      k0z=dsqrt(k0*k0-k0x*k0x-k0y*k0y)
      if (k0z.ge.0.99999d0*k0) then

         E0x=E0p
         E0y=E0s
         E0z=0.d0

      else
         nsx=-k0y/dsqrt(k0x*k0x+k0y*k0y)
         nsy=k0x/dsqrt(k0x*k0x+k0y*k0y)
         npx=k0z*nsy
         npy=-k0z*nsx
         npz=-k0x*nsy+k0y*nsx
         np=dsqrt(npx*npx+npy*npy+npz*npz)
         npx=npx/np
         npy=npy/np
         npz=npz/np
         
         E0x=E0s*nsx+E0p*npx
         E0y=E0s*nsy+E0p*npy
         E0z=E0p*npz

      endif
      
      exparg=cdexp(icomp*(x*k0x+y*k0y+z*k0z))

      Ex=E0x*exparg
      Ey=E0y*exparg
      Ez=E0z*exparg
      end
