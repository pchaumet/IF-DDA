      subroutine ondeplanemulti(x,y,z,k0,E0m,ssm,ppm,thetam,phim,nbinc
     $     ,Ext,Eyt,Ezt,nstop,infostr)
c     programme d'une onde plane
c     l axe z vertical sert de reference
c     pour la polarisation le plan (x,y) sert de reference
c     par default, phi=0,theta>0, k positif pour z croissant et x croissant
c     phi>0 ky>0
      implicit none
      integer nstop,nbinc,i
      double precision x,y,z,k0,ssm(10),ppm(10),ss,pp,thetam(10)
     $     ,phim(10),theta,phi
      double complex E0m(10),Ex,Ey,Ez,Ext,Eyt,Ezt,E0

      character(64) infostr

      if (nbinc.gt.10) then
         infostr='too many plane wave'
         nstop=1
         return
      endif
      Ext=0.d0
      Eyt=0.d0
      Ezt=0.d0

      do i=1,nbinc

         theta=thetam(i)
         phi=phim(i)
         ss=ssm(i)
         pp=ppm(i)
         E0=E0m(i)
         call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi,Ex,Ey,Ez,nstop
     $        ,infostr)
         Ext=Ext+Ex
         Eyt=Eyt+Ey
         Ezt=Ezt+Ez
         
      enddo
      
      end
