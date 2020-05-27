      subroutine dipoleinc(xdip,ydip,zdip,thetat,phit,x,y,z,aretecube,k0
     $     ,E0,Ex,Ey,Ez,nstop,infostr)
      implicit none
      integer nstop
      double precision x,y,z,k0,pi,xdip,ydip,zdip,theta,phi,thetat,phit
     $     ,aretecube,dist,vol
      double complex E0,Ex,Ey,Ez,icomp,uncomp,propaesplibre(3,3),p(3)
     $     ,const

      character(64) infostr
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      
      theta=thetat
      theta=theta*pi/180.d0
      phi=phit
      phi=phi*pi/180.d0

c     calcul orientation et amplitude du dipole
      
      p(1)=E0*dsin(theta)*dcos(phi)
      p(2)=E0*dsin(theta)*dsin(phi)
      p(3)=E0*dcos(theta)

    
      
      dist=dsqrt((xdip-x)*(xdip-x)+(ydip-y)*(ydip-y) +(zdip-z)*(zdip-z))
      vol=1/(4.d0*pi)

      if (dist.le.aretecube/100.d0) then
c     vol=aretecube*aretecube*aretecube/(4.d0*pi) mais pour ne pas
c     dependre de la maille pour l'amplitude je divise par le volume du
c     cube
       

         const=-1.d0/3.d0/(aretecube*aretecube*aretecube)+icomp*2.d0
     $        /3.d0*vol*k0*k0*k0

         Ex=const*p(1)
         Ey=const*p(2)
         Ez=const*p(3)
      else
      
         call propa_espace_libre(x,y,z,xdip,ydip,zdip,k0, propaesplibre)
    
         Ex=(propaesplibre(1,1)*p(1)+propaesplibre(1,2)*p(2)
     $        +propaesplibre(1,3)*p(3))*vol

         Ey=(propaesplibre(2,1)*p(1)+propaesplibre(2,2)*p(2)
     $        +propaesplibre(2,3)*p(3))*vol

         Ez=(propaesplibre(3,1)*p(1)+propaesplibre(3,2)*p(2)
     $        +propaesplibre(3,3)*p(3))*vol
      endif
      
      end
c****************************************************
c****************************************************
c****************************************************
      subroutine dipoleinside(xdip,ydip,zdip,xs,ys,zs,aretecube,nmax
     $     ,nbsphere)
      implicit none
      integer i,nmax,nbsphere,ic,jc,kc

      double precision xdip,ydip,zdip,aretecube
      double precision xs(nmax),ys(nmax),zs(nmax)


c     test si le dipole eclairant est au moins a une demi maille de tous
c     les dipoles (dans ce cas exterieur a l'objet), si proche alors mis
c     a la position du dipole le plus proche.


      do i=1,nbsphere
         ic=nint((xdip-xs(i))/aretecube)
         jc=nint((ydip-ys(i))/aretecube)
         kc=nint((zdip-zs(i))/aretecube)
        
         
         if (ic.eq.0.and.jc.eq.0.and.kc.eq.0) then
            xdip=xs(i)
            ydip=ys(i)
            zdip=zs(i)
            write(*,*) 'antenna inside the near field domain',xdip,ydip
     $           ,zdip
            return
         endif
      enddo

      write(*,*) 'antenna outside the near field domain'
      
      end
