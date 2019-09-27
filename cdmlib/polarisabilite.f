      subroutine poladiff(rayon1,eps,eps0,k0,dddis,
     *     methode,polarisabilite)
      implicit none
      integer dddis
      double precision rayon,volume,pi,k0,k03,eps0,SLDR,rayon1,indice
     $     ,maille
      double complex eps,rapport,pola,icomp
      double complex polarisabilite
      character*2 methode

c     dis: 0 si c est une sphere et 1 c est un element de discretisation
c     rayon=rayonon si dis=0 sinon facteur pour le cube -->sphere
c     methode: LA:lakthakia; CM: claussius-Mossotti;
c     RR: CM avec reaction de rayononnement;DB: dungey borhen;
c     DR: draine;DD: Draine 2

c      write(*,*) rayon,eps,eps0,k0,methode,dis

      pi=dacos(-1.d0)
      rapport=(eps-eps0)/(eps+2.d0*eps0)
      icomp=(0.d0,1.d0)
      indice=dsqrt(eps0)
      k03=k0*k0*k0*indice*indice*indice
      rayon=rayon1
      
      
      if (dddis.eq.0) then
         volume=rayon*rayon*rayon
      elseif (dddis.eq.1) then
         maille=rayon
         rayon=rayon*((0.75d0/pi)**(1.d0/3.d0))
         volume=rayon*rayon*rayon
      else
         write(*,*) 'mauvaise valeur de dddis',dddis
         stop
      endif

c     polarisabilite de CM
      pola=rapport*volume

      if (methode.eq.'CM') then
         polarisabilite=pola*eps0
      elseif (methode.eq.'RR'.or.methode.eq.'PS') then
         polarisabilite=eps0*pola/(1.d0-2.d0/3.d0*icomp*k03*pola)     
      elseif (methode.eq.'GB') then
         polarisabilite=pola/(1.d0-2.d0/3.d0*icomp*k03*pola-k0*k0*pola
     $        *indice*indice/rayon)
      elseif (methode.eq.'LA') then         
         polarisabilite=pola/(1.d0-2.d0*rapport*((1.d0-icomp*k0*indice
     $        *rayon)*cdexp(icomp*k0*indice*rayon)-1.d0))        
      elseif (methode.eq.'LR') then
         SLDR=1.d0/5.d0
         polarisabilite=eps0*pola/( 1.d0+pola*( (-1.8915316d0+eps/eps0
     $        *0.1648469d0-eps/eps0*1.7700004d0*SLDR)*k0*k0*indice
     $        *indice/maille-2.d0/3.d0*icomp*k03*indice*indice*indice) )
      else
         write(*,*) 'incorrect method'
         write(*,*) 'method',methode
         stop
      endif

      end
c********************************************************************
c********************************************************************
c********************************************************************
c     calcul de la pola par integration sur la maille
      subroutine intmaille(aretecube,k0,polarisabilite)
      implicit none

      integer  KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
      parameter (nw=40000000,ndim=2,nf=2)
      double precision A(NDIM), B(NDIM), WRKSTR(NW) 
      double precision  ABSEST(NF), FINEST(NF), ABSREQ, RELREQ,pi

      integer i
      double precision k0,aretecube,aa,polar  
      double complex icomp,polam,polarisabilite

      external integrationmaillecentrale1
      external integrationmaillecentrale2

      common/k0aretecube/aa

      aa=k0*aretecube
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)

c      write(*,*) 'eee',k0,aretecube,aa

c     mode propagatif
      MINCLS = 4000
      MAXCLS = 10000000
      KEY = 0
      ABSREQ = 0.d0
      RELREQ = 1.d-10

      A(1)=0.d0
      A(2)=0.d0
      B(1)=1.d0
      B(2)=pi/2.d0

      call  DCUHRE(NDIM,NF,A,B, MINCLS, MAXCLS
     $     ,integrationmaillecentrale1, ABSREQ, RELREQ,KEY,NW,0,finest
     $     ,ABSEST,NEVAL,IFAIL, WRKSTR) 
      
      do N = 1,NF         
c         write(*,*) N, ABSEST(N), FINEST(N) ,IFAIL
         FINEST(N)=FINEST(N)
         if (IFAIL.ne.0) stop
      enddo
      
      polam=16/pi*(FINEST(1)+icomp*FINEST(2))
c      write(*,*) 'polam',polam,2.d0/3.d0*k0*k0*k0*aretecube*aretecube
c     $     *aretecube
c     mode evanescent
      polar=0.d0

      do i=0,1000000
         MINCLS = 10000
         MAXCLS = 40000000
         KEY = 0
         ABSREQ = 0.d0
         RELREQ = 1.d-8
         
         A(1)=100.d0*dble(i)
         A(2)=0.d0
         B(1)=100.d0*dble(i+1)
         B(2)=pi/2.d0
      
         call  DCUHRE(NDIM,NF,A,B, MINCLS, MAXCLS
     $        ,integrationmaillecentrale2, ABSREQ, RELREQ,KEY,NW,0
     $        ,finest ,ABSEST,NEVAL,IFAIL, WRKSTR) 
         
         do N = 1,NF         
c            write(*,*) N, ABSEST(N), FINEST(N) ,IFAIL
            FINEST(N)=FINEST(N)
            if (IFAIL.ne.0) then
               write(*,*) 'IFAIL',IFAIL,NEVAL
            endif
         enddo

c         write(*,*) 'essai',16/pi*FINEST(1),polar*RELREQ/10.d0
         if (i.ne.0.and.dabs(16/pi*FINEST(1)).le.dabs(polar*RELREQ/10.d0
     $        )) goto 10
         polar=polar+16/pi*FINEST(1)
      enddo

 10   polarisabilite=polar+polam+4.d0/3.d0*pi

      end
c*************************************************
      subroutine integrationmaillecentrale1(ndim,zz,nfun,f)
      implicit none
      integer ndim,nfun
      double precision zz(ndim),f(nfun)
      
      double precision aa,pi,racine,sina,cosa,x,alpha,u
      double complex ff,expixu,icomp

      common/k0aretecube/aa
      
      x=zz(1)
      alpha=zz(2)
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      racine=dsqrt(1.d0-x*x)
      u=aa/2.d0
      sina=dsin(alpha)
      cosa=dcos(alpha)
      expixu=cdexp(icomp*x*u)

      if (x.eq.0.d0) then
         ff=icomp*u
      else
         ff=(expixu-1.d0)/x-x*expixu
      endif

      if (x.eq.1.d0) then
         ff=ff*u*u
      else
         if (alpha.eq.0.d0) then
            ff=ff*u*dsin(racine*u*cosa)/cosa/racine
         elseif (alpha.eq.pi/2.d0) then
            ff=ff*u*dsin(racine*u*sina)/sina/racine
         else
            ff=ff*dsin(racine*u*sina)/sina*dsin(racine*u*cosa)/cosa
     $           /racine/racine
         endif
      endif

      f(2)=dimag(ff)
      f(1)=dreal(ff)

      end
c*************************************************
      subroutine integrationmaillecentrale2(ndim,zz,nfun,f)
      implicit none
      integer ndim,nfun
      double precision zz(ndim),f(nfun)
      
      double precision aa,pi,racine,sina,cosa,x,alpha,u,ff
      double precision fac
      common/k0aretecube/aa
      
      x=zz(1)
      alpha=zz(2)
      pi=dacos(-1.d0)
      fac=1.d0+x*x
      racine=dsqrt(fac)
      u=aa/2.d0
      sina=dsin(alpha)
      cosa=dcos(alpha)

      if (x.eq.0.d0) then
         ff=u
      else
         ff=(1.d0-fac*dexp(-x*u))/x/fac
      endif

      if (alpha.eq.0.d0) then
         ff=ff*racine*u*dsin(racine*u*cosa)/cosa
      elseif (alpha.eq.pi/2.d0) then
         ff=ff*racine*u*dsin(racine*u*sina)/sina
      else
         ff=ff*dsin(racine*u*sina)/sina*dsin(racine*u*cosa)/cosa
      endif

      f(1)=ff
      f(2)=0.d0

      end

