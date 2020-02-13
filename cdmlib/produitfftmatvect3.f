c     FFB=FFB-A*FFX==> FFB=XI et FFX=XI*pola==> FFB=(I-A*pola)*XI
      subroutine produitfftmatvect3(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $     ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz,Tabdip
     $     ,ntotalm,ntotal,nmax,ndipole,nxm,nym ,nzm ,nx,ny,nz,nx2,ny2
     $     ,nxy2,nz2,FFX,FFB,planb,planf)
      implicit none
      integer ii,jj,i,j,k,indice,ntotal,ntotalm,nxm,nym,nzm,nx,ny,nz
     $     ,nmax,nx2,ny2,nxy2,nz2,ndipole,Tabdip(nxm*nym*nzm)
      double complex FFTTENSORxx(ntotalm),FFTTENSORxy(ntotalm)
     $     ,FFTTENSORxz(ntotalm),FFTTENSORyy(ntotalm)
     $     ,FFTTENSORyz(ntotalm),FFTTENSORzz(ntotalm),FFX(3*nmax),FFB(3
     $     *nmax),vectx(ntotalm) ,vecty(ntotalm),vectz(ntotalm)
      double complex ctmpx,ctmpy,ctmpz
      integer n,nxy
      double precision dntotal
      integer*8 planf,planb
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD
c     double precision t1,t2

      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      nxy=nx*ny
      dntotal=dble(ntotal)

c     calcul FFT du vecteur B
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice)
!$OMP DO SCHEDULE(STATIC)    
      do indice=1,ntotal
         vectx(indice)=0.d0
         vecty(indice)=0.d0
         vectz(indice)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii,jj,i,j,k,indice) 
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(3)      
      do k=1,nz
         do j=1,ny
            do i=1,nx
c     position du dipole
               ii=i+nx*(j-1)+nxy*(k-1)
               if (Tabdip(ii).ne.0) then
                  indice=i+nx2*(j-1)+nxy2*(k-1)
                  jj=3*Tabdip(ii)
                  vectx(indice)=FFX(jj-2)
                  vecty(indice)=FFX(jj-1)
                  vectz(indice)=FFX(jj)
               endif

            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

#ifdef USE_FFTW
      call dfftw_execute_dft(planb, vectx, vectx)   
      call dfftw_execute_dft(planb, vecty, vecty)     
      call dfftw_execute_dft(planb, vectz, vectz)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz3d(vectx,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(vecty,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(vectz,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL           
#endif

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice,ctmpx,ctmpy,ctmpz)   
!$OMP DO  SCHEDULE(STATIC) 
      do indice=1,ntotal
         ctmpx=vectx(indice)
         ctmpy=vecty(indice)
         ctmpz=vectz(indice)
         vectx(indice)=FFTTENSORxx(indice)*ctmpx+FFTTENSORxy(indice)
     $        *ctmpy+ FFTTENSORxz(indice)*ctmpz
         vecty(indice)=FFTTENSORxy(indice)*ctmpx+FFTTENSORyy(indice)
     $        *ctmpy+ FFTTENSORyz(indice)*ctmpz         
         vectz(indice)=FFTTENSORxz(indice)*ctmpx+FFTTENSORyz(indice)
     $        *ctmpy+ FFTTENSORzz(indice)*ctmpz
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     FFT inverse (deconvolution)
#ifdef FFTW
      call dfftw_execute_dft(planf, vectx, vectx)
      call dfftw_execute_dft(planf, vecty, vecty)
      call dfftw_execute_dft(planf, vectz, vectz)
#endif

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice,jj,i,j,k,n) 
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)        
      do k=1,nz
         do j=1,ny
            do i=1,nx
               n=i+nx*(j-1)+nxy*(k-1)
               if (Tabdip(n).ne.0) then
                  indice=i+nx2*(j-1)+nxy2*(k-1)
                  jj=3*Tabdip(n)
                  FFB(jj-2)=FFB(jj-2)-vectx(indice)/dntotal
                  FFB(jj-1)=FFB(jj-1)-vecty(indice)/dntotal
                  FFB(jj)=FFB(jj)-vectz(indice)/dntotal
               endif             
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      end
   

c*************************************************
