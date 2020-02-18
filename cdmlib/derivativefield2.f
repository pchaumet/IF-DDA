c***************************************************************
c     Compute the derivative of the local field
c***************************************************************
      subroutine derivativefield2(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
     $     ,ntotal,ntotalm,nmax,DFFTTENSORxx,DFFTTENSORxy,DFFTTENSORxz
     $     ,DFFTTENSORyy,DFFTTENSORyz,DFFTTENSORzz,vectx,vecty,vectz
     $     ,test,FF0,planb,planf)
      implicit none
      integer i,j,k,ii,jj,kk,nx,ny,nz,nxy2,nx2,ny2,nz2,indice,nmax
     $     ,ntotal,ntotalm,test
      double precision k0,x0,y0,z0,xx0,yy0,zz0,aretecube,dntotal
      double complex DFFTTENSORxx(ntotalm),DFFTTENSORxy(ntotalm)
     $     ,DFFTTENSORxz(ntotalm),DFFTTENSORyy(ntotalm)
     $     ,DFFTTENSORyz(ntotalm),DFFTTENSORzz(ntotalm),vectx(ntotalm)
     $     ,vecty(ntotalm),vectz(ntotalm),FF0(3*nmax),Stenseurd(3,3,3)
      double complex ctmpx,ctmpy,ctmpz
      integer*8 planb,planf
      integer FFTW_BACKWARD,FFTW_FORWARD
      FFTW_BACKWARD=+1
      FFTW_FORWARD=-1      
      
      xx0=0.d0
      yy0=0.d0
      zz0=0.d0
      dntotal=dble(ntotal)
      
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(ii,jj,kk,indice,i,j,k,x0,y0,z0,Stenseurd)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)    
      do kk=1,nz2
         do jj=1,ny2
            do ii=1,nx2
               indice=ii+nx2*(jj-1)+nxy2*(kk-1)
               if (ii.eq.nx+1.or.jj.eq.ny+1.or.kk.eq.nz+1) then
                  DFFTTENSORxx(indice)=0.d0
                  DFFTTENSORxy(indice)=0.d0
                  DFFTTENSORxz(indice)=0.d0
                  DFFTTENSORyy(indice)=0.d0
                  DFFTTENSORyz(indice)=0.d0
                  DFFTTENSORzz(indice)=0.d0
               else
                  if (ii.gt.nx) then
                     i=(ii-1)-nx2
                  else
                     i=ii-1
                  endif
                  if (jj.gt.ny) then
                     j=(jj-1)-ny2
                  else
                     j=jj-1
                  endif
                  if (kk.gt.nz) then
                     k=(kk-1)-nz2
                  else
                     k=kk-1
                  endif
                  x0=dble(i)*aretecube
                  y0=dble(j)*aretecube
                  z0=dble(k)*aretecube
                  call propesplibdermat(x0,y0,z0,xx0,yy0,zz0,k0,test
     $                 ,Stenseurd)
                  DFFTTENSORxx(indice)=Stenseurd(1,1,test)
                  DFFTTENSORxy(indice)=Stenseurd(1,2,test)
                  DFFTTENSORxz(indice)=Stenseurd(1,3,test)
                  DFFTTENSORyy(indice)=Stenseurd(2,2,test)
                  DFFTTENSORyz(indice)=Stenseurd(2,3,test)
                  DFFTTENSORzz(indice)=Stenseurd(3,3,test)
               endif
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c     Compute FFT

#ifdef USE_FFTW
      call dfftw_execute_dft(planb,DFFTTENSORxx,DFFTTENSORxx)
      call dfftw_execute_dft(planb,DFFTTENSORxy,DFFTTENSORxy) 
      call dfftw_execute_dft(planb,DFFTTENSORxz,DFFTTENSORxz)
      call dfftw_execute_dft(planb,DFFTTENSORyy,DFFTTENSORyy)
      call dfftw_execute_dft(planb,DFFTTENSORyz,DFFTTENSORyz)
      call dfftw_execute_dft(planb,DFFTTENSORzz,DFFTTENSORzz)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORxx,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORxy,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORxz,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORyy,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORyz,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORzz,NX2,NY2,NZ2,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif

c     product of the FFT
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice,ctmpx,ctmpy,ctmpz)
!$OMP DO SCHEDULE(STATIC) 
      do indice=1,ntotal
         ctmpx=vectx(indice)
         ctmpy=vecty(indice)
         ctmpz=vectz(indice)
         DFFTTENSORxx(indice)=DFFTTENSORxx(indice)*ctmpx
     $        +DFFTTENSORxy(indice)*ctmpy+DFFTTENSORxz(indice) *ctmpz
         
         DFFTTENSORyy(indice)=DFFTTENSORxy(indice)*ctmpx
     $        +DFFTTENSORyy(indice)*ctmpy+DFFTTENSORyz(indice) *ctmpz
         
         DFFTTENSORzz(indice)=DFFTTENSORxz(indice)*ctmpx
     $        +DFFTTENSORyz(indice)*ctmpy+DFFTTENSORzz(indice) *ctmpz
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL        

c     FFT inverse (deconvolution)

#ifdef USE_FFTW
      call dfftw_execute_dft(planf,DFFTTENSORxx,DFFTTENSORxx)
      call dfftw_execute_dft(planf,DFFTTENSORyy,DFFTTENSORyy)
      call dfftw_execute_dft(planf,DFFTTENSORzz,DFFTTENSORzz)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORxx,NX2,NY2,NZ2,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORyy,NX2,NY2,NZ2,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz3d(DFFTTENSORzz,NX2,NY2,NZ2,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL

#endif

c     remet le vecteur resultat dans FF0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,indice,ii)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)+nxy2*(k-1)
               ii=3*(i+nx*(j-1)+nx*ny*(k-1))
               FF0(ii-2)=DFFTTENSORxx(indice)/dntotal
               FF0(ii-1)=DFFTTENSORyy(indice)/dntotal
               FF0(ii)=DFFTTENSORzz(indice)/dntotal
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
      end      
