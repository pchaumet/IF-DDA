      subroutine diffractefft2dlens(nx,ny,nz,nxm,nym,nzm,nfft2d,tabfft2
     $     ,k0,xs,ys,zs,aretecube,Eloinx,Eloiny,Eloinz,FF,imax,deltakx
     $     ,deltaky,Ediffkzpos,numaper,nside,plan2f,plan2b,nstop
     $     ,infostr)
      implicit none
      integer nx,ny,nz,nxm,nym,nzm,nfft2d,nstop,nside
      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
     $     ,aretecube,k0,numaper,dside
      double complex FF(3*nxm*nym*nzm),Ediffkzpos(nfft2d,nfft2d,3)


      integer nfft2d2,imax,i,j,k,tabfft2(nfft2d),indice,kk,ii,jj
      double precision deltakx,deltaky,var1,var2,kx,ky,kz,fac,pi
      double complex ctmp,ctmp1,icomp,Eloinx(nfft2d*nfft2d)
     $     ,Eloiny(nfft2d*nfft2d),Eloinz(nfft2d*nfft2d)
      character(64) infostr
      integer*8 plan2f,plan2b
      integer FFTW_FORWARD

      FFTW_FORWARD=-1   
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      dside=dble(nside)
      if (nfft2d.gt.16384) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif

      if (deltakx.ge.numaper) then
         nstop=1
         infostr='In FFT lens nfft2d too small'
         return
      endif

      nfft2d2=nfft2d/2  
      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltakx
      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltaky
c     fac=dble(nfft2d*nfft2d)
      fac=1.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
      do i=1,nfft2d
         if (i-nfft2d2-1.ge.0) then
            tabfft2(i)=i-nfft2d2
         else
            tabfft2(i)=nfft2d2+i
         endif
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      

      imax=nint(k0/deltakx)+1
      write(*,*) 'Number of point in NA : ',2*imax+1
      if (2*imax+1.gt.nfft2d) then
         write(*,*) 'Size',2*imax+1,nfft2d
         infostr='In FFT diffraction for lens nfft2d too small!'
         nstop = 1
         return
      endif
      if (2*imax+1.lt.7) then
         write(*,*) 'Size',2*imax+1,nfft2d
         infostr='In FFT diffract nfft2d too small'
         nstop = 1
         return
      endif

      
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i,j) 
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(2)     
      do i=1,nfft2d
         do j=1,nfft2d
            Ediffkzpos(i,j,1)=0.d0
            Ediffkzpos(i,j,2)=0.d0
            Ediffkzpos(i,j,3)=0.d0         
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

    
      do k=1,nz

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)  
!$OMP DO SCHEDULE(STATIC)           
         do i=1,nfft2d*nfft2d
            Eloinx(i)=0.d0
            Eloiny(i)=0.d0
            Eloinz(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kk,indice)  
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(2)     
         do i=1,nx
            do j=1,ny
               kk=i+nx*(j-1)+nx*ny*(k-1)
               indice=tabfft2(i)+nfft2d*(tabfft2(j)-1)
               Eloinx(indice)=FF(3*(kk-1)+1)                
               Eloiny(indice)=FF(3*(kk-1)+2)    
               Eloinz(indice)=FF(3*(kk-1)+3)  
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

#ifdef USE_FFTW
         call dfftw_execute_dft(plan2f,Eloinx,Eloinx)
         call dfftw_execute_dft(plan2f,Eloiny,Eloiny)
         call dfftw_execute_dft(plan2f,Eloinz,Eloinz)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Eloinx,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Eloiny,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Eloinz,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif
         
         kk=1+nx*ny*(k-1)

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp,ctmp1)
!$OMP DO  SCHEDULE(DYNAMIC) COLLAPSE(2)              
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               if (kx*kx+ky*ky.lt.numaper*numaper) then
                  kz=dsqrt(k0*k0-kx*kx-ky*ky) 
                  indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                 +nfft2d2+1)-1)

                  ctmp=(kx*Eloinx(indice)+ky*Eloiny(indice)+kz*dside
     $                 *Eloinz(indice))/k0
                  ctmp1=cdexp(-icomp*dside*kz*zs(kk))*cdexp(-icomp*(var1
     $                 *dble(i)+var2*dble(j)))/(-2.d0*pi*icomp*kz)

                  Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)
     $                 +(Eloinx(indice)-kx*ctmp/k0)*ctmp1
                  Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)
     $                 +(Eloiny(indice)-ky*ctmp/k0)*ctmp1
                  Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)
     $                 +(Eloinz(indice)-kz*dside*ctmp/k0)*ctmp1

                endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
      enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,ii,jj)
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(2)   
      do i=-imax,imax
         do j=-imax,imax
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
            ii=imax+i+1
            jj=imax+j+1
            Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)*fac*k0*k0
            Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)*fac*k0*k0
            Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)*fac*k0*k0
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      end
