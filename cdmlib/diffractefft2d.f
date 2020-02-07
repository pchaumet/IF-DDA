      subroutine diffractefft2d(nx,ny,nz,nxm,nym,nzm,nfft2d,k0,xs,ys,zs
     $     ,aretecube,Eloinx,Eloiny,Eloinz,FF,imax,deltakx,deltaky
     $     ,Ediffkzpos,Ediffkzneg,plan2f,plan2b,nstop,infostr)
      implicit none
      integer nx,ny,nz,nxm,nym,nzm,nfft2d,nstop,IFO
      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
     $     ,aretecube,k0
      double complex FF(3*nxm*nym*nzm),Ediffkzpos(nfft2d,nfft2d,3)
     $     ,Ediffkzneg(nfft2d,nfft2d,3)

      integer nfft2d2,imax,i,j,k,tabfft2(4096),indice,kk,ii,jj
      double precision deltakx,deltaky,var1,var2,kx,ky,kz,fac,pi
      double complex ctmp,ctmp1,icomp,Eloinx(nfft2d*nfft2d)
     $     ,Eloiny(nfft2d*nfft2d),Eloinz(nfft2d*nfft2d)
      character(64) infostr
      integer*8 plan2f,plan2b
      write(*,*) 'Quick method with FFT2',nfft2d
      
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltakx
      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltaky
c     fac=dble(nfft2d*nfft2d)
      fac=1.d0
      if (nfft2d.gt.4096) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      IFO=2
    
      if (deltakx.ge.k0) then
         nstop=1
         infostr='In FFT Poynting nfft2d too small'
         return
      endif
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
      write(*,*) 'deltakx',deltakx,imax,k0
      if (2*imax+1.gt.nfft2d) then
         write(99,*) '2*imax+1',imax,2*imax+1,nfft2d
         infostr='In FFT diffract nfft2d too small'
         nstop = 1
      return
      endif
      if (2*imax+1.lt.7) then
         write(99,*) '2*imax+1',imax,2*imax+1,nfft2d
         infostr='In FFT diffract nfft2d too small'
         nstop = 1
      return
      endif
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)     
      do i=1,nfft2d
         do j=1,nfft2d
            Ediffkzpos(i,j,1)=0.d0
            Ediffkzpos(i,j,2)=0.d0
            Ediffkzpos(i,j,3)=0.d0
            Ediffkzneg(i,j,1)=0.d0
            Ediffkzneg(i,j,2)=0.d0
            Ediffkzneg(i,j,3)=0.d0
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
#endif

         kk=1+nx*ny*(k-1)
         
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,ii,jj,kx,ky,kz,indice,ctmp,ctmp1)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
         do i=-imax,imax
            do j=-imax,imax
               ii=imax+i+1
               jj=imax+j+1
               kx=dble(i)*deltakx
               ky=dble(j)*deltaky
               if (k0*k0-kx*kx-ky*ky.gt.0.d0) then
                  kz=dsqrt(k0*k0-kx*kx-ky*ky) 
                  indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                 +nfft2d2+1)-1)

                  ctmp=(kx*Eloinx(indice)+ky*Eloiny(indice)+kz
     $                 *Eloinz(indice))/k0
                  ctmp1=cdexp(-icomp*kz*zs(kk))*cdexp(-icomp*(var1
     $                 *dble(i)+var2*dble(j)))

                  Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)
     $                 +(Eloinx(indice)-kx*ctmp/k0)*ctmp1
                  Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)
     $                 +(Eloiny(indice)-ky*ctmp/k0)*ctmp1
                  Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)
     $                 +(Eloinz(indice)-kz*ctmp/k0)*ctmp1

                  ctmp=(kx*Eloinx(indice)+ky*Eloiny(indice)-kz
     $                 *Eloinz(indice))/k0
                  ctmp1=cdexp(icomp*kz*zs(kk))*cdexp(-icomp*(var1
     $                 *dble(i)+var2*dble(j)))

                  Ediffkzneg(ii,jj,1) =Ediffkzneg(ii,jj,1)
     $                 +(Eloinx(indice)-kx*ctmp/k0)*ctmp1
                  Ediffkzneg(ii,jj,2) =Ediffkzneg(ii,jj,2)
     $                 +(Eloiny(indice)-ky*ctmp/k0)*ctmp1
                  Ediffkzneg(ii,jj,3) =Ediffkzneg(ii,jj,3)
     $                 +(Eloinz(indice)+kz*ctmp/k0)*ctmp1
c                  write(*,*) Ediffkzpos(ii,jj,1)*fac, Ediffkzpos(ii,jj
c     $                 ,2)*fac,Ediffkzpos(ii,jj,3)*fac,Ediffkzneg(ii,jj
c     $                 ,1)*fac,Ediffkzneg(ii ,jj,2)*fac,Ediffkzneg(ii,jj
c     $                 ,3)*fac,ii,jj
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      enddo


      write(*,*) 'end diffracte'
      end
