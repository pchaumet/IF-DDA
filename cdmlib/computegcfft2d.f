      subroutine computegcfft2d(imax,deltakx,deltaky,k0,I0,theta,phi
     $     ,nfft2d,Ediffkzpos,Ediffkzneg,gasym,Cscai)
      implicit none
      integer i,j,imax,ii,jj,nfft2d
      double precision Cscai,gasym,kx,ky,kz,deltakx,deltaky,k0,Emod,I0
     $     ,k03,theta,phi,pi
      double complex Ediffkzpos(nfft2d,nfft2d,3),Ediffkzneg(nfft2d
     $     ,nfft2d,3)
      
c      write(*,*) 'coucou11',imax,deltakx,deltaky,k0,I0,theta,phi ,nfft2d
      Cscai=0.d0
      gasym=0.d0
      k03=k0*k0*k0
      pi=dacos(-1.d0)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,kz,ii,jj,Emod)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2) REDUCTION(+:Cscai,gasym)
      do i=-imax,imax
         do j=-imax,imax

            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
            if (k0*k0-kx*kx-ky*ky.gt.0.d0) then 
                           
               kz=dsqrt(k0*k0-kx*kx-ky*ky)
               ii=imax+i+1
               jj=imax+j+1

               
               Emod=cdabs(Ediffkzpos(ii,jj,1))**2+cdabs(Ediffkzpos(ii,jj
     $              ,2))**2+cdabs(Ediffkzpos(ii,jj,3))**2
               Cscai=Cscai+deltakx*deltaky*Emod/kz
               gasym=gasym+deltakx*deltaky*Emod/kz*(kx *dsin(theta
     $              *pi/180.d0)*dcos(phi*pi /180.d0) +ky
     $              *dsin(theta*pi/180.d0)* dsin(phi *pi /180.d0)
     $              +kz*dcos(theta*pi/180.d0))/k0
               
               Emod=cdabs(Ediffkzneg(ii,jj,1))**2+cdabs(Ediffkzneg(ii,jj
     $              ,2))**2+cdabs(Ediffkzneg(ii,jj,3))**2
               
               Cscai=Cscai+deltakx*deltaky*Emod/kz
               gasym=gasym+deltakx*deltaky*Emod/kz*(kx *dsin(theta
     $              *pi/180.d0)*dcos(phi*pi /180.d0) +ky
     $              *dsin(theta*pi/180.d0)*dsin(phi*pi /180.d0)-kz
     $              *dcos(theta*pi/180.d0))/k0

            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      gasym=gasym/cscai
      Cscai=Cscai*k03/I0

      end
