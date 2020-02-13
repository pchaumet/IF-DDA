      subroutine fouriertoimage(deltakx,deltaky,gross,Eimagex,Eimagey
     $     ,Eimagez,Eimageincx,Eimageincy,Eimageincz,nfft2D,nfft2d2
     $     ,plan2b,plan2f)
      
      implicit none
      integer nfft2D,nfft2d2
      double precision deltakx,deltaky,gross
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Eimageincx(nfft2d*nfft2d)
     $     ,Eimageincy(nfft2d *nfft2d),Eimageincz(nfft2d*nfft2d)

      integer i,j,indicex,indicey,indice,kk
      double precision tmp
      double complex ctmp
      integer*8 plan2f,plan2b

      tmp=deltakx*deltaky/gross
c     FFT comme dans la diffraction a cause du grossissement negatif qui
c     repasse la FFT de inverse à directe.
      
#ifdef USE_FFTW
      call dfftw_execute_dft(plan2f,Eimageincx,Eimageincx)
      call dfftw_execute_dft(plan2f,Eimageincy,Eimageincy)
      call dfftw_execute_dft(plan2f,Eimageincz,Eimageincz)
      call dfftw_execute_dft(plan2f,Eimagex,Eimagex)
      call dfftw_execute_dft(plan2f,Eimagey,Eimagey)
      call dfftw_execute_dft(plan2f,Eimagez,Eimagez)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz2d(Eimageincx,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimageincy,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimageincz,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagex,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagey,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
      call fftsingletonz2d(Eimagez,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,indicex,indicey,indice,kk,ctmp)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)           
      do i=-nfft2d2,nfft2d2-1
         do j=-nfft2d2,-1
            if (i.ge.0) then
               indicex=i+1
            else
               indicex=nfft2d+i+1
            endif
            if (j.ge.0) then
               indicey=j+1
            else
               indicey=nfft2d+j+1
            endif

            indice=indicex+nfft2d*(indicey-1)
            kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
            
            ctmp=Eimagex(kk)
            Eimagex(kk)=Eimagex(indice)*tmp
            Eimagex(indice)=ctmp*tmp
            ctmp=Eimagey(kk)
            Eimagey(kk)=Eimagey(indice)*tmp
            Eimagey(indice)=ctmp*tmp
            ctmp=Eimagez(kk)
            Eimagez(kk)=Eimagez(indice)*tmp
            Eimagez(indice)=ctmp*tmp

            ctmp=Eimageincx(kk)
            Eimageincx(kk)=Eimageincx(indice)*tmp
            Eimageincx(indice)=ctmp*tmp
            ctmp=Eimageincy(kk)
            Eimageincy(kk)=Eimageincy(indice)*tmp
            Eimageincy(indice)=ctmp*tmp
            ctmp=Eimageincz(kk)
            Eimageincz(kk)=Eimageincz(indice)*tmp
            Eimageincz(indice)=ctmp*tmp
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      
c$$$!$OMP PARALLEL DEFAULT(SHARED)
c$$$!$OMP& PRIVATE(i,j,indice,kk,ctmp)
c$$$!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)           
c$$$      do i=-nfft2d2,nfft2d2-1
c$$$         do j=-nfft2d2,-1
c$$$            kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
c$$$            indice=nfft2d*nfft2d+1-kk
c$$$
c$$$            ctmp=Eimagex(kk)
c$$$            Eimagex(kk)=Eimagex(indice)
c$$$            Eimagex(indice)=ctmp
c$$$            ctmp=Eimagey(kk)
c$$$            Eimagey(kk)=Eimagey(indice)
c$$$            Eimagey(indice)=ctmp
c$$$            ctmp=Eimagez(kk)
c$$$            Eimagez(kk)=Eimagez(indice)
c$$$            Eimagez(indice)=ctmp
c$$$
c$$$            ctmp=Eimageincx(kk)
c$$$            Eimageincx(kk)=Eimageincx(indice)
c$$$            Eimageincx(indice)=ctmp
c$$$            ctmp=Eimageincy(kk)
c$$$            Eimageincy(kk)=Eimageincy(indice)
c$$$            Eimageincy(indice)=ctmp
c$$$            ctmp=Eimageincz(kk)
c$$$            Eimageincz(kk)=Eimageincz(indice)
c$$$            Eimageincz(indice)=ctmp
c$$$         enddo
c$$$      enddo
c$$$!$OMP ENDDO 
c$$$!$OMP END PARALLEL


      end
