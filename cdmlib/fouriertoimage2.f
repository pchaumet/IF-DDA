      subroutine fouriertoimage2(deltakx,deltaky,gross,Eimagex,Eimagey
     $     ,Eimagez,nfft2D,nfft2d2 ,plan2b,plan2f)
      
      implicit none
      integer nfft2D,nfft2d2
      double precision deltakx,deltaky,gross
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d)

      integer i,j,indicex,indicey,indice,kk
      double precision tmp
      double complex ctmp
      integer*8 plan2f,plan2b

      tmp=deltakx*deltaky/gross
c     FFT comme dans la diffraction a cause du grossissement negatif qui
c     repasse la FFT de inverse à directe.
      
#ifdef USE_FFTW
      call dfftw_execute_dft(plan2f,Eimagex,Eimagex)
      call dfftw_execute_dft(plan2f,Eimagey,Eimagey)
      call dfftw_execute_dft(plan2f,Eimagez,Eimagez)
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

         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c$$$
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
c$$$         enddo
c$$$      enddo
c$$$!$OMP ENDDO 
c$$$!$OMP END PARALLEL




      end
