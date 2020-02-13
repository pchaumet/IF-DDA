      subroutine gaussianchampfft(xs,ys,zs,aretecube,k0,w0,E0,ss,pp
     $     ,theta,phi,xgaus,ygaus,zgaus,beam,ndipole,nx,ny,nz,nxm,nym
     $     ,nzm ,nmax,nfft2d,imagex,imagey,imagez,FF0,tolinit,plan2f
     $     ,plan2b,nstop ,infostr)
      
      implicit none

      integer i,j,k,l,kk,nfft2d,ii,jj,nstop,kkm,jjm,tabfft(nfft2d
     $     *nfft2d),comparaison,imaxk0,indicex,indicey,nloin,nxy
      double complex imagex(nfft2d*nfft2d),imagey(nfft2d*nfft2d)
     $     ,imagez(nfft2d*nfft2d),ctmp,E0
      double precision xmax,xmin,ymax,ymin,x,y,z,xgaus,ygaus,zgaus
     $     ,lambda,tol,tolinit
      character(64) infostr
c     input data
      integer nx,ny,nz,nxm,nym,nzm,nmax,ndipole,indice
      double precision xs(nmax),ys(nmax),zs(nmax),k0,w0,ss,pp,theta,phi
     $     ,aretecube,kz,kx,ky,deltakx,deltaky,pi,kzc,c,quatpieps0

      double complex FF0(3*nmax),icomp,Ex,Ey,Ez,fac
c     output data
      character(64) beam
      integer*8 plan2f ,plan2b
      
#ifndef USE_FFTW
      integer FFTW_BACKWARD, FFTW_FORWARD
      FFTW_BACKWARD=+1
      FFTW_FORWARD=-1
#endif      
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      lambda=2.d0*pi/k0
      nxy=nx*ny
      fac=1.d0/dble(nfft2d*nfft2d)
      
c     arret si theta trop fort car ne prend que les kz positifs
      if (dabs(theta).ge.45.d0) then
         infostr='Angle too large for fft gaussian'
         nstop=1
      endif

      
c     verifie taille fft
      if (nfft2d.le.nx) then
         nstop=1
         infostr='nfft2d too small'
      endif
      if (nfft2d.le.ny) then
         nstop=1
         infostr='nfft2d too small'
      endif
      
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=2.d0*pi/(dble(nfft2d)*aretecube)

      if (deltakx.ge.k0/2.d0) then
         nstop=1
         infostr='nfft2d too small'
      endif
      if (deltaky.ge.k0/2.d0) then
         nstop=1
         infostr='nfft2d too small'
      endif
      if (nstop.eq.1) return

      xmax=-1.d300
      xmin=1.d300
      ymax=-1.d300
      ymin=1.d300

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:xmax,ymax)
!$OMP& REDUCTION(min:xmin,ymin)
      do i=1,ndipole
         xmax=max(xmax,xs(i))
         xmin=min(xmin,xs(i))
         ymax=max(ymax,ys(i))
         ymin=min(ymin,ys(i))     
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      kkm=(nfft2d-nx)/2
      jjm=(nfft2d-ny)/2

      l=1
      z=zs(1)
      nloin=0
      tol=tolinit
        
c     intialise le premier z au champ incident avec la phase
      do j=1,nfft2d
         do k=1,nfft2d       
            ii=k+(j-1)*nfft2d
            x=xmin+dble(k-kkm-1)*aretecube
            y=ymin+dble(j-jjm-1)*aretecube
 
            if (beam(1:14).eq.'gfftwavelinear') then
               call gaussianchamp(x,y,z,xgaus,ygaus,zgaus,theta ,phi,w0
     $              ,k0,ss,pp,E0,Imagex(ii),Imagey(ii),Imagez(ii),tol
     $              ,nloin,nstop ,infostr)
            elseif (beam(1:16).eq.'gfftwavecircular') then
               call gaussianchampcirc(x,y,z,xgaus,ygaus,zgaus ,theta
     $              ,phi,w0,k0,ss,E0,Imagex(ii),Imagey(ii),Imagez(ii)
     $              ,tol,nloin,nstop,infostr)
            endif
            
            if (comparaison(x,y,z,xs(l),ys(l),z,aretecube) .eq.1 .and.
     $           l.le.ndipole)  then
               tabfft(ii)=l
               l=l+1
            else
               tabfft(ii)=0
            endif
            
         enddo
      enddo

c     calcul premiere couche en filtrant les frequences hors propagatives
#ifdef USE_FFTW
      call dfftw_execute_dft(plan2b,Imagex,Imagex)
      call dfftw_execute_dft(plan2b,Imagey,Imagey)
      call dfftw_execute_dft(plan2b,Imagez,Imagez)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      call fftsingletonz2d(Imagex,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Imagey,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
      call fftsingletonz2d(Imagez,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,ii,jj,kx,ky,kz)
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(2)          
         do j=1,nfft2d
            do i=1,nfft2d
               indice=i+nfft2d*(j-1)
               if (i.gt.nfft2d/2) then
                  ii=i-nfft2d-1
               else
                  ii=i-1
               endif
               if (j.gt.nfft2d/2) then
                  jj=j-nfft2d-1
               else
                  jj=j-1
               endif              
               kx=deltakx*dble(ii)
               ky=deltaky*dble(jj)
               kz=k0*k0-kx*kx-ky*ky
               if (kz.le.0.d0) then         
                  Imagex(indice)=0.d0
                  Imagey(indice)=0.d0
                  Imagez(indice)=0.d0
               else
                  Imagex(indice)=Imagex(indice)*fac
                  Imagey(indice)=Imagey(indice)*fac
                  Imagez(indice)=Imagey(indice)*fac
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

#ifdef USE_FFTW
         call dfftw_execute_dft(plan2f,Imagex,Imagex)
         call dfftw_execute_dft(plan2f,Imagey,Imagey)
         call dfftw_execute_dft(plan2f,Imagez,Imagez)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Imagex,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagey,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagez,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif
         
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i,j,indice,x,y,ii,kk)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2) 
         do j=1,nfft2d
            do i=1,nfft2d         
               indice=i+(j-1)*nfft2d
               x=xmin+dble(i-kkm-1)*aretecube
               y=ymin+dble(j-jjm-1)*aretecube
               
               if (tabfft(indice).ne.0)  then
                  
                  ii=tabfft(indice)
                  kk=3*(ii-1)                
                  FF0(kk+1)=Imagex(indice)
                  FF0(kk+2)=Imagey(indice)
                  FF0(kk+3)=Imagez(indice)
                  
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         
c     commence boucle sur les z
      do k=2,nz

c     FFT à deux dimensions
#ifdef USE_FFTW
         call dfftw_execute_dft(plan2b,Imagex,Imagex)
         call dfftw_execute_dft(plan2b,Imagey,Imagey)
         call dfftw_execute_dft(plan2b,Imagez,Imagez)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Imagex,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagey,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagez,nfft2d,nfft2d,FFTW_BACKWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,ii,jj,kx,ky,kz,ctmp)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)      
         do j=1,nfft2d
            do i=1,nfft2d
               indice=i+nfft2d*(j-1)
               if (i.gt.nfft2d/2) then
                  ii=i-nfft2d-1
               else
                  ii=i-1
               endif
               if (j.gt.nfft2d/2) then
                  jj=j-nfft2d-1
               else
                  jj=j-1
               endif              
               kx=deltakx*dble(ii)
               ky=deltaky*dble(jj)
               kz=k0*k0-kx*kx-ky*ky
               if (kz.ge.0.d0) then
                  kz=dsqrt(kz)
                  ctmp=cdexp(icomp*kz*aretecube)
                  Imagex(indice)=Imagex(indice)*ctmp*fac
                  Imagey(indice)=Imagey(indice)*ctmp*fac
                  Imagez(indice)=Imagez(indice)*ctmp*fac
               else
                  Imagex(indice)=0.d0
                  Imagey(indice)=0.d0
                  Imagez(indice)=0.d0
                  
               endif

            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

#ifdef USE_FFTW
         call dfftw_execute_dft(plan2f,Imagex,Imagex)
         call dfftw_execute_dft(plan2f,Imagey,Imagey)
         call dfftw_execute_dft(plan2f,Imagez,Imagez)
#else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
         call fftsingletonz2d(Imagex,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagey,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP SECTION   
         call fftsingletonz2d(Imagez,nfft2d,nfft2d,FFTW_FORWARD)
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,x,y,ii,kk)
!$OMP DOSCHEDULE(STATIC) COLLAPSE(2) 
         do j=1,nfft2d
            do i=1,nfft2d         
               indice=i+(j-1)*nfft2d
               x=xmin+dble(i-kkm-1)*aretecube
               y=ymin+dble(j-jjm-1)*aretecube
               
               if (tabfft(indice).ne.0)  then
                  
                  ii=nxy*(k-1)+tabfft(indice)
                  kk=3*(ii-1)
                
                  FF0(kk+1)=Imagex(indice)
                  FF0(kk+2)=Imagey(indice)
                  FF0(kk+3)=Imagez(indice)
                  
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
c     stop
      enddo

      write(*,*) 'fin fft gaussian'
      end
