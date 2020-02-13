      subroutine diffractefft2dinc(nfft2d,k0,E0,ss,pp,theta,phi, thetam,
     $     phim, ppm, ssm,E0m,nbinc,xdip,ydip,zdip,xgaus,ygaus,zgaus,w0
     $     ,aretecube,tol,Eloinx,Eloiny,Eloinz ,imax,deltakx ,deltaky
     $     ,Ediffkzpos,numaper,beam,plan2f,plan2b,nstop,infostr)
      implicit none
      integer nfft2d,nstop
      double precision k0,numaper,aretecube
      double complex Ediffkzpos(nfft2d,nfft2d,3)


      integer nfft2d2,imax,i,j,k,indice,ii,jj
      double precision deltakx,deltaky,deltax,deltay,kx,ky,kz,fac,pi
      double complex icomp,Eloinx(nfft2d*nfft2d) ,Eloiny(nfft2d *nfft2d)
     $     ,Eloinz(nfft2d*nfft2d),Ex,Ey,Ez
      character(64) infostr

c     variable pour champ incident
      integer nloin,ikxinc,jkyinc
      double precision ss,pp,theta,phi,xdip ,ydip,zdip,xgaus,ygaus
     $     ,zgaus,w0,x,y,z,tol,tolc,Emod,Emodmax,kxinc,kyinc
      double complex E0,Em(3)
      character(64) beam
      integer nbinc
      double precision thetam(10), phim(10), ppm(10), ssm(10)
      double complex E0m(10)
      integer indicex,indicey
      integer*8 plan2f,plan2b
#ifndef USE_FFTW
      integer FFTW_FORWARD
      FFTW_FORWARD=-1
#endif
      write(*,*) 'Compute the FFT of the incident field',nfft2d
      
      nloin=0
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      x=0.d0
      y=0.d0
      z=0.d0


      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj)
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(2)      
      do i=-imax,imax
         do j=-imax,imax
            ii=imax+i+1
            jj=imax+j+1
            Ediffkzpos(ii,jj,1)=0.d0
            Ediffkzpos(ii,jj,2)=0.d0
            Ediffkzpos(ii,jj,3)=0.d0
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            
      if (nfft2d.gt.4096) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      
      fac=1.d0/deltakx/deltakx/dble(nfft2d*nfft2d)
      if (deltakx.ge.numaper) then
         nstop=1
         infostr='In FFT lens inc nfft2d too small'
         return
      endif
      nfft2d2=nfft2d/2    
      
      deltax=2.d0*pi/(dble(nfft2d)*deltakx)
      deltay=2.d0*pi/(dble(nfft2d)*deltaky)


      if (beam(1:11).eq.'pwavelinear') then
         kxinc=dsin(theta*pi/180.d0)*dcos(phi*pi/180.d0)*k0
         kyinc=dsin(theta*pi/180.d0)*dsin(phi*pi/180.d0)*k0
         call deltakroutine(kxinc,kyinc,deltakx,deltaky,k0,ikxinc
     $        ,jkyinc)
         kx=deltakx*dble(ikxinc)
         ky=deltakx*dble(jkyinc)
         call   ondeplanekxky(x,y,z,k0,E0,ss,pp,kx,ky,Ex,Ey,Ez,nstop
     $        ,infostr)
         ii=imax+ikxinc+1
         jj=imax+jkyinc+1
         Ediffkzpos(ii,jj,1)=Ex/deltakx/deltakx
         Ediffkzpos(ii,jj,2)=Ey/deltakx/deltakx
         Ediffkzpos(ii,jj,3)=Ez/deltakx/deltakx

         return
         
      elseif (beam(1:13).eq.'pwavecircular') then
         kxinc=dsin(theta*pi/180.d0)*dcos(phi*pi/180.d0)*k0
         kyinc=dsin(theta*pi/180.d0)*dsin(phi*pi/180.d0)*k0
         call deltakroutine(kxinc,kyinc,deltakx,deltaky,k0,ikxinc
     $        ,jkyinc)
         kx=deltakx*dble(ikxinc)
         ky=deltakx*dble(jkyinc)
         call ondecircekxky(x,y,z,k0,E0,ss,kx,ky,Ex,Ey,Ez)
         ii=imax+ikxinc+1
         jj=imax+jkyinc+1
         Ediffkzpos(ii,jj,1)=Ex/deltakx/deltakx
         Ediffkzpos(ii,jj,2)=Ey/deltakx/deltakx
         Ediffkzpos(ii,jj,3)=Ez/deltakx/deltakx

         return
         
      elseif (beam(1:15).eq.'wavelinearmulti') then
         
         do i=1,nbinc

            theta=thetam(i)
            phi=phim(i)
            kxinc=dsin(theta*pi/180.d0)*dcos(phi*pi/180.d0)*k0
            kyinc=dsin(theta*pi/180.d0)*dsin(phi*pi/180.d0)*k0
            call deltakroutine(kxinc,kyinc,deltakx,deltaky,k0,ikxinc
     $           ,jkyinc)
            kx=deltakx*dble(ikxinc)
            ky=deltakx*dble(jkyinc)
            ss=ssm(i)
            pp=ppm(i)
            E0=E0m(i)
            call   ondeplanekxky(x,y,z,k0,E0,ss,pp,kx,ky,Ex,Ey,Ez,nstop
     $           ,infostr)
            ii=imax+ikxinc+1
            jj=imax+jkyinc+1
            Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)+Ex/deltakx/deltakx
            Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)+Ey/deltakx/deltakx
            Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)+Ez/deltakx/deltakx
         enddo
         return
      endif

      
c     fac=dble(nfft2d*nfft2d)

c     calcul du champ incident à z=0 avec la maille deltax

c     attention si beam gaussien alors passe par le calcul du paraxial
c     pour accélerer le cacul
      if (beam(1:11).eq.'gwavelinear'.or.beam(1:14).eq
     $     .'gfftwavelinear') then
         Emodmax=0.d0
         do i=-nfft2d2,nfft2d2-1
            x=deltax*dble(i)
            do j=-nfft2d2,nfft2d2-1
               y=deltax*dble(j)
               call gaussianparalinear(x,y,z,xgaus,ygaus,zgaus,theta,
     $              phi,w0,k0,ss,pp,E0,Em(1),Em(2),Em(3),nstop,infostr)
               Emodmax=max(dreal(Em(1)*dconjg(Em(1))+Em(2)
     $              *dconjg(Em(2))+Em(3)*dconjg(Em(3))),Emodmax)
            enddo
         enddo
      endif
      if (beam(1:11).eq.'gwavecircular'.or.beam(1:16).eq
     $     .'gfftwavecircular') then
         Emodmax=0.d0
         do i=-nfft2d2,nfft2d2-1
            x=deltax*dble(i)
            do j=-nfft2d2,nfft2d2-1
               y=deltax*dble(j)
               call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus,theta,phi
     $              ,w0,k0,ss,E0,Em(1),Em(2),Em(3),nstop,infostr)
               Emodmax=max(dreal(Em(1)*dconjg(Em(1))+Em(2)
     $              *dconjg(Em(2))+Em(3)*dconjg(Em(3))),Emodmax)
            enddo
         enddo
      endif
      
      do i=-nfft2d2,nfft2d2-1
         x=deltax*dble(i)
         if (i.ge.0) then
            indicex=i+1
         else
            indicex=nfft2d+i+1
         endif
         do j=-nfft2d2,nfft2d2-1
            if (j.ge.0) then
               indicey=j+1
            else
               indicey=nfft2d+j+1
            endif
            
            y=deltax*dble(j)
            indice=indicex+nfft2d*(indicey-1)

            if (beam(1:11).eq.'pwavelinear') then              
               call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi,Em(1)
     $              ,Em(2),Em(3),nstop,infostr)            
            elseif (beam(1:13).eq.'pwavecircular') then
               call ondecirce(x,y,z,k0,E0,ss,theta,phi,Em(1),Em(2)
     $              ,Em(3))
               write(*,*) 'rrrr'
            elseif (beam(1:15).eq.'wavelinearmulti') then
               call ondeplanemulti(x,y,z,k0,E0m,ssm,ppm,thetam,phim
     $              ,nbinc,Em(1),Em(2),Em(3),nstop,infostr)               
            elseif (beam(1:7).eq.'antenna') then
c     ne rentre jamais ici
               call dipoleinc(xdip,ydip,zdip,theta,phi,x,y,z
     $              ,aretecube,k0,E0,Em(1),Em(2),Em(3),nstop
     $              ,infostr)                     
            elseif (beam(1:11).eq.'gwavelinear'.or.beam(1:14).eq
     $              .'gfftwavelinear') then
               call gaussianparalinear(x,y,z,xgaus,ygaus,zgaus,theta
     $              ,phi,w0,k0,ss,pp,E0,Em(1),Em(2),Em(3),nstop,infostr)
               Emod=dreal(Em(1)*dconjg(Em(1))+Em(2)*dconjg(Em(2))
     $              +Em(3)*dconjg(Em(3)))
               tolc=dsqrt(Emod/Emodmax)*10.d0
               if (tolc.ge.tol) then
                  call gaussianchamp(x,y,z,xgaus,ygaus,zgaus,theta ,phi
     $                 ,w0,k0,ss,pp,E0,Em(1),Em(2),Em(3),tolc,nloin
     $                 ,nstop,infostr)                  
               endif
            elseif (beam(1:13).eq.'gwavecircular'.or.beam(1:16).eq
     $              .'gfftwavecircular') then
               call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus,theta ,phi
     $              ,w0,k0,ss,E0,Em(1),Em(2),Em(3),nstop ,infostr)
               Emod=dreal(Em(1)*dconjg(Em(1))+Em(2)*dconjg(Em(2))
     $              +Em(3)*dconjg(Em(3)))
               tolc=dsqrt(Emod/Emodmax)*10.d0
               if (tolc.ge.tol) then
                  call gaussianchampcirc(x,y,z,xgaus,ygaus,zgaus,theta
     $                 ,phi,w0,k0,ss,E0,Em(1),Em(2),Em(3),tol,nloin)
               endif
            elseif (beam(1:15).eq.'gparawavelinear') then
               call gaussianparalinear(x,y,z,xgaus,ygaus,zgaus
     $              ,theta,phi,w0,k0,ss,pp,E0,Em(1),Em(2),Em(3)
     $              ,nstop,infostr)
            elseif (beam(1:17).eq.'gparawavecircular') then
               call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus,theta
     $              ,phi,w0,k0,ss,E0,Em(1),Em(2),Em(3),nstop
     $              ,infostr)
            endif
            Eloinx(indice)=Em(1)*fac
            Eloiny(indice)=Em(2)*fac
            Eloinz(indice)=Em(3)*fac
         enddo
      enddo

c     calcul de la FFT
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
         
      
c     sauvegarde dans le tableau que pour les kz interessant et ajouts
c     du -2 i pi gamma

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,kx,ky,indicex,indicey,indice,kz)
!$OMP DO  SCHEDULE(DYNAMIC) COLLAPSE(2) 
      do i=-imax,imax
         do j=-imax,imax
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky

            ii=imax+i+1
            jj=imax+j+1
            if (kx*kx+ky*ky.lt.numaper*numaper) then


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
               
               kz=dsqrt(k0*k0-kx*kx-ky*ky) 
               indice=indicex+nfft2d*(indicey-1)
               Ediffkzpos(ii,jj,1)=Eloinx(indice)
               Ediffkzpos(ii,jj,2)=Eloiny(indice)
               Ediffkzpos(ii,jj,3)=Eloinz(indice)
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      end
