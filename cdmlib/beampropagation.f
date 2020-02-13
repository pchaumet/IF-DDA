      subroutine beampropagation(xs,ys,zs,aretecube,k0,w0,E0,ss,pp,theta
     $     ,phi,xgaus,ygaus,zgaus,beam,epsilon,ndipole,nx,ny,nz,nxm,nym
     $     ,nzm ,nmax,nfft2d,imagex,imagey,imagez,FF0,FFloc,FF,plan2f
     $     ,plan2b,nstop ,infostr)
      
      implicit none

      integer i,j,k,l,kk,nfft2d,ii,jj,nstop,kkm,jjm,tabfft(nfft2d
     $     *nfft2d),comparaison,imaxk0,indicex,indicey
      double complex imagex(nfft2d*nfft2d),imagey(nfft2d*nfft2d)
     $     ,imagez(nfft2d*nfft2d),ctmp,E0,expik0a
      double precision xmax,xmin,ymax,ymin,x,y,z,xgaus,ygaus,zgaus
     $     ,lambda
      character(64) infostr
      character(64) beam
c     input data
      integer nx,ny,nz,nxm,nym,nzm,nmax,ndipole,indice
      double precision xs(nmax),ys(nmax),zs(nmax),k0,w0,ss,pp,theta,phi
     $     ,aretecube,kz,kx,ky,deltakx,deltaky,pi,kzc,c,quatpieps0,var1
     $     ,var2

      double complex epsilon(nmax,3,3),FFloc(3*nmax),FF(3*nmax),FF0(3
     $     *nmax),icomp,Ex,Ey,Ez,fac
      integer*8 plan2f,plan2b


      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      c=299792458.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      lambda=2.d0*pi/k0
      expik0a=cdexp(icomp*k0*aretecube)
      fac=1.d0/dble(nfft2d*nfft2d)
      
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
      var1=(dble(nfft2d/2)*aretecube)*deltakx-aretecube*deltakx/2.d0
      var2=(dble(nfft2d/2)*aretecube)*deltaky-aretecube*deltakx/2.d0

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

      write(*,*) 'int',k0/deltakx,theta*180.d0/pi
      l=1
      z=zs(1)
c     intialise le premier z au champ incident avec la phase
      do j=1,nfft2d
         do k=1,nfft2d       
            ii=k+(j-1)*nfft2d
            x=xmin+dble(k-kkm-1)*aretecube
            y=ymin+dble(j-jjm-1)*aretecube
            if (comparaison(x,y,z,xs(l),ys(l),z
     $           ,aretecube).eq.1.and.l.le.ndipole)  then
               tabfft(ii)=l
               kk=3*(l-1)
               FFloc(kk+1)=FF0(kk+1)*cdexp(icomp*k0*(cdsqrt(epsilon(l,1
     $              ,1))-1.d0)*aretecube)
               FFloc(kk+2)=FF0(kk+2)*cdexp(icomp*k0*(cdsqrt(epsilon(l,2
     $              ,2))-1.d0)*aretecube)
               FFloc(kk+3)=FF0(kk+3)*cdexp(icomp*k0*(cdsqrt(epsilon(l,3
     $              ,3))-1.d0)*aretecube)        
               Imagex(ii)=FFloc(kk+1)
               Imagey(ii)=FFloc(kk+2)
               Imagez(ii)=FFloc(kk+3)
               l=l+1
            else
               tabfft(ii)=0
               if (beam(1:11).eq.'pwavelinear') then
                  call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi ,Imagex(ii)
     $                 ,Imagey(ii),Imagez(ii),nstop,infostr)
               elseif (beam(1:13).eq.'pwavecircular') then
                  call ondecirce(x,y,z,k0,E0,ss,theta,phi,Imagex(ii)
     $                 ,Imagey(ii),Imagez(ii))
               elseif (beam(1:15).eq.'gparawavelinear') then
                  call gaussianparalinear(x,y,z,xgaus,ygaus,zgaus,theta
     $                 ,phi,w0,k0,ss,pp,E0,Imagex(ii),Imagey(ii)
     $                 ,Imagez(ii),nstop,infostr)
               elseif (beam(1:17).eq.'gparawavecircular') then
                  call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus ,theta
     $                 ,phi,w0,k0,ss,E0,Imagex(ii),Imagey(ii)
     $                 ,Imagez(ii),nstop,infostr)
               endif
            endif

         enddo
      enddo

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
                  ctmp=cdexp(-icomp*(k0-kz)*aretecube)
                  Imagex(indice)=Imagex(indice)*ctmp
                  Imagey(indice)=Imagey(indice)*ctmp
                  Imagez(indice)=Imagez(indice)*ctmp
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indice,x,y,ii,kk)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)  
         do j=1,nfft2d
            do i=1,nfft2d         
               indice=i+(j-1)*nfft2d
               x=xmin+dble(i-kkm-1)*aretecube
               y=ymin+dble(j-jjm-1)*aretecube
               if (tabfft(indice).ne.0)  then
                  
                  ii=nx*ny*(k-1)+tabfft(indice)
                  kk=3*(ii-1)
                  
                  
                  FFloc(kk+1)=Imagex(indice)*cdexp(icomp*k0
     $                 *cdsqrt(epsilon(ii,1,1))*aretecube)*fac
                  FFloc(kk+2)=Imagey(indice)*cdexp(icomp*k0
     $                 *cdsqrt(epsilon(ii,2,2))*aretecube)*fac
                  FFloc(kk+3)=Imagez(indice)*cdexp(icomp*k0
     $                 *cdsqrt(epsilon(ii,3,3))*aretecube)*fac

                  Imagex(indice)=FFloc(kk+1)
                  Imagey(indice)=FFloc(kk+2)
                  Imagez(indice)=FFloc(kk+3)                 
               else
                 Imagex(indice)=Imagex(indice)*expik0a*fac
                 Imagey(indice)=Imagey(indice)*expik0a*fac
                 Imagez(indice)=Imagez(indice)*expik0a*fac
               endif
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
c     stop
      enddo



c     calcul du champ diffracte lointain a partir de la derniere couche
c     sans l'incident en faisant du kirchhoff
c$$$
c$$$      do j=1,nfft2d
c$$$         do i=1,nfft2d         
c$$$            indice=i+(j-1)*nfft2d
c$$$            x=xmin+dble(i-kkm-1)*aretecube
c$$$            y=ymin+dble(j-jjm-1)*aretecube
c$$$            z=zs(ndipole)
c$$$            if (beam(1:11).eq.'pwavelinear') then
c$$$               call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi ,Ex,Ey,Ez
c$$$     $              ,nstop,infostr)
c$$$            elseif (beam(1:13).eq.'pwavecircular') then
c$$$               call ondecirce(x,y,z,k0,E0,ss,theta,phi,Ex,Ey,Ez)
c$$$            elseif (beam(1:15).eq.'gparawavelinear') then
c$$$               call gaussianparalinear(x,y,z,xgaus,ygaus,zgaus,theta
c$$$     $              ,phi,w0,k0,ss,pp,E0,Ex,Ey,Ez)
c$$$            elseif (beam(1:17).eq.'gparawavecircular') then
c$$$               call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus ,theta
c$$$     $              ,phi,w0,k0,ss,E0,Ex,Ey,Ez,nstop,infostr)
c$$$            endif
c$$$            Imagex(indice)=Imagex(indice)-Ex
c$$$            Imagey(indice)=Imagey(indice)-Ey
c$$$            Imagez(indice)=Imagez(indice)-Ez
c$$$         enddo
c$$$      enddo
c$$$
c$$$
c$$$!     $OMP PARALLEL 
c$$$!     $OMP SECTIONS 
c$$$!     $OMP SECTION   
c$$$      CALL ZFFT2D(Imagex,nfft2d,nfft2d,-1)
c$$$!     $OMP SECTION 
c$$$      CALL ZFFT2D(Imagey,nfft2d,nfft2d,-1)
c$$$!     $OMP SECTION 
c$$$      CALL ZFFT2D(Imagez,nfft2d,nfft2d,-1)
c$$$!     $OMP END SECTIONS
c$$$!     $OMP END PARALLEL      
c$$$
c$$$      fac=dble(nfft2d*nfft2d)*aretecube*aretecube/icomp/4.d0/pi
c$$$     $     *cdexp(-icomp*k0*zs(ndipole))
c$$$      imaxk0=nint(k0/deltakx)+1
c$$$      do i=-imaxk0,imaxk0               
c$$$         if (i.ge.0) then
c$$$            indicex=i+1
c$$$         else
c$$$            indicex=nfft2d+i+1
c$$$         endif
c$$$
c$$$         do j=-imaxk0,imaxk0
c$$$            if (j.ge.0) then
c$$$               indicey=j+1
c$$$            else
c$$$               indicey=nfft2d+j+1
c$$$            endif
c$$$            indice=indicex+nfft2d*(indicey-1)
c$$$            kx=deltakx*dble(i)
c$$$            ky=deltaky*dble(j)
c$$$            kz=k0*k0-kx*kx-ky*ky
c$$$            ctmp=cdexp(icomp*(var1*dble(i)+var2*dble(j)))
c$$$            if (kz.ge.0.d0) then
c$$$               kz=dsqrt(kz)
c$$$               Ex=Imagex(indice)*fac*ctmp
c$$$               Ey=Imagey(indice)*fac*ctmp
c$$$               Ez=Imagez(indice)*fac*ctmp
c$$$               ctmp=kx*Ex+ky*Ey+kz*Ez
c$$$               Imagex(indice)=2.d0*kz*Ex-ctmp*kx*kz/k0/k0
c$$$               Imagey(indice)=2.d0*kz*Ey-ctmp*ky*kz/k0/k0
c$$$               Imagez(indice)=2.d0*kz*Ez-ctmp*(1.d0+kz*kz/k0/k0)
c$$$
c$$$            else
c$$$               Imagex(indice)=0.d0
c$$$               Imagey(indice)=0.d0
c$$$               Imagez(indice)=0.d0   
c$$$            endif
c$$$c            write(1001,*) (cdabs(Imagex(indice))**2.d0
c$$$c     $           +cdabs(Imagey(indice))**2.d0 +cdabs(Imagez(indice))
c$$$c     $           **2.d0)*c/8.d0/pi*quatpieps0
c$$$               
c$$$         enddo
c$$$      enddo    

      end
