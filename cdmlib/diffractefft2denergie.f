c     Cette fonction calcul le champ diffracté sur un plan, d'ou la
c     correction en -2i pi gamma.
      subroutine diffractefft2denergie(nx,ny,nz,nxm,nym,nzm,nfft2d,k0,xs
     $     ,ys,zs,E0,ss,pp,theta,phi,thetam ,phim, ppm, ssm,E0m,nbinc
     $     ,xdip ,ydip,zdip,xgaus,ygaus,zgaus ,w0,aretecube,tol,Eloinx
     $     ,Eloiny,Eloinz,FF,imax,deltakx ,deltaky,Ediffkzpos,Ediffkzneg
     $     ,beam,efficacite,efficaciteref ,efficacitetrans,nsectionsca
     $     ,nquickdiffracte,plan2f,plan2b,nstop ,infostr)
      implicit none
      integer nx,ny,nz,nxm,nym,nzm,nfft2d,nstop,nsectionsca
     $     ,nquickdiffracte
      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
     $     ,aretecube,k0,k02
      double complex FF(3*nxm*nym*nzm),Ediffkzpos(nfft2d,nfft2d,3)
     $     ,Ediffkzneg(nfft2d,nfft2d,3)

      integer nfft2d2,imax,i,j,k,tabfft2(4096),indice,kk,ii,jj
      double precision deltakx,deltaky,var1,var2,kx,ky,kz,fac,pi
      double complex ctmp,ctmp1,icomp,Eloinx(nfft2d*nfft2d)
     $     ,Eloiny(nfft2d*nfft2d),Eloinz(nfft2d*nfft2d)
      character(64) infostr

c     variable pour champ incident
      integer nloin,indicex,indicey
      double precision ss,pp,theta,phi,xdip ,ydip,zdip,xgaus,ygaus
     $     ,zgaus,w0,x,y,z,tol,deltax,deltay,tolc,Emod,Emodmax
      double complex E0,Em(3)
      character(64) beam
      integer nbinc
      double precision thetam(10), phim(10), ppm(10), ssm(10)
      double complex E0m(10)
c     energie
      double precision fluxinc,fluxref,fluxtra,efficacite,efficaciteref
     $     ,efficacitetrans
      double precision tmp1,tmp2
      integer*8 plan2f,plan2b
      
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      deltakx=2.d0*pi/(dble(nfft2d)*aretecube)
      deltaky=deltakx
      nfft2d2=nfft2d/2  
      var1=(xs(1)+dble(nfft2d2)*aretecube)*deltakx
      var2=(ys(1)+dble(nfft2d2)*aretecube)*deltaky
      k02=k0*k0
      write(*,*) 'delta k',deltakx
      if (nfft2d.gt.4096) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      
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
      write(*,*) 'Number of point in the numerical aperture',imax
      write(*,*) 'Size of FFT',nfft2d
      
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


      if (nsectionsca*nquickdiffracte.eq.0) then
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j) 
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)          
         do i=1,nfft2d
            do j=1,nfft2d  
               Ediffkzneg(i,j,1)=0.d0
               Ediffkzneg(i,j,2)=0.d0
               Ediffkzneg(i,j,3)=0.d0
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
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2) 
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
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)             
            do i=-imax,imax
               do j=-imax,imax
                  ii=imax+i+1
                  jj=imax+j+1
                  kx=dble(i)*deltakx
                  ky=dble(j)*deltaky
                  if (k0*k0-kx*kx-ky*ky.gt.0.d0) then

                     kz=dsqrt(k0*k0-kx*kx-ky*ky) 
                     
                     indice=tabfft2(i+nfft2d2+1)+nfft2d*(tabfft2(j
     $                    +nfft2d2+1)-1)

                     ctmp=(kx*Eloinx(indice)+ky*Eloiny(indice)+kz
     $                    *Eloinz(indice))/k0
                     ctmp1=cdexp(-icomp*kz*zs(kk))*cdexp(-icomp*(var1
     $                    *dble(i)+var2*dble(j)))

                     Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)
     $                    +(Eloinx(indice)-kx*ctmp/k0)*ctmp1
                     Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)
     $                    +(Eloiny(indice)-ky*ctmp/k0)*ctmp1
                     Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)
     $                    +(Eloinz(indice)-kz*ctmp/k0)*ctmp1

                     ctmp=(kx*Eloinx(indice)+ky*Eloiny(indice)-kz
     $                    *Eloinz(indice))/k0
                     ctmp1=cdexp(icomp*kz*zs(kk))*cdexp(-icomp*(var1
     $                    *dble(i)+var2*dble(j)))

                     Ediffkzneg(ii,jj,1) =Ediffkzneg(ii,jj,1)
     $                    +(Eloinx(indice)-kx*ctmp/k0)*ctmp1
                     Ediffkzneg(ii,jj,2) =Ediffkzneg(ii,jj,2)
     $                    +(Eloiny(indice)-ky*ctmp/k0)*ctmp1
                     Ediffkzneg(ii,jj,3) =Ediffkzneg(ii,jj,3)
     $                    +(Eloinz(indice)+kz*ctmp/k0)*ctmp1
                     
                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         enddo 
      endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,kz,ctmp,ii,jj)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
      do i=-imax,imax
         do j=-imax,imax
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
            if (k0*k0-kx*kx-ky*ky.gt.0.d0) then               
               kz=dsqrt(k0*k0-kx*kx-ky*ky)
               ctmp=-2.d0*pi*icomp*kz
               ii=imax+i+1
               jj=imax+j+1
               Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)*k02/ctmp
               Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)*k02/ctmp
               Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)*k02/ctmp
               Ediffkzneg(ii,jj,1)=Ediffkzneg(ii,jj,1)*k02/ctmp
               Ediffkzneg(ii,jj,2)=Ediffkzneg(ii,jj,2)*k02/ctmp
               Ediffkzneg(ii,jj,3)=Ediffkzneg(ii,jj,3)*k02/ctmp      
            
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL


c     calcul pour rajouter le champ incident
      deltax=2.d0*pi/(dble(nfft2d)*deltakx)
      deltay=2.d0*pi/(dble(nfft2d)*deltaky)
      write(*,*) 'Size of the pixel',deltax
      z=0.d0
      nloin=0
      fac=1.d0/deltakx/deltakx/dble(nfft2d*nfft2d)

c     attention si beam gaussien alors passe par le calcul du paraxial
c     pour accélerer le calcul
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
            elseif (beam(1:15).eq.'wavelinearmulti') then
               call ondeplanemulti(x,y,z,k0,E0m,ssm,ppm,thetam,phim
     $              ,nbinc,Em(1),Em(2),Em(3),nstop,infostr) 
            elseif (beam(1:7).eq.'antenna') then
c     arret avant ne rentre jamais ici
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
#endif
         
      fluxinc=0.d0
      fluxref=0.d0
      fluxtra=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,kx,ky,ii,jj,indicex,indicey,indice,kz)
!$OMP DO SCHEDULE(DYNAMIC)  COLLAPSE(2)
!$OMP& REDUCTION(+:fluxinc,fluxref,fluxtra)         
      do i=-imax,imax
         do j=-imax,imax
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky

            ii=imax+i+1
            jj=imax+j+1
            if (kx*kx+ky*ky.lt.k0*k0) then


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
c               write(*,*) 'ij',i,j,imax,(dreal(Eloinx(indice)
c     $              *dconjg(Eloinx(indice)))+dreal(Eloiny(indice)
c     $              *dconjg(Eloiny(indice)))+dreal(Eloinz(indice)
c     $              *dconjg(Eloinz(indice)))),cdabs(E0)/deltakx/deltaky
c     $              ,kz
               
               fluxinc=fluxinc+(dreal(Eloinx(indice)
     $              *dconjg(Eloinx(indice)))+dreal(Eloiny(indice)
     $              *dconjg(Eloiny(indice)))+dreal(Eloinz(indice)
     $              *dconjg(Eloinz(indice))))*kz
               
c               Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)+Eloinx(indice)
c               Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)+Eloiny(indice)
c               Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)+Eloinz(indice)
c     write(*,*) 'Ediff2',Ediffkzpos(ii,jj,1) ,Ediffkzpos(ii,jj
c     $              ,2),Ediffkzpos(ii,jj,3)
               fluxref=fluxref+(dreal(Ediffkzneg(ii,jj,1)
     $              *dconjg(Ediffkzneg(ii,jj,1))) + dreal(Ediffkzneg(ii
     $              ,jj,2) *dconjg(Ediffkzneg(ii,jj,2))) +
     $              dreal(Ediffkzneg(ii,jj,3)*dconjg(Ediffkzneg(ii,jj
     $              ,3))))*kz
               fluxtra=fluxtra+(dreal((Ediffkzpos(ii,jj,1)
     $              +Eloinx(indice))*dconjg(Ediffkzpos(ii,jj,1)
     $              +Eloinx(indice))) + dreal((Ediffkzpos(ii,jj,2)
     $              +Eloiny(indice))*dconjg(Ediffkzpos(ii,jj,2)
     $              +Eloiny(indice)))+dreal((Ediffkzpos(ii,jj,3)
     $              +Eloinz(indice))*dconjg(Ediffkzpos(ii,jj,3)
     $              +Eloinz(indice))))*kz
               
            endif
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     

      write(*,*) 'Flux reflexion   : ',fluxref
      write(*,*) 'Flux trasmission : ',fluxtra
      write(*,*) 'Flux incident    : ',fluxinc
      efficacite=(fluxref+fluxtra)/fluxinc
      efficaciteref=fluxref/fluxinc
      efficacitetrans=fluxtra/fluxinc
      
      tmp1=0.d0
      tmp2=0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,indicex,indicey,indice)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2) REDUCTION(+:tmp1,tmp2)
      do i=-nfft2d2,nfft2d2-1
         do j=-nfft2d2,nfft2d2-1
            kx=dble(i)*deltakx
            ky=dble(j)*deltaky
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
            
            if (kx*kx+ky*ky.lt.k0*k0) then
               tmp1=tmp1+dreal(Eloinx(indice) *dconjg(Eloinx(indice)))
     $              +dreal(Eloiny(indice) *dconjg(Eloiny(indice)))
     $              +dreal(Eloinz(indice) *dconjg(Eloinz(indice)))
            endif
            tmp2=tmp2+dreal(Eloinx(indice) *dconjg(Eloinx(indice)))
     $           +dreal(Eloiny(indice) *dconjg(Eloiny(indice)))
     $           +dreal(Eloinz(indice) *dconjg(Eloinz(indice)))
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      write(*,*) tmp1,tmp2
      write(*,*) 'Energy outside the NA (%)',dabs(tmp1-tmp2)/tmp2*100.d0
      if (dabs(tmp1-tmp2)/tmp2*100.d0.ge.1.d0) then
      infostr='Energy: Beam too inclined or waist or nfft too small'
         nstop=1
      endif
      
      end
