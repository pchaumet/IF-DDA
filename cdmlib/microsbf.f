      subroutine microsbf(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $     ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $     ,ntotalm,ntotal,ldabi,nlar,nmax,nxm,nym,nzm,nx,ny,nz ,nx2,ny2
     $     ,nxy2,nz2,ndipole,nbsphere,nbsphere3,nproche,nrig,nfft2d
     $     ,tabdip,XI ,XR,wrk ,FF,FF0 ,FFloc ,polarisa ,epsilon
     $     ,methodeit,tolinit,tol1 ,nloop ,ncompte,xs ,ys,zs ,aretecube
     $     ,numaper,numaperinc,npolainc,nquicklens,eps0,k0,P0 ,irra,w0
     $     ,gross,zlens,Eimagex ,Eimagey,Eimagez ,Efourierx,Efouriery
     $     ,Efourierz,Eimageincx ,Eimageincy ,Eimageincz ,Efourierincx
     $     ,Efourierincy ,Efourierincz ,Ediffkzpos,Ediffkzneg,kxy,xy
     $     ,nside,planf ,planb ,plan2f ,plan2b,nmat,file_id ,group_idmic
     $     ,nstop ,infostr)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
c     variables en argument      
      integer ntotalm,ntotal,ldabi,nlar,nmax,nxm,nym,nzm,nx,ny,nz ,nx2
     $     ,ny2,nxy2,nz2,nbsphere,nbsphere3,nloop,ncompte,nstop,nfft2d
     $     ,ndipole,nproche,nrig,nmat,ipol,npol,npolainc,nquicklens
     $     ,nside,niter,niterii
      double precision tol,tolinit,tol1,aretecube,eps0,k0,k02,P0,irra,w0
     $     ,I0,numaper,numaperinc,gross,deltax,zlens
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xs,ys,zs
      double complex, dimension(8*nxm*nym*nzm) :: FFTTENSORxx,
     $     FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz, FFTTENSORzz
     $     ,vectx,vecty,vectz
      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex E0,Eloc(3),Em(3),epsani(3,3)
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Efourierx(nfft2d*nfft2d)
     $     ,Efouriery(nfft2d*nfft2d) ,Efourierz(nfft2d*nfft2d)
     $     ,Eimageincx(nfft2d*nfft2d),Eimageincy(nfft2d*nfft2d)
     $     ,Eimageincz(nfft2d*nfft2d),Efourierincx(nfft2d*nfft2d)
     $     ,Efourierincy(nfft2d*nfft2d),Efourierincz(nfft2d*nfft2d)
      double complex Ediffkzpos(nfft2d,nfft2d,3),Ediffkzneg(nfft2d
     $     ,nfft2d,3)
      double precision kxy(nfft2d),xy(nfft2d)
      integer, dimension(nxm*nym*nzm) :: Tabdip
      character(64) infostr,beam
      character(12) methodeit
      
      integer*8 planf,planb,plan2f,plan2b,planfn,planbn
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD

c     variable pour le bright field
      integer i,j,ii,jj,k,kk,ideltam,idelta,jdelta,nsens,imaxk0,indice
     $     ,indicex,indicey,nfft2d2,ikxinc,jkyinc,imul,imaxinc
      double precision deltak,numaperk,phi,theta,zero,pi,kx,ky,kz,ss,pp
     $     ,u1,u2,tmp,normal(3),xmin,xmax,ymin,ymax,deltakx,deltaky
     $     ,deltakxy,kxinc,kyinc,u(3),v(3),sintmp,costmp,deltakm
      double complex tmpx,tmpy,tmpz,Ex,Ey,Ez,Emx,Emy,Emz,ctmp,ctmp1
     $     ,icomp,zfocus

      character(LEN=100) :: datasetname
#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif

      integer(hid_t) :: file_id
      integer(hid_t) :: group_idmic
      integer :: dim(4)
      integer error


      write(*,*) '***** Bright field microscope ******'

c     initialise
      zero=0.d0
      pi=dacos(-1.d0)
      beam='pwavelinear'
      nfft2d2=nfft2d/2
      npolainc=0
      icomp=(0.d0,1.d0)
      numaperk=k0*numaperinc
      k02=k0*k0
      if (numaperinc.ge.1.d0.or.numaper.le.0.d0) then
         infostr='NA inc strictly between 0 and 1'
         nstop=1
         return
      endif

      
c     calcul de deltak      
      xmax=-1.d300
      xmin=1.d300
      ymax=-1.d300
      ymin=1.d300
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:xmax,ymax)
!$OMP& REDUCTION(min:xmin,ymin)      
      do i=1,nbsphere
         xmax=max(xmax,xs(i))
         xmin=min(xmin,xs(i))
         ymax=max(ymax,ys(i))
         ymin=min(ymin,ys(i))     
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      deltakm=pi/max(xmax-xmin,ymax-ymin)
      deltax=aretecube
      
      if (nquicklens.eq.1) then
         
         deltakx=2.d0*pi/(aretecube*dble(nfft2d))

         if (deltakx.ge.numaperk) then
            nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 1'
            return
         endif        
         imul=idint(deltakm/deltakx)
         if (imul.eq.0) then
            nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 2'
            return            
         endif
 223     deltak=deltakx*dble(imul)
         write(*,*) 'change delta k incident:',deltak,'m-1',imul
         imaxinc=nint(k0/deltak)+1
         if (imaxinc.le.3) then
            imul=imul-1
            if (imul.eq.0) then
               nstop=1
          infostr='In FFT lens nfft2d too small:Increase size of FFT 3'
               return
            endif
            goto 223
         endif
         write(*,*) 'Step size delta k diffracted: ',deltakx,'m-1'
      else
         k=0
 224     deltakx=2.d0*pi/(dble(nfft2d)*aretecube)/dble(2**k)
         imaxk0=nint(k0/deltakx)+1            
         if (imaxk0.le.5) then
            k=k+1
            write(*,*) 'Change delta k diffracted:',deltakx,'m-1',k
            goto 224
         endif
         write(*,*) 'Final delta k diffracted',deltakx,'m-1'

         k=0
 222     deltak=deltakx*dnint(deltakm/deltakx)/dble(2**k)
         imaxinc=nint(k0/deltak)+1
         if (imaxinc.le.2) then
            k=k+1
            write(*,*) 'Change delta k incident:',deltak,'m-1',k
            goto 222
         endif
         write(*,*) 'Final delta k incident : ',deltak,'m-1'
      endif

      ideltam=imaxinc

      ii=0
      do idelta=-ideltam,ideltam
         do jdelta=-ideltam,ideltam
            kxinc=idelta*deltak
            kyinc=jdelta*deltak
            if (kxinc*kxinc+kyinc*kyinc.le.numaperk*numaperk) then
               ii=ii+1
            endif
         enddo
      enddo

      if (nfft2d.gt.4096) then
         nstop=1
         infostr='nfft2d too large'
         return
      endif
      if (deltak.ge.numaperk) then
         nstop=1
         infostr='In FFT lens nfft2d too small'
         return
      endif

c     initalise
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
      do i=1,nfft2d*nfft2d
         Eimageincx(i)=0.d0
         Eimageincy(i)=0.d0
         Eimageincz(i)=0.d0
         Eimagex(i)=0.d0
         Eimagey(i)=0.d0
         Eimagez(i)=0.d0                          
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
      
      npol=1
      niter=ii
      if (npolainc.eq.0) then
         npol=2
         niter=2*ii
      endif
c     calcul puissance
      P0=P0/dble(npol*ii)
c      P0=1.d0
c      npol=1
      call irradiance(P0,w0,E0,irra)
      I0=cdabs(E0)**2
      niterii=0
      do ipol=1,npol
       
         if (npolainc.eq.1) then
            ss=1.d0
            pp=0.d0
         elseif (npolainc.eq.2) then
            ss=0.d0
            pp=1.d0
         else
            if (ipol.eq.1) then
               ss=1.d0
               pp=0.d0
            else
               ss=0.d0
               pp=1.d0
            endif
         endif

         
c     sommation 
         do idelta=-ideltam,ideltam
            do jdelta=-ideltam,ideltam
c         do idelta=0,0
c            do jdelta=0,0
               kxinc=idelta*deltak
               kyinc=jdelta*deltak
               if (kxinc*kxinc+kyinc*kyinc.le.numaperk*numaperk) then
                  niterii=niterii+1
                  write(*,*) '*** incidence',niterii,'/',niter,' *****'
c     calcul champ incident
                  
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
                  do i=1,nbsphere           
                     call ondeplanekxky(xs(i),ys(i),zs(i),k0,E0,ss,pp
     $                    ,kxinc,kyinc,FF0(3*i-2),FF0(3*i-1),FF0(3*i)
     $                    ,nstop,infostr)
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     calcul champ local
                  if (nrig.eq.0) then
                     tol=tolinit
                     ncompte=0
                     nloop=0
                     if (nproche.eq.-1) then
                        call inverserig(FFTTENSORxx,FFTTENSORxy
     $                       ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz
     $                       ,FFTTENSORzz,vectx,vecty,vectz ,Tabdip
     $                       ,ntotalm,ntotal,ldabi,nlar,nmax,ndipole,nxm
     $                       ,nym ,nzm,nx,ny,nz,nx2,ny2,nxy2,nz2
     $                       ,nbsphere,nbsphere3,XI ,XR,wrk,FF,FF0,FFloc
     $                       ,polarisa,methodeit,tol,tol1,nloop,ncompte
     $                       ,planf,planb,nstop,infostr)
                        if (nstop.eq.1) return
                     else

                        call inverserigopt(FFTTENSORxx,FFTTENSORxy
     $                       ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz
     $                       ,FFTTENSORzz,vectx,vecty,vectz,ntotalm
     $                       ,ntotal ,ldabi,nlar,nmax,nxm,nym,nzm,nx,ny
     $                       ,nz,nx2,ny2 ,nxy2,nz2,nbsphere,nbsphere3,XI
     $                       ,XR,wrk,FF,FF0 ,FFloc,polarisa,methodeit
     $                       ,tol,tol1,nloop ,ncompte,planf,planb,nstop
     $                       ,infostr)
                        if (nstop.eq.1) return
                     endif
                  elseif (nrig.eq.1) then 
c     **************************************        
c     Renormalized Born approximation field
c     **************************************
                     write(*,*) 'Renormalized Born approximation'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere3
                        FFloc(i)=FF0(i)
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL          

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,jj)    
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere
                        k=3*(i-1)
                        do ii=1,3
                           FF(k+ii)=0.d0
                           do jj=1,3
                              FF(k+ii)=FF(k+ii)+polarisa(i,ii,jj)
     $                             *FFloc(k+jj)
                           enddo
                        enddo
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL
                  elseif (nrig.eq.2) then 
c     **************************************        
c     Born approximation field
c     **************************************
                     write(*,*) 'Born approximation'
                     nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,Em,Eloc,ii,jj,epsani)   
!$OMP DO SCHEDULE(STATIC) 
                     do k=1,nbsphere
                        kk=3*(k-1)
                        Em(1)=FF0(kk+1)
                        Em(2)=FF0(kk+2)
                        Em(3)=FF0(kk+3)
                        do ii=1,3
                           do jj=1,3
                              epsani(ii,jj)=epsilon(k,ii,jj)
                           enddo
                        enddo 
                        call local_macro(Eloc,Em,epsani,aretecube,k0
     $                       ,nsens)
                        FFloc(kk+1)=Eloc(1)
                        FFloc(kk+2)=Eloc(2)
                        FFloc(kk+3)=Eloc(3)
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL          

c     dipole
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,jj)    
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere
                        k=3*(i-1)
                        do ii=1,3
                           FF(k+ii)=0.d0
                           do jj=1,3
                              FF(k+ii)=FF(k+ii)+polarisa(i,ii,jj)
     $                             *FFloc(k+jj)
                           enddo
                        enddo
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL
                  elseif (nrig.eq.3) then
c     renormalized Born ordre 1
                     write(*,*) 'Born series order 1'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
                     do i=1,nbsphere3
                        xr(i)=FF0(i)           
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL
                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)  
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere
                        k=3*(i-1)
                        xi(k+1)=-polarisa(i,1,1)*xr(k+1)-polarisa(i,1,2)
     $                       *xr(k+2)-polarisa(i,1,3)*xr(k+3)
                        xi(k+2)=-polarisa(i,2,1)*xr(k+1)-polarisa(i,2,2)
     $                       *xr(k+2)-polarisa(i,2,3)*xr(k+3)
                        xi(k+3)=-polarisa(i,3,1)*xr(k+1)-polarisa(i,3,2)
     $                       *xr(k+2)-polarisa(i,3,3)*xr(k+3)
                     enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
                     
                     call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy
     $                    ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz
     $                    ,FFTTENSORzz ,vectx,vecty,vectz,Tabdip,ntotalm
     $                    ,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny,nz
     $                    ,nx2,ny2,nxy2,nz2 ,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)       
                     do i=1,nbsphere
                        k=3*(i-1)
                        FFloc(k+1)=xr(k+1)
                        FFloc(k+2)=xr(k+2)
                        FFloc(k+3)=xr(k+3)
                        FF(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)
     $                       *xr(k+2)+polarisa(i,1,3)*xr(k+3)
                        FF(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)
     $                       *xr(k+2)+polarisa(i,2,3)*xr(k+3)
                        FF(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)
     $                       *xr(k+2)+polarisa(i,3,3)*xr(k+3)
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      

                     
                  elseif (nrig.eq.4) then
c     Rytov renormalize
                     write(*,*) 'Rytov approximation on local field'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere3
                        xr(i)=0.d0       
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL
                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere
                        k=3*(i-1)
                        xi(k+1)=-polarisa(i,1,1)*FF0(k+1)-polarisa(i,1
     $                       ,2)*FF0(k+2)-polarisa(i,1,3)*FF0(k+3)
                        xi(k+2)=-polarisa(i,2,1)*FF0(k+1)-polarisa(i,2
     $                       ,2)*FF0(k+2)-polarisa(i,2,3)*FF0(k+3)
                        xi(k+3)=-polarisa(i,3,1)*FF0(k+1)-polarisa(i,3
     $                       ,2)*FF0(k+2)-polarisa(i,3,3)*FF0(k+3)
                     enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
                     
                     call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy
     $                    ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz
     $                    ,FFTTENSORzz ,vectx,vecty,vectz,Tabdip,ntotalm
     $                    ,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny,nz
     $                    ,nx2,ny2,nxy2,nz2 ,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp)
!$OMP DO SCHEDULE(STATIC)      
                     do i=1,nbsphere
                        k=3*(i-1)
                        tmp=dsqrt(dreal(FF0(k+1)*dconjg(FF0(k+1))+FF0(k
     $                       +2)*dconjg(FF0(k+2))+FF0(k+3)*dconjg(FF0(k
     $                       +3))))*1.d-6
                        if (cdabs(FF0(k+1)).le.tmp) then
                           FFloc(k+1)=0.d0
                        else
                           FFloc(k+1)=FF0(k+1)*cdexp(xr(k+1)/FF0(k+1))
                        endif
                        if (cdabs(FF0(k+2)).le.tmp) then
                           FFloc(k+2)=0.d0
                        else
                           FFloc(k+2)=FF0(k+2)*cdexp(xr(k+2)/FF0(k+2))
                        endif
                        if (cdabs(FF0(k+3)).le.tmp) then
                           FFloc(k+3)=0.d0
                        else
                           FFloc(k+3)=FF0(k+3)*cdexp(xr(k+3)/FF0(k+3))
                        endif
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
                     do i=1,nbsphere
                        k=3*(i-1)  
                        FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1
     $                       ,2)*FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
                        FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2
     $                       ,2)*FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
                        FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3
     $                       ,2)*FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3)
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      

                  elseif (nrig.eq.5) then
c     Rytov
                     write(*,*)
     $                    'Rytov approximation on macroscopic field'
c     a3=aretecube*aretecube*aretecube
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC) 
                     do i=1,nbsphere
                        k=3*(i-1)
                        xr(k+1)=-FF0(k+1)*(epsilon(i,1,1)-eps0)/3.d0
                        xr(k+2)=-FF0(k+2)*(epsilon(i,2,2)-eps0)/3.d0
                        xr(k+3)=-FF0(k+3)*(epsilon(i,3,3)-eps0)/3.d0
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL
                     nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,Em,Eloc,ii,jj,epsani)
!$OMP DO SCHEDULE(STATIC)
                     do i=1,nbsphere
                        k=3*(i-1)
                        Em(1)=FF0(k+1)
                        Em(2)=FF0(k+2)
                        Em(3)=FF0(k+3)
                        do ii=1,3
                           do jj=1,3
                              epsani(ii,jj)=epsilon(i,ii,jj)
                           enddo
                        enddo 
                        call local_macro(Eloc,Em,epsani,aretecube,k0
     $                       ,nsens)
                        FF(k+1)=Eloc(1)
                        FF(k+2)=Eloc(2)
                        FF(k+3)=Eloc(3)
                        xi(k+1)=-polarisa(i,1,1)*FF(k+1)-polarisa(i,1,2)
     $                       *FF(k+2)-polarisa(i,1,3)*FF(k+3)
                        xi(k+2)=-polarisa(i,2,1)*FF(k+1)-polarisa(i,2,2)
     $                       *FF(k+2)-polarisa(i,2,3)*FF(k+3)
                        xi(k+3)=-polarisa(i,3,1)*FF(k+1)-polarisa(i,3,2)
     $                       *FF(k+2)-polarisa(i,3,3)*FF(k+3)
                     enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
                     
                     call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy
     $                    ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz
     $                    ,FFTTENSORzz ,vectx,vecty,vectz,Tabdip,ntotalm
     $                    ,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny,nz
     $                    ,nx2,ny2,nxy2,nz2 ,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp)
!$OMP DO SCHEDULE(STATIC)      
                     do i=1,nbsphere
                        k=3*(i-1)
                        tmp=dsqrt(dreal(FF(k+1)*dconjg(FF(k+1))+FF(k+2)
     $                       *dconjg(FF(k+2))+FF(k+3)*dconjg(FF(k+3))))
     $                       *1.d -6
                        if (cdabs(FF(k+1)).le.tmp) then
                           FFloc(k+1)=0.d0
                        else
                           FFloc(k+1)=FF(k+1)*cdexp(xr(k+1)/FF(k+1))
                        endif
                        if (cdabs(FF(k+2)).le.tmp) then
                           FFloc(k+2)=0.d0
                        else
                           FFloc(k+2)=FF(k+2)*cdexp(xr(k+2)/FF(k+2))
                        endif
                        if (cdabs(FF(k+3)).le.tmp) then
                           FFloc(k+3)=0.d0
                        else
                           FFloc(k+3)=FF(k+3)*cdexp(xr(k+3)/FF(k+3))
                        endif
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)      
                     do i=1,nbsphere
                        k=3*(i-1)  
                        FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1
     $                       ,2)*FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
                        FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2
     $                       ,2)*FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
                        FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3
     $                       ,2)*FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3)
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      
                     
                  elseif (nrig.eq.6) then
c     Beam propagation method
                     write(*,*) 'Beam propagation method'

                     call beampropagationmacro(xs,ys,zs,aretecube,k0,w0
     $                    ,E0,ss,pp,theta,phi,zero,zero,zero,beam
     $                    ,epsilon,ndipole ,nx ,ny,nz,nxm,nym ,nzm,nmax
     $                    ,nfft2d,Eimagex,Eimagey,Eimagez,FF0,FFloc,FF
     $                    ,plan2f,plan2b ,nstop,infostr)

                     if (nstop.eq.1) return

                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
                     do i=1,nbsphere
                        k=3*(i-1)  
                        FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1
     $                       ,2)*FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
                        FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2
     $                       ,2)*FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
                        FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3
     $                       ,2)*FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3)
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL    
                  elseif (nrig.eq.7) then
c     Beam propagation method
                     call beampropagation(xs,ys,zs,aretecube,k0,w0,E0,ss
     $                    ,pp,theta,phi,zero,zero,zero,beam,epsilon
     $                    ,ndipole,nx ,ny,nz,nxm,nym ,nzm,nmax,nfft2d
     $                    ,Eimagex,Eimagey,Eimagez,FF0 ,FFloc,FF,plan2f
     $                    ,plan2b,nstop,infostr)

                     if (nstop.eq.1) return

                     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
                     do i=1,nbsphere
                        k=3*(i-1)  
                        FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1
     $                       ,2)*FFloc(k+2)+polarisa(i,1,3)*FFloc(k+3)
                        FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2
     $                       ,2)*FFloc(k+2)+polarisa(i,2,3)*FFloc(k+3)
                        FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3
     $                       ,2)*FFloc(k+2)+polarisa(i,3,3)*FFloc(k+3)
                     enddo            
!$OMP ENDDO 
!$OMP END PARALLEL    
                  endif
c     *********************************************************
c     fin calcul champ local
c     *********************************************************
c     *********************************************************            
c     calcul champ diffracte
c     *********************************************************
c                  write(*,*) 'fff local',FF
                  if (nquicklens.eq.1) then
c                    write(*,*) 'fff local',nx,ny,nz,nxm,nym,nzm,nfft2d
c    $                    ,k0,imaxk0,deltakx,deltaky
                     call diffractefft2dlens(nx,ny,nz,nxm,nym,nzm,nfft2d
     $                    ,k0,xs,ys,zs,aretecube,Efourierx,Efouriery
     $                    ,Efourierz,FF,imaxk0,deltakx,deltaky
     $                    ,Ediffkzpos,numaper,nside,plan2f,plan2b ,nstop
     $                    ,infostr)
c                     write(*,*) 'fff diff',Ediffkzpos
                     if (nstop.eq.1) return
                  else
                     
                     
                     deltaky=deltakx
                     do i=-imaxk0,imaxk0
                        kx=deltakx*dble(i)
                        do j=-imaxk0,imaxk0
                           ky=deltaky*dble(j)
                           
                           if (dsqrt(kx*kx+ky*ky).le.numaper) then
                              kz=dsqrt(k02-kx*kx-ky*ky)
                              normal(1)=kx/k0
                              normal(2)=ky/k0
                              normal(3)=dsqrt(1.d0-normal(1)*normal(1)
     $                             -normal(2)*normal(2))
                              Emx=0.d0
                              Emy=0.d0
                              Emz=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,kk,ctmp,ctmp1)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Emx,Emy,Emz)
                              do ii=1,nbsphere
                                 kk=3*(ii-1)
                                 ctmp=cdexp(-icomp*k0*(normal(1)*xs(ii)
     $                                +normal(2)*ys(ii)+normal(3)
     $                                *zs(ii)))
                                 ctmp1=normal(1)*FF(kk+1)+normal(2)
     $                                *FF(kk+2)+normal(3)*FF(kk+3)
                                 Emx=Emx+ctmp*(FF(kk+1)-ctmp1*normal(1))
                                 Emy=Emy+ctmp*(FF(kk+2)-ctmp1*normal(2))
                                 Emz=Emz+ctmp*(FF(kk+3)-ctmp1*normal(3))
                              enddo
!$OMP ENDDO 
!$OMP END PARALLEL                     
                              ctmp=-2.d0*pi*icomp*kz                   
                              kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                              ii=imaxk0+i+1
                              jj=imaxk0+j+1
                              Ediffkzpos(ii,jj,1)=Emx*k02/ctmp
                              Ediffkzpos(ii,jj,2)=Emy*k02/ctmp
                              Ediffkzpos(ii,jj,3)=Emz*k02/ctmp
                           endif
                        enddo
                     enddo
                  endif   
c     *********************************************************            
c     fin champ diffracte
c     *********************************************************
c     *********************************************************            
c     calcul image
c     *********************************************************  
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
                  do i=1,nfft2d*nfft2d
                     Efourierincx(i)=0.d0
                     Efourierincy(i)=0.d0
                     Efourierincz(i)=0.d0
                     Efourierx(i)=0.d0
                     Efouriery(i)=0.d0
                     Efourierz(i)=0.d0                          
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL

                  call deltakroutine(kxinc,kyinc,deltak,deltak,k0
     $                 ,ikxinc,jkyinc)
                  deltakxy=deltakx*deltaky
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice)
!$OMP& PRIVATE(kx,ky,kz,u1,u2,normal,tmp,tmpx,tmpy,tmpz,ii,jj)
!$OMP& PRIVATE(phi,theta,Ex,Ey,Ez,zfocus,u,v,costmp,sintmp)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
                  do i=-imaxk0,imaxk0 
                     do j=-imaxk0,imaxk0
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

                        kx=deltakx*dble(i)
                        ky=deltaky*dble(j)
                        tmp=kx*kx+ky*ky
                        if (tmp.le.numaper*numaper) then
                           kz=dsqrt(k02-tmp)
                           zfocus=cdexp(icomp*kz*zlens)

                           ii=imaxk0+i+1
                           jj=imaxk0+j+1
                           Efourierx(indice)=Ediffkzpos(ii,jj,1)*zfocus
                           Efouriery(indice)=Ediffkzpos(ii,jj,2)*zfocus
                           Efourierz(indice)=Ediffkzpos(ii,jj,3)*zfocus
                           Efourierincx(indice)=Ediffkzpos(ii,jj,1)
     $                          *zfocus
                           Efourierincy(indice)=Ediffkzpos(ii,jj,2)
     $                          *zfocus
                        
                           Efourierincz(indice)=Ediffkzpos(ii,jj,3)
     $                          *zfocus

                           if (i.eq.ikxinc.and.j.eq.jkyinc) then

                              call   ondeplanekxky(0.d0,0.d0,0.d0,k0,E0
     $                             ,ss,pp,kx,ky,Ex,Ey,Ez,nstop,infostr)
                              Efourierincx(indice)=Efourierincx(indice)
     $                             +Ex/deltakxy*zfocus
                              Efourierincy(indice)=Efourierincy(indice)
     $                             +Ey/deltakxy*zfocus
                              Efourierincz(indice)=Efourierincz(indice)
     $                             +Ez/deltaky*zfocus
                           endif
c                           write(*,*) 'fff diff',Efourierincy(indice)
c     $                          ,Ey/deltakx/deltaky*icomp,indice

                           if (tmp.ne.0.d0) then      

                              u(1)=kx/k0
                              u(2)=ky/k0
                              u(3)=dsqrt(1.d0-u(1)*u(1)-u(2)*u(2))
                              v(1)=-kx/k0/gross
                              v(2)=-ky/k0/gross
                              v(3)=dsqrt(1.d0-v(1)*v(1)-v(2)*v(2))
                              u1=u(2)*v(3)-u(3)*v(2)
                              u2=-u(1)*v(3)+u(3)*v(1)
                              costmp=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
                              sintmp=dsqrt(u1*u1+u2*u2)
                              u1=u1/sintmp
                              u2=u2/sintmp
                              tmp=dsqrt(u(3)/v(3))

                              tmpx=(u1*u1+(1.d0-u1*u1)*costmp)
     $                             *Efourierincx(indice)+u1*u2*(1.d0
     $                             -costmp)*Efourierincy(indice)+u2
     $                             *sintmp*Efourierincz(indice)
                              tmpy=u1*u2*(1.d0-costmp)
     $                             *Efourierincx(indice)+(u2*u2+(1.d0-u2
     $                             *u2)*costmp)*Efourierincy(indice) -u1
     $                             *sintmp*Efourierincz(indice)
                              tmpz=-u2*sintmp*Efourierincx(indice)+u1
     $                             *sintmp*Efourierincy(indice)+costmp
     $                             *Efourierincz(indice)
                              Efourierincx(indice)=tmpx*tmp
                              Efourierincy(indice)=tmpy*tmp
                              Efourierincz(indice)=tmpz*tmp
                              tmpx=(u1*u1+(1.d0-u1*u1)*costmp)
     $                             *Efourierx(indice) +u1*u2*(1.d0
     $                             -costmp)*Efouriery(indice)+u2 *sintmp
     $                             *Efourierz(indice)
                              tmpy=u1*u2*(1.d0-costmp) *Efourierx(indice
     $                             )+(u2*u2+(1.d0-u2 *u2)*costmp)
     $                             *Efouriery(indice)-u1*sintmp
     $                             *Efourierz(indice)
                              tmpz=-u2*sintmp*Efourierx(indice)+u1
     $                             *sintmp*Efouriery(indice)+costmp
     $                             *Efourierz(indice)
                              Efourierx(indice)=tmpx*tmp
                              Efouriery(indice)=tmpy*tmp
                              Efourierz(indice)=tmpz*tmp
                           endif
                        endif
                     enddo
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL
c                  do i=1,nfft2D*nfft2D
c                     write(998,*) 'f',Efouriery(i),i
c                  enddo
                  call fouriertoimage(deltakx,deltaky,gross,Efourierx
     $                 ,Efouriery,Efourierz,Efourierincx,Efourierincy
     $                 ,Efourierincz,nfft2D ,nfft2d2 ,plan2b,plan2f)
c                  do i=1,nfft2D*nfft2D
c                     write(998,*) 'f',Efouriery(i),i
c                  enddo
c                  write(*,*) 'fff Eimage',Efouriery
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice)   
!$OMP DO SCHEDULE(STATIC)    
                  do indice=1,nfft2d*nfft2d
                     Eimagex(indice)=Efourierx(indice)
     $                    *dconjg(Efourierx(indice))+Eimagex(indice)
                     Eimagey(indice)=Efouriery(indice)
     $                    *dconjg(Efouriery(indice))+Eimagey(indice)
                     Eimagez(indice)=Efourierz(indice)
     $                    *dconjg(Efourierz(indice))+Eimagez(indice)
                     Eimageincx(indice)=Efourierincx(indice)
     $                    *dconjg(Efourierincx(indice))
     $                    +Eimageincx(indice)
                     Eimageincy(indice)=Efourierincy(indice)
     $                    *dconjg(Efourierincy(indice))
     $                    +Eimageincy(indice)
                     Eimageincz(indice)=Efourierincz(indice)
     $                    *dconjg(Efourierincz(indice))
     $                    +Eimageincz(indice)
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               endif

               if (nstop == -1) then
                  infostr =
     $                 'Calculation cancelled during iterative method'
                  return
               endif
            enddo
         enddo
      enddo

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice)   
!$OMP DO SCHEDULE(STATIC)        
      do indice=1,nfft2d*nfft2d
         Eimagex(indice)=cdsqrt(Eimagex(indice))
         Eimagey(indice)=cdsqrt(Eimagey(indice))
         Eimagez(indice)=cdsqrt(Eimagez(indice))
         Eimageincx(indice)=cdsqrt(Eimageincx(indice))
         Eimageincy(indice)=cdsqrt(Eimageincy(indice))
         Eimageincz(indice)=cdsqrt(Eimageincz(indice))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      deltax=2.d0*pi/dble(nfft2d)/deltakx

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)   
      do i=-nfft2d2,nfft2d2-1
         xy(i+nfft2d2+1)=deltax*dble(i)*gross
         kxy(i+nfft2d2+1)=deltakx*dble(i)/k0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      
      if (nmat.eq.0) then
         open(400,file='kxfourier.mat')
         open(401,file='ximage.mat')
         do i=-nfft2d2,nfft2d2-1
            write(401,*) xy(i+nfft2d2+1)
         enddo
         do i=-imaxk0,imaxk0
            kx=deltakx*dble(i)
            write(400,*) kx/k0
         enddo
         close(400)
         close(401)

         open(310,file='kxincidentbf.mat')
         open(311,file='kyincidentbf.mat')
         do idelta=-ideltam,ideltam
            do jdelta=-ideltam,ideltam
               kxinc=idelta*deltak
               kyinc=jdelta*deltak
               if (kxinc*kxinc+kyinc*kyinc.le.numaperk*numaperk) then
                  write(310,*) kxinc/k0
                  write(311,*) kyinc/k0
               endif
            enddo
         enddo
         close (310)
         close (311)

         open(301,file='imagebf.mat')
         open(302,file='imagebfx.mat')
         open(303,file='imagebfy.mat')
         open(304,file='imagebfz.mat')
         open(305,file='imageincbf.mat')
         open(306,file='imageincbfx.mat')
         open(307,file='imageincbfy.mat')
         open(308,file='imageincbfz.mat')

         do i=1,nfft2D*nfft2D
            write(301,*) dsqrt(dreal(Eimagex(i)**2.d0+Eimagey(i)**2.d0
     $           +Eimagez(i)**2.d0))
            write(302,*) dreal(Eimagex(i))
            write(303,*) dreal(Eimagey(i))
            write(304,*) dreal(Eimagez(i))
            write(305,*) dsqrt(dreal(Eimageincx(i)**2.d0+Eimageincy(i)
     $           **2.d0+Eimageincz(i)**2.d0))
            write(306,*) dreal(Eimageincx(i))
            write(307,*) dreal(Eimageincy(i))
            write(308,*) dreal(Eimageincz(i))
         enddo
         close(301)
         close(302)
         close(303)
         close(304)
         close(305)
         close(306)
         close(307)
         close(308)
      elseif (nmat.eq.2) then
         dim(1)=nfft2d
         dim(2)=nfft2d
         datasetname='x Image'
         call hdf5write1d(group_idmic,datasetname,xy,dim)
         k=0
         datasetname='Image bright field kz>0'
         call writehdf5mic(Eimagex,Eimagey,Eimagez,nfft2d ,imaxk0
     $        ,Ediffkzpos,k,datasetname,group_idmic)
         datasetname='Image+incident bright field kz>0'
         call writehdf5mic(Eimageincx,Eimageincy,Eimageincz
     $        ,nfft2d,imaxk0,Ediffkzpos,k,datasetname,group_idmic)

         ii=0
         do idelta=-ideltam,ideltam
            do jdelta=-ideltam,ideltam
               kxinc=idelta*deltak
               kyinc=jdelta*deltak
               if (kxinc*kxinc+kyinc*kyinc.le.numaperk*numaperk) then
                  ii=ii+1
                  kxy(ii)=kxinc/k0
               endif
            enddo
         enddo
         dim(1)=ii
         dim(2)=nfft2d
         datasetname='kx incident bf'
         call hdf5write1d(group_idmic,datasetname,kxy,dim)

         ii=0
         do idelta=-ideltam,ideltam
            do jdelta=-ideltam,ideltam
               kxinc=idelta*deltak
               kyinc=jdelta*deltak
               if (kxinc*kxinc+kyinc*kyinc.le.numaperk*numaperk) then
                  ii=ii+1
                  kxy(ii)=kyinc/k0
               endif
            enddo
         enddo
         dim(1)=ii
         dim(2)=nfft2d
         datasetname='ky incident bf'
         call hdf5write1d(group_idmic,datasetname,kxy,dim)
      endif


      end
