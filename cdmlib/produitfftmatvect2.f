c     FFB=FFB-A*FFX==> FFB=XI et FFX=XI*pola==> FFB=(I-A*pola)*XI
      subroutine produitfftmatvect2(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $     ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz,vectbx
     $     ,vectby,vectbz,FF,Tabdip,ntotalm,ntotal ,nmax,ndipole,nxm,nym
     $     ,nzm ,nx,ny,nz,nx2,ny2,nxy2,nz2,FFX,FFB)
      implicit none
      integer ii,jj,i,j,k,indice,ntotal,ntotalm,nxm,nym,nzm,nx,ny,nz
     $     ,nmax,nx2,ny2,nxy2,nz2,ndipole,Tabdip(nxm*nym*nzm)
      double complex FFTTENSORxx(ntotalm),FFTTENSORxy(ntotalm)
     $     ,FFTTENSORxz(ntotalm),FFTTENSORyy(ntotalm)
     $     ,FFTTENSORyz(ntotalm),FFTTENSORzz(ntotalm),FFX(3*nmax),FFB(3
     $     *nmax),vectx(ntotalm) ,vecty(ntotalm),vectz(ntotalm)
     $     ,vectbx(ntotalm),vectby(ntotalm),vectbz(ntotalm),FF(ntotalm)

c     calcul FFT du vecteur B
c      write(*,*) 'ttt',Tabdip
      do k=1,nz
         do j=1,ny
            do i=1,nx
c     position du dipole
               ii=i+nx*(j-1)+nx*ny*(k-1)
               indice=i+nx2*(j-1)+nxy2*(k-1)
c               write(*,*) 'rrr',Tabdip(ii)
               if (Tabdip(ii).eq.0) then
                  vectx(indice)=0.d0
                  vecty(indice)=0.d0
                  vectz(indice)=0.d0
               else
                  jj=3*Tabdip(ii)
                  vectx(indice)=FFX(jj-2)
                  vecty(indice)=FFX(jj-1)
                  vectz(indice)=FFX(jj)
c                  write(*,*) 'FF',FFX(jj-2),jj,ii,Tabdip(ii)
               endif
               vectx(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
               vectx(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vectx(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vectx(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vectx(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vectx(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
               vectx(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
               
               vecty(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
               vecty(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vecty(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vecty(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vecty(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vecty(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
               vecty(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
               
               vectz(i+nx+nx2*(j-1)+nxy2*(k-1))=0.d0
               vectz(i+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vectz(i+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vectz(i+nx+nx2*(j+ny-1)+nxy2*(k-1))=0.d0
               vectz(i+nx+nx2*(j-1)+nxy2*(k+nz-1))=0.d0
               vectz(i+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
               vectz(i+nx+nx2*(j+ny-1)+nxy2*(k+nz-1))=0.d0
            enddo
         enddo
      enddo
      CALL ZFFT3D(vectx,NX2,NY2,NZ2,1)
      CALL ZFFT3D(vecty,NX2,NY2,NZ2,1)
      CALL ZFFT3D(vectz,NX2,NY2,NZ2,1)
      
c     produit des TF
      do indice=1,ntotal
         vectbx(indice)=FFTTENSORxx(indice)*vectx(indice)+
     *        FFTTENSORxy(indice)*vecty(indice)+
     *        FFTTENSORxz(indice)*vectz(indice)
         
         vectby(indice)=FFTTENSORxy(indice)*vectx(indice)+
     *        FFTTENSORyy(indice)*vecty(indice)+
     *        FFTTENSORyz(indice)*vectz(indice)
         
         vectbz(indice)=FFTTENSORxz(indice)*vectx(indice)+
     *        FFTTENSORyz(indice)*vecty(indice)+
     *        FFTTENSORzz(indice)*vectz(indice)
      enddo
      
c     FFT inverse (deconvolution)
      CALL ZFFT3D(vectbx,nx2,ny2,nz2,-1)
      CALL ZFFT3D(vectby,nx2,ny2,nz2,-1)
      CALL ZFFT3D(vectbz,nx2,ny2,nz2,-1)
      
c     write(*,*) 'COLB',vecs(3*nbsphere,COLB)
c     remet le vecteur resultat dans FF
      do k=1,nz
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)+nxy2*(k-1)
               ii=3*(i+nx*(j-1)+nx*ny*(k-1))
               FF(ii-2)=-vectbx(indice)
               FF(ii-1)=-vectby(indice)
               FF(ii)=-vectbz(indice)
            enddo
         enddo
      enddo
      
c     calcul le B final en tenant compte de la diago faite au debut
c     et sans oblier le fac
      do i=1,ndipole
         ii=3*i
         if (Tabdip(i).ne.0) then
            jj=3*Tabdip(i)
            FFB(jj-2)=FFB(jj-2)+FF(ii-2)
            FFB(jj-1)=FFB(jj-1)+FF(ii-1)
            FFB(jj)=FFB(jj)+FF(ii)
         endif
      enddo

      end
   

c*************************************************
