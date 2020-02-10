      subroutine objetsphereinhomo(eps,eps0,xs,ys,zs,k0 ,aretecube
     $     ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,methode ,na
     $     ,epsilon,polarisa,rayon,lc,hc,ng,epsb,nmat,file_id
     $     ,group_iddip,infostr,nstop)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,i,j,k,test
     $     ,IP(3),nnnr,dddis,inv,na,nstop,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,lc,hc,x,y,z
     $     ,eps0,aretecube,pi
      double complex eps,polarisa(nmax,3,3) ,epsilon(nmax,3 ,3),ctmp
     $     ,epsb(nmax),icomp

      integer nk,ngraine,ng,iret,omp_get_max_threads
      double precision rayon,kx,ky,kz,dkx,dky,dkz,coeff1,coeff2,coeff3
     $     ,phase,spec,moyenne,ecartype,lx,ly,lz,sunif
      character*2 methode
      character(64) infostr
      integer*8 planb
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD

      character(LEN=100) :: datasetname

#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif

      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error

      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64

      if (lc.le.0.d0) then
         infostr='coherence length equal to zero'
         nstop=1
         return
      endif      
      if (hc.le.0.d0) then
         infostr='standard deviation equal to zero'
         nstop=1
         return
      endif
      
c     Initialization
      nbsphere=0
      ndipole=0 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
      do i=1,nmax
         Tabdip(i)=0
         polarisa(i,1,1)=0.d0
         polarisa(i,1,2)=0.d0
         polarisa(i,1,3)=0.d0
         polarisa(i,2,1)=0.d0
         polarisa(i,2,2)=0.d0
         polarisa(i,2,3)=0.d0
         polarisa(i,3,1)=0.d0
         polarisa(i,3,2)=0.d0
         polarisa(i,3,3)=0.d0
         epsilon(i,1,1)=0.d0
         epsilon(i,1,2)=0.d0
         epsilon(i,1,3)=0.d0
         epsilon(i,2,1)=0.d0
         epsilon(i,2,2)=0.d0
         epsilon(i,2,3)=0.d0
         epsilon(i,3,1)=0.d0
         epsilon(i,3,2)=0.d0
         epsilon(i,3,3)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      dddis=1
      inv=1
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      ngraine=4*ng+1

      lc=lc*1.d-9
      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      rayon=rayon*1.d-9

      if (rayon.le.1.d-12) then
         nstop=1
         infostr='object sphere: radius=0!'
         return
      endif
c      write (*,*) 'na = ',na
c      write (*,*) 'nnnr = ',nnnr

c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.

      nx=nnnr
      ny=nnnr
      nz=nnnr
c     size of the subunit
      aretecube=2.d0*rayon/dble(nnnr)

#ifdef USE_FFTW
      call dfftw_plan_dft_3d(planb, nx,ny,nz,epsb,epsb,FFTW_BACKWARD
     $     ,FFTW_ESTIMATE)
#else 
      call fftsingletonz3d(epsb, nx,ny,nz,FFTW_BACKWARD)
#endif

      lx=dble(nx)*aretecube
      ly=dble(ny)*aretecube
      lz=dble(nz)*aretecube
      dkx=2.d0*pi/lx
      dky=2.d0*pi/ly
      dkz=2.d0*pi/lz
      nk=0
      
      coeff1=0.125d0*hc*hc*lc*lc*lc/(pi*dsqrt(pi))
      coeff2=lc*lc
      coeff3=4.d0*dkx*dky*dkz
      do i=1,nz
         do j=1,ny
            do k=1,nx
               nk=nk+1
               if (k.ge.1 .and. k.le.(nx/2+1)) then
                  kx=dble(k-1)*dkx
               else
                  kx=dble(k-1-nx)*dkx
               endif
               if (j.ge.1 .and. j.le.(ny/2+1)) then
                  ky=dble(j-1)*dky
               else
                  ky=dble(j-1-ny)*dky
               endif
               if (i.ge.1 .and. i.le.(nz/2+1)) then
                  kz=dble(i-1)*dkz
               else
                  kz=dble(i-1-nz)*dkz
               endif
c     calcul du spec
               spec=coeff1*dexp(-0.25d0*coeff2*(kx*kx+ky*ky+kz*kz))
c     phase aleatoire
               phase=2.d0*pi*SUNIF(ngraine)
c     Amplitude complexe
               epsb(nk)=dsqrt(spec*coeff3)*cdexp(icomp *phase)
c     Spectre symetrique pour diffraction harmonique
               if (nk.eq.1) epsb(nk)=0.d0
               if (nk.ge.nx/2+2.and.nk.le.nx) epsb(nk)=0.d0
               if (nk.ge.(ny/2+1)*nx+1.and.nk.le.nx*ny) epsb(nk)=0.d0
               if (nk.ge.(nz/2+1)*nx*ny+1.and.nk.le.nx*ny*nz) epsb(nk)
     $              =0.d0

            enddo
         enddo
      enddo     
      
c     Profil des hauteurs
#ifdef USE_FFTW
      call dfftw_execute_dft(planb, epsb, epsb)   
#endif
      
      moyenne=0.d0
      ecartype=0.d0
      do i=1,nk
         epsb(i)=dreal(epsb(i))+eps
         moyenne=moyenne+dreal(epsb(i))
         ecartype=ecartype+cdabs(epsb(i))**2.d0      
      enddo
        
      
      moyenne=moyenne/dble(nk)
      ecartype=ecartype/dble(nk)
    
      
      if (na.eq.-1) then
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x

                  ndipole=ndipole+1
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then
                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z
                     eps=epsb(ndipole)
                     call poladiff(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)  
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                  endif
               enddo
            enddo
         enddo
      elseif (na.eq.0) then
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then                    
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z
                     eps=epsb(ndipole)

                     call poladiff(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)  
c     write (*,*) 'nbsphere = ', nbsphere
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                  else
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z
                     epsilon(nbsphere,1,1)=(1.d0,0.d0)
                     epsilon(nbsphere,2,2)=(1.d0,0.d0)
                     epsilon(nbsphere,3,3)=(1.d0,0.d0)
                  endif
               enddo
            enddo
         enddo
      else
         infostr='na should be equal to -1 or 0'
         nstop=1
         return
      endif
      if (ndipole.gt.nmax) then
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif

      if (nmat.eq.0) then
         do i=1,nbsphere
            write(10,*) xs(i)
            write(11,*) ys(i)
            write(12,*) zs(i)
         enddo
      elseif (nmat.eq.2) then
         dim(1)=nbsphere
         dim(2)=nmax
         datasetname="Dipole position x"
         call hdf5write1d(group_iddip,datasetname,xs,dim)
         datasetname="Dipole position y"
         call hdf5write1d(group_iddip,datasetname,ys,dim)
         datasetname="Dipole position z"
         call hdf5write1d(group_iddip,datasetname,zs,dim)
      endif


      close(10)
      close(11)
      close(12)
      close(20)
      close(21)
      close(22)
      write (*,*) ' OBJECT SPHERE INHOMOGENEOUS FINISHED'
      write (*,*) ' SEED',ng,'Average',moyenne,'Standard deviation'
     $     ,dsqrt(ecartype-moyenne*moyenne)

      write (99,*) ' SEED',ng
      write(99,*) 'Average',moyenne
      write(99,*) 'Standard deviation',dsqrt(ecartype-moyenne*moyenne)

      end
