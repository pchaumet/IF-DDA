      subroutine objetarbitrary(trope,eps,epsani,eps0,xs,ys,zs,k0
     $     ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,methode
     $     ,namefile,na,epsilon,polarisa,nmat,file_id ,group_iddip
     $     ,infostr ,nstop)

      use HDF5

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,nx1,ny1,nz1,na,nstop,ierror,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,x,y,z,eps0
     $     ,aretecube,xc,yc,zc
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp

      character*2 methode
      character*3 trope
      character(64) namefile
      character(64) infostr

      character(LEN=100) :: datasetname
      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error

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

c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  
c     read the input file
      open(15,file=namefile,status='old',iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary object: name of file does not exist'
         nstop=1
         return
      endif
      read(15,*) nx,ny,nz
      read(15,*) aretecube
      aretecube=aretecube*1.d-9

      ndipole=nx*ny*nz
      if (ndipole.gt.nmax) then
         write(*,*) 'Size of the parameter too small',ndipole
     $        ,'should be smaller that',nmax
         write(*,*) 'Increase nxm nym and nzm'
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif
      ii=0
      do i=1,nz
         do j=1,ny
            do k=1,nx
               ii=ii+1
               read(15,*) xs(ii),ys(ii),zs(ii)
               xs(ii)=xs(ii)*1.d-9
               ys(ii)=ys(ii)*1.d-9
               zs(ii)=zs(ii)*1.d-9
            enddo
         enddo
      enddo
      
      write(*,*)'You choose an arbitrary object that you define'
      write(*,*)'The code assumes that you define correctly' 
      write(*,*)'a cubic lattice for your object'
      


      ndipole=0
      do i=1,nz
         do j=1,ny
            do k=1,nx
               
               ndipole=ndipole+1
               nbsphere=nbsphere+1
               Tabdip(ndipole)=nbsphere

               if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) zs(i)
               if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) ys(j)
               if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) xs(k)

               
               if (trope.eq.'iso') then
                  read(15,*) eps
                  call poladiff(aretecube,eps,eps0,k0,dddis,methode
     $                 ,ctmp)
                  polarisa(nbsphere,1,1)=ctmp
                  polarisa(nbsphere,2,2)=ctmp
                  polarisa(nbsphere,3,3)=ctmp
                  epsilon(nbsphere,1,1)=eps
                  epsilon(nbsphere,2,2)=eps
                  epsilon(nbsphere,3,3)=eps
               else
                  do ii=1,3
                     do jj=1,3
                        read(15,*) epsani(ii,jj)
                     enddo
                  enddo
                  call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                 ,methode,inv,polaeps)
                  do ii=1,3
                     do jj=1,3
                        polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                        epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      if (ndipole.gt.nmax) then
         write(*,*) 'Size of the parameter too small',ndipole
     $        ,'should be smaller that',nmax
         write(*,*) 'Increase nxm nym and nzm'
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
      close(15)
      close(20)
      close(21)
      close(22)

      end
