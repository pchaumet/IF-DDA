      subroutine objetparanxnynz(trope,eps,epsani,eps0,xs,ys,zs,k0
     $     ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym
     $     ,nzm,nxmp,nymp,nzmp,methode,epsilon ,polarisa,sidex,sidey
     $     ,sidez,xg,yg,zg,na ,nmat,file_id,group_iddip,infostr,nstop)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,na,nstop,nxm,nym,nzm,nmat,nxmp
     $     ,nymp,nzmp
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg,yg,zg,x,y,z,eps0
     $     ,aretecube,sidex,sidey,sidez,side,pi,mat(3,3) ,xr,yr,zr,x1,x2
     $     ,y1,y2,z1,z2
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp

      character*2 methode
      character*3 trope
      character(64) infostr

      character(LEN=100) :: datasetname

#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif
      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error

      write(*,*) 'object2'
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
      if (nmat.eq.0) then
c     mesh 
         open(20,file='x.mat')
         open(21,file='y.mat')
         open(22,file='z.mat')  
c     discretization of the object under study
         open(10,file='xc.mat')
         open(11,file='yc.mat')
         open(12,file='zc.mat')  
      endif

      nx=nxm-2*nxmp
      ny=nym-2*nymp
      nz=nzm-2*nzmp
      aretecube=aretecube*1.d-9
      sidex=aretecube*dble(nx)
      sidey=aretecube*dble(ny)
      sidez=aretecube*dble(nz)
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      

      write(*,*) 'object4',nx,nxm,nxmp
      if (sidex.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidex=0!'
         return
      elseif (sidey.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidey=0!'
         return
      elseif (sidez.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidez=0!'
         return
      endif
      write(*,*) 'object'

      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif
      
      do i=1,nz
         do j=1,ny
            do k=1,nx

               x=aretecube*(dble(k)-0.5d0)-sidex/2.d0
               y=aretecube*(dble(j)-0.5d0)-sidey/2.d0
               z=aretecube*(dble(i)-0.5d0)-sidez/2.d0               

               if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(22,*) z+zg
               if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(21,*) y+yg
               if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(20,*) x+xg

               ndipole=ndipole+1                   
               nbsphere=nbsphere+1
               
               Tabdip(ndipole)=nbsphere
               xs(nbsphere)=x+xg
               ys(nbsphere)=y+yg
               zs(nbsphere)=z+zg                   
               
               if (trope.eq.'iso') then
                  call poladiff(aretecube,eps,eps0,k0,dddis
     $                 ,methode,ctmp)  
                  polarisa(nbsphere,1,1)=ctmp
                  polarisa(nbsphere,2,2)=ctmp
                  polarisa(nbsphere,3,3)=ctmp
                  epsilon(nbsphere,1,1)=eps
                  epsilon(nbsphere,2,2)=eps
                  epsilon(nbsphere,3,3)=eps
               else
                  call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                 ,methode,inv,polaeps)
                  do ii=1,3
                     do jj=1,3
                        epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                        polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo 

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

      write(*,*) 'Cuboid::nbsphere=',nbsphere,nx,ny,nz
      end
