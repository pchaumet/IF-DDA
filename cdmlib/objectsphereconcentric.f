      subroutine objetsphereconcentric(trope,epss,epsanis,numbersphere
     $     ,numberspheremax,xg,yg,zg,rayons,eps0,xs,ys,zs,k0,aretecube
     $     ,tabdip,tabnbs,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,methode
     $     ,na,epsilon,polarisa,nmat,file_id,group_iddip,infostr,nstop)
      use HDF5
      implicit none
      integer nmax,tabdip(nmax),tabnbs(nmax),nbsphere,ndipole,nx,ny,nz
     $     ,na,ii,jj,i,j,k,l,test,IP(3),nnnr,dddis,inv,is,numbersphere
     $     ,numberspheremax,nstop,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg,yg,zg,x,y,z
     $     ,eps0,aretecube ,rayons(numberspheremax),ray,centre
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),epss(numberspheremax)
     $     ,epsanis(3,3,numberspheremax),ctmp
      character*2 methode
      character*3 trope
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
         tabnbs(i)=0
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

      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      if (rayons(1).eq.0.d0) then
         infostr='spheres concentrics: 1st rayon=0'
         nstop=1
         return
      endif
      do is=1,numbersphere
         rayons(is)=rayons(is)*1.d-9        
         if (is.ge.2.and.rayons(is)-rayons(is-1).le.0.d0) then
            infostr='spheres concentrics: rayons non croissants'
            nstop=1
            return
         endif
      enddo
      aretecube=2.d0*rayons(numbersphere)/dble(nnnr)

      nx=nnnr
      ny=nnnr
      nz=nnnr
      centre=aretecube*dble(nnnr)/2.d0
      write(*,*) 'eps concentric',epss
      write(*,*) 'Box including the  spheres concentriques',nx,ny,nz
     $     ,aretecube,na,numbersphere    
      if (na.eq.-1) then
         do i=1,nz
            do j=1,ny
               do k=1,nx
                  x=dble(k-1)*aretecube+aretecube/2.d0-centre
                  y=dble(j-1)*aretecube+aretecube/2.d0-centre
                  z=dble(i-1)*aretecube+aretecube/2.d0-centre
                  ray=dsqrt(x*x+y*y+z*z)

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x
                  
                  ndipole=ndipole+1
                  is=0
                  if (ray.le.rayons(1)) then
                     is=1
                  else
                     do l=2,numbersphere
                        if (ray.le.rayons(l).and.ray.ge.rayons(l-1))is=l
                     enddo
                  endif
c                  write(*,*) 'is sortant',is,ray,i,j,k
                  if (is.ne.0) then
                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere
                     Tabnbs(nbsphere)=is
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     if (trope.eq.'iso') then
                        eps=epss(is)
c                        write(*,*) 'eps sortant',is,eps,ray,i,j,k
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        do ii=1,3
                           do jj=1,3                              
                              epsani(ii,jj)=epsanis(ii,jj,is)
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                           enddo
                        enddo
                        call polaepstens(aretecube,epsani,eps0,k0
     $                       ,dddis,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                  endif
               enddo             
            enddo
         enddo
      elseif (na.eq.0) then
         do i=1,nz
            do j=1,ny
               do k=1,nx
                  x=dble(k-1)*aretecube+aretecube/2.d0-centre
                  y=dble(j-1)*aretecube+aretecube/2.d0-centre
                  z=dble(i-1)*aretecube+aretecube/2.d0-centre
                  ray=dsqrt(x*x+y*y+z*z)

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  is=0
                  if (ray.le.rayons(1)) then
                     is=1
                  else
                     do l=2,numbersphere
                        if (ray.le.rayons(l).and.ray.ge.rayons(l-1))is=l
                     enddo
                  endif                  
                  if (is.ne.0) then
                     Tabnbs(nbsphere)=is
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     if (trope.eq.'iso') then
                        eps=epss(is)
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        do ii=1,3
                           do jj=1,3
                              epsani(ii,jj)=epsanis(ii,jj,is)
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                           enddo
                        enddo
                        call polaepstens(aretecube,epsani,eps0,k0
     $                       ,dddis,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                  else                     
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
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
      close(15)
      close(20)
      close(21)
      close(22)

      end
