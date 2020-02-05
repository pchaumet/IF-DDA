      subroutine objetnspheres(trope,epss,epsanis,numbersphere
     $     ,numberspheremax,xg,yg,zg,rayons,eps0,xs,ys,zs ,k0,aretecube
     $     ,tabdip,tabnbs,nnnr,nmax,nbsphere,ndipole,nx ,ny,nz,methode
     $     ,na,epsilon,polarisa,nmat,file_id ,group_iddip,infostr,nstop)

#ifdef USE_HDF5
      use HDF5
#endif
      
      implicit none
      integer nmax,tabdip(nmax),tabnbs(nmax),nbsphere,ndipole,nx,ny,nz
     $     ,na,ii,jj,i,j,k,test,IP(3),nnnr,dddis,inv,is,numbersphere
     $     ,numberspheremax,nstop,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg(numberspheremax)
     $     ,yg(numberspheremax),zg(numberspheremax),x,y,z,eps0,aretecube
     $     ,rayons(numberspheremax),xmin,xmax,ymin,ymax,zmin,zmax
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),epss(numberspheremax)
     $     ,epsanis(3,3,numberspheremax),ctmp
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
      
      write(*,*) 'Begin of nspheres'
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

      xmax=-1.d300
      xmin=1.d300
      ymax=-1.d300
      ymin=1.d300
      zmax=-1.d300
      zmin=1.d300
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      
      
      do is=1,numbersphere
         rayons(is)=rayons(is)*1.d-9
         xg(is)=xg(is)*1.d-9
         yg(is)=yg(is)*1.d-9
         zg(is)=zg(is)*1.d-9
         xmax=max(xmax,xg(is)+rayons(is))
         xmin=min(xmin,xg(is)-rayons(is))
         ymax=max(ymax,yg(is)+rayons(is))
         ymin=min(ymin,yg(is)-rayons(is))
         zmax=max(zmax,zg(is)+rayons(is))
         zmin=min(zmin,zg(is)-rayons(is))
         if (rayons(is).eq.0.d0) then
            nstop=1
            infostr='object nspheres: radius=0!'
            return
         endif
      enddo
      if (numbersphere.gt.numberspheremax) then
         infostr='number of spheres too large: limited to 20'
         nstop=1
         return
      endif

      aretecube=max(zmax-zmin,ymax-ymin,xmax-xmin)/dble(nnnr)


      nx=idnint((xmax-xmin)/aretecube)
      ny=idnint((ymax-ymin)/aretecube)
      nz=idnint((zmax-zmin)/aretecube)

      write(*,*) 'Box including the N spheres',nx,ny,nz,'meshsize'
     $     ,aretecube
      write(*,*) 'Number of sphere',numbersphere
      write(*,*) 'X',xmin*1.d9,xmax*1.d9
      write(*,*) 'Y',ymin*1.d9,ymax*1.d9
      write(*,*) 'Z',zmin*1.d9,zmax*1.d9
      if (na.eq.-1) then
         write(*,*) 'na',na
         do i=1,nz
            do j=1,ny
               do k=1,nx
                  x=xmin+dble(k-1)*aretecube+aretecube/2.d0
                  y=ymin+dble(j-1)*aretecube+aretecube/2.d0
                  z=zmin+dble(i-1)*aretecube+aretecube/2.d0

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.eq.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x

                  ndipole=ndipole+1
                  test=0
                  do is=1,numbersphere
                     if ((x-xg(is))**2+(y-yg(is))**2+(z-zg(is))
     $                    **2.le.rayons(is)**2) then
                        test=test+1
                        nbsphere=nbsphere+1
                        Tabdip(ndipole)=nbsphere
                        Tabnbs(nbsphere)=is
                        xs(nbsphere)=x
                        ys(nbsphere)=y
                        zs(nbsphere)=z
                        if (trope.eq.'iso') then
                           eps=epss(is)
                           call poladiff(aretecube,eps,eps0,k0,dddis
     $                          ,methode,ctmp)  
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
     $                          ,dddis,methode,inv,polaeps)
                           do ii=1,3
                              do jj=1,3
                                 polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                              enddo
                           enddo
                        endif
                     endif
                  enddo
                  if (test.ge.2) then
                     infostr='sphere are not distincts'
                     nstop=1
                     return
                  endif
               enddo
            enddo
         enddo
      elseif (na.eq.0) then
         write(*,*) 'na',na
         do i=1,nz
            do j=1,ny
               do k=1,nx
                  x=xmin+dble(k-1)*aretecube+aretecube/2.d0
                  y=ymin+dble(j-1)*aretecube+aretecube/2.d0
                  z=zmin+dble(i-1)*aretecube+aretecube/2.d0
                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  test=0
                  do is=1,numbersphere
                     if ((x-xg(is))**2+(y-yg(is))**2+(z-zg(is))
     $                    **2.le.rayons(is)**2) then
                        test=test+1
c                        write(*,*) 'test0',test,nbsphere
                        Tabnbs(nbsphere)=is
                        xs(nbsphere)=x
                        ys(nbsphere)=y
                        zs(nbsphere)=z
                        if (trope.eq.'iso') then
                           eps=epss(is)
                           call poladiff(aretecube,eps,eps0,k0,dddis
     $                          ,methode,ctmp)  
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
     $                          ,dddis,methode,inv,polaeps)
                           do ii=1,3
                              do jj=1,3
                                 polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                              enddo
                           enddo
                        endif
                     endif
                  enddo
                  if (test.eq.0) then
c                     write(*,*) 'test1',test,nbsphere
                     xs(nbsphere)=x
                     ys(nbsphere)=y
                     zs(nbsphere)=z  
                     epsilon(nbsphere,1,1)=(1.d0,0.d0)
                     epsilon(nbsphere,2,2)=(1.d0,0.d0)
                     epsilon(nbsphere,3,3)=(1.d0,0.d0)
                  elseif (test.ge.2) then
                     infostr='sphere are not distincts'
                     nstop=1
                     return
                  endif
               enddo
            enddo
         enddo

      else
         infostr='na should be equal to -1 or 0'
         nstop=1
         return
      endif
      write(*,*) 'End of nspheres'
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
