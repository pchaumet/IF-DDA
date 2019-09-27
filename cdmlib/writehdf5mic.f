      subroutine writehdf5mic(Ex,Ey,Ez,nfft2d,imaxk0,Ediff,ntest ,name
     $     ,group_idmic)
      use HDF5
      
      implicit none

      integer nfft2d,nfft2d2,imaxk0,ntest,i,j,k,nimaxk0
      double complex Ex(nfft2d*nfft2d),Ey(nfft2d*nfft2d),Ez(nfft2d
     $     *nfft2d),Ediff(nfft2d*nfft2d,3)

      character(40) :: name
      character(LEN=100) :: datasetname
      integer(hid_t) :: group_idmic
      integer :: dim(4)

      integer ii,jj,indice
     
      nimaxk0=2*imaxk0+1
      
      if (ntest.eq.0) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nfft2d*nfft2d
            Ediff(i,1)= dsqrt(dreal(Ex(i)*dconjg(Ex(i))+Ey(i)
     $           *dconjg(Ey(i))+Ez(i)*dconjg(Ez(i))))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         
         dim(1)=nfft2d*nfft2d
         dim(2)=nfft2d*nfft2d
         datasetname=trim(name)//" field modulus"
         call hdf5write1d(group_idmic,datasetname,dreal(Ediff(:,1)),dim)
         datasetname=trim(name)//" field x component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ex),dim)
         datasetname=trim(name)//" field x component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ex),dim)
         datasetname=trim(name)//" field y component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ey),dim)
         datasetname=trim(name)//" field y component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ey),dim)
         datasetname=trim(name)//" field z component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ez),dim)
         datasetname=trim(name)//" field z component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ez),dim)


      else
         nfft2d2=nfft2d/2
         dim(1)=nimaxk0*nimaxk0
         dim(2)=nfft2d*nfft2d
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nfft2d*nfft2d
            Ediff(i,1)=0.d0
            Ediff(i,2)=0.d0
            Ediff(i,3)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,indice,k)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         do j=-imaxk0,imaxk0
            do i=-imaxk0,imaxk0
               ii=i+imaxk0+1
               jj=j+imaxk0+1
               k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
               indice=ii+nimaxk0*(jj-1)                 
               Ediff(indice,1)=dsqrt(dreal(Ex(k)*dconjg(Ex(k))+Ey(k)
     $              *dconjg(Ey(k))+Ez(k)*dconjg(Ez(k))))
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL


         datasetname=trim(name)//" field modulus"        
         call hdf5write1d(group_idmic,datasetname,dreal(Ediff(:,1)),dim)

         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,indice,k)   
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         do j=-imaxk0,imaxk0
            do i=-imaxk0,imaxk0
               ii=i+imaxk0+1
               jj=j+imaxk0+1
               k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
               indice=ii+nimaxk0*(jj-1)
               Ediff(indice,1)=Ex(k)
               Ediff(indice,2)=Ey(k)
               Ediff(indice,3)=Ez(k)
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         datasetname=trim(name)//" field x component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ediff(:,1)),dim)
         datasetname=trim(name)//" field x component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ediff(:,1)),dim)
         datasetname=trim(name)//" field y component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ediff(:,2)),dim)
         datasetname=trim(name)//" field y component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ediff(:,2)),dim)
         datasetname=trim(name)//" field z component real part"
         call hdf5write1d(group_idmic,datasetname,dreal(Ediff(:,3)),dim)
         datasetname=trim(name)//" field z component imaginary part"
         call hdf5write1d(group_idmic,datasetname,dimag(Ediff(:,3)),dim)

         
      endif



      end
