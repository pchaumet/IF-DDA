      !ecrit une sous matrice dont les dimensions sont definies par
      !dim(1)=nx dim(2)=nxm 
      subroutine hdf5write1d(file_id, datasetname, data0, dim0)
      use HDF5
      implicit none
      integer(hid_t) :: file_id
      integer, dimension(4) :: dim0
      double precision, dimension(dim0(2)) :: data0
      
      character(LEN=100) :: datasetname
      integer(hid_t) :: dataset_id
      integer(hid_t) :: dataspace_id
      integer :: rank           !nb axes
      integer(hsize_t), dimension(1) :: dims
      
      integer(hid_t) :: memspace_id
      integer(hsize_t), dimension(1) :: memcount, memoffset, memstride,
     $     memblock
      
      integer(hid_t) :: group_id
      integer :: error
      
      integer i, j
      
      dims = dim0(1)
      rank = 1
      
                                !cree le dataspace du dataset
      call h5screate_simple_f(rank, dims, dataspace_id, error)
c      write(*,*) 'error',error
                                !cree le dataset avec les info du dataspace
      call h5dcreate_f(file_id, datasetname,  H5T_IEEE_F64BE,
     $     dataspace_id, dataset_id, error)
c      write(*,*) 'error1',error
      dims = dim0(2)
                                !cree espace memoire avec la sous matrice selon l hyperslab
      call h5screate_simple_f(rank, dims, memspace_id,error) 
c      write(*,*) 'error2',error
      memoffset(1) = 0
      memcount(1) = dim0(1)
      memstride(1) = 1
      memblock(1) = 1
                                !selectionne la sous matrice dans data0
      call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F,
     $     memoffset, memcount, error, memstride, memblock)
c     write(*,*) 'error3',error
                                !ecrit le dataset
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data0, dims, error,
     $     memspace_id, dataspace_id)
c      write(*,*) 'error4',error
                                !ferme le dataset
      call h5dclose_f(dataset_id, error)
c      write(*,*) 'error5',error
                                !ferme le dataspace du dataset
      call h5sclose_f(dataspace_id, error)
c      write(*,*) 'error',error
      end subroutine hdf5write1d
      
