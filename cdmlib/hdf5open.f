      ! ouvre un fichier hdf5 existant
      ! ne pas oublier de fermer avec hdf5close
      subroutine hdf5open(filename, file_id)
#ifdef USE_HDF5

         use HDF5
         implicit none
         character(LEN=100) :: filename
         integer(hid_t) :: file_id
         integer :: error 

         call h5open_f(error)
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
#endif
      end subroutine hdf5open
