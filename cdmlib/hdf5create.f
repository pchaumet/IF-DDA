      ! cree un fichier hdf5
      ! ne pas oublier de fermer avec hdf5close
      subroutine hdf5create(filename, file_id)
         use HDF5
         implicit none
         character(LEN=100) :: filename
         integer(hid_t) :: file_id
         integer :: error 

         call h5open_f(error)
         ! H5F_ACC_TRUNC_F cree si n existe pas sinon ecrase
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      end subroutine hdf5create
