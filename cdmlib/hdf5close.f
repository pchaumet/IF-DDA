      ! ferme un fichier hdf5 et l interface fortran
      subroutine hdf5close(file_id)
#ifdef USE_HDF5
         use HDF5
         implicit none
         integer(hid_t) :: file_id
         integer :: error

         call h5fclose_f(file_id, error)
         call h5close_f(error)
#endif
      end subroutine hdf5close
