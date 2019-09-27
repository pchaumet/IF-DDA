      subroutine comparaisonreel(a,aref,test)
      implicit none
      integer test
      double precision a,aref,rap
      test=1
      if (aref.eq.0.d0) then
         if (a.eq.0.d0) test=0
         return
      else
         rap=dabs((a-aref)/aref)
         if (rap.le.1.d-12) test=0         
         return
      endif
      end
c**************************************************
      subroutine comparaisoncomplexe(a,aref,test)
      implicit none
      integer test
      double precision rap1,rap2,arefr,arefi
      double complex a,aref
      test=1
     
      arefr=dreal(aref)
      arefi=dimag(aref)
      if (arefr.eq.0.d0.and.dreal(a).eq.0.d0)  then
         rap1=0.d0
      else
         rap1=dabs((dreal(a)-arefr)/arefr)
      endif

      if (arefi.eq.0.d0.and.dimag(a).eq.0.d0)  then
         rap2=0.d0
      else
         rap2=dabs((dimag(a)-arefi)/arefi)
      endif
      
      if (rap1.le.1.d-12.and.rap2.le.1.d-12) test=0

      return

      end
