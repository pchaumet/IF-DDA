      subroutine calculatedate(valuesf,valuesi,tf,ti,message)
      implicit none
      integer valuesf(8),valuesi(8)
      double precision tf,ti,temps,tempsj
      character(64) message
      integer nlong
      nlong=len(trim(message))

      write(*,*) 'CPU time ',message(1:nlong),' in second  : ',tf-ti
      
      temps=(valuesf(5)-valuesi(5))*3600.d0+(valuesf(6)-valuesi(6))*60
     $     .d0+(valuesf(7)-valuesi(7))+(valuesf(8)-valuesi(8))/1000.d0
      if (valuesf(3).gt.valuesi(3)) then
         tempsj=(valuesf(3)-valuesi(3))*3600*24.d0
         write(*,*) 'Real time ',message(1:nlong),' in second : ',temps
     $        +tempsj
      elseif (valuesf(3).lt.valuesi(3)) then
         if (valuesi(3).eq.1) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.2) tempsj=dble(valuesf(3)+28-valuesi(3))
         if (valuesi(3).eq.3) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.4) tempsj=dble(valuesf(3)+30-valuesi(3))
         if (valuesi(3).eq.5) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.6) tempsj=dble(valuesf(3)+30-valuesi(3))
         if (valuesi(3).eq.7) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.8) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.9) tempsj=dble(valuesf(3)+30-valuesi(3))
         if (valuesi(3).eq.10) tempsj=dble(valuesf(3)+31-valuesi(3))
         if (valuesi(3).eq.11) tempsj=dble(valuesf(3)+30-valuesi(3))
         if (valuesi(3).eq.12) tempsj=dble(valuesf(3)+31-valuesi(3))
         tempsj=tempsj*3600*24.d0
         write(*,*) 'Real time ',message(1:nlong),' in second : ',temps
     $        +tempsj
      else
         write(*,*) 'Real time ',message(1:nlong),' in second : ',temps
      endif
      
      end
