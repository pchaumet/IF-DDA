c     proramme qui calcul la derivee
      subroutine derivative(nordre,nderiv,f,x,h,xderiv,deriv)
      implicit none
      integer nordre,nderiv
      double precision h,heps,x(0:4),xderiv
      double complex f(0:4),deriv

      heps=h/10.d0
c     derivee premiere
      if (nderiv.eq.1) then
         if (nordre.eq.3) then

            if (dabs(xderiv-x(0)).le.heps) then
               deriv=(-3.d0*f(0)+4.d0*f(1)-f(2))/2.d0/h
            elseif (dabs(xderiv-x(1)).le.heps) then
               deriv=(f(2)-f(0))/2.d0/h            
            elseif (dabs(xderiv-x(2)).le.heps) then
               deriv=(f(0)-4.d0*f(1)+3.d0*f(2))/2.d0/h
            endif

         elseif (nordre.eq.5) then
             if (dabs(xderiv-x(0)).le.heps) then
                deriv=(-50.d0*f(0)+96.d0*f(1)-72.d0*f(2)+32.d0*f(3)-6.d0
     $               *f(4))/24.d0/h
            elseif (dabs(xderiv-x(1)).le.heps) then
               deriv=(-6.d0*f(0)-20.d0*f(1)+36.d0*f(2)-12.d0*f(3)+2.d0
     $              *f(4))/24.d0/h            
            elseif (dabs(xderiv-x(2)).le.heps) then
               deriv=(2.d0*f(0)-16.d0*f(1)+16.d0*f(3)-2.d0*f(4))/24.d0/h
            elseif (dabs(xderiv-x(3)).le.heps) then
               deriv=(-2.d0*f(0)+12.d0*f(1)-36.d0*f(2)+20.d0*f(3)+6.d0
     $              *f(4))/24.d0/h     
            elseif (dabs(xderiv-x(4)).le.heps) then
               deriv=(6.d0*f(0)-32.d0*f(1)+72.d0*f(2)-96.d0*f(3)+50.d0
     $              *f(4))/24.d0/h
            endif
         else
            write(*,*) 'nordre not defined',nordre
         endif
      else
         write(*,*) 'nderiv not defined',nderiv
      endif
      return
      end
