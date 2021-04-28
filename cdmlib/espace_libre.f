c     propagateur d espace libre
      subroutine propa_espace_libre(x,y,z,x0,y0,z0,k0,
     *     propaesplibre)
      implicit none
      integer i,j
      double precision x,y,z,x0,y0,z0,k0,Id,Rab,Rtenseur,Rvect,pi,Rab2
      double complex propaesplibre,const1,const2
      dimension propaesplibre(3,3),Id(3,3),Rtenseur(3,3),Rvect(3)

      pi=dacos(-1.d0)
      Rab2=0.d0
      Rvect(1)=(x-x0)
      Rvect(2)=(y-y0)
      Rvect(3)=(z-z0)

      do i=1,3
         do j=1,3
            Id(i,j)=0.d0
            if (i.eq.j) Id(i,i)=1.d0
            Rtenseur(i,j)=Rvect(i)*Rvect(j)
         enddo
         Rab2=Rab2+Rvect(i)*Rvect(i)
      enddo
      Rab=dsqrt(Rab2)

      if (Rab.eq.0.d0) then
         do i=1,3
            do j=1,3
               propaesplibre(i,j)=(0.d0,0.d0) 
	   enddo
         enddo
      else
         const1=((1.d0,0.d0)/Rab-(0.d0,1.d0)*k0)/Rab2
         const2=k0*k0/Rab*(1.d0,0.d0)
         do i=1,3
            do j=1,3
               propaesplibre(i,j)=((3.d0*Rtenseur(i,j)/Rab2-Id(i,j))
     $              *const1+(Id(i,j)-Rtenseur(i,j)/Rab2)*const2)
     $              *cdexp((0.d0,1.d0)*k0*Rab)              
            enddo
         enddo
      endif
      end
