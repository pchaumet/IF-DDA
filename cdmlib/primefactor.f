c      integer n,ncomp,test
c      read(*,*) n
c      call  PRIMEFACTOR(N,1,test)
c
c      END
c     ***********************************************************
      SUBROUTINE PRIMEFACTOR(N,ncomp,test)
      implicit none
      integer N,n2,test,IP(10),ncomp
C
c     si test est egal 1 non multiple de 2 ou 3 ou 5
      IP=0
      N2=N
      test=0
      
 10   IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      ELSE IF (n2.gt.1) then
 20      IF (MOD(N2,7) .EQ. 0) THEN
            IP(4)=IP(4)+1
            N2=N2/7
            GO TO 20
         ELSEIF (MOD(N2,11) .EQ. 0) THEN
            IP(5)=IP(5)+1
            N2=N2/11
            GO TO 20
         ELSEIF (MOD(N2,13) .EQ. 0) THEN
            IP(6)=IP(6)+1
            N2=N2/13
            GO TO 20
         ELSEIF (MOD(N2,17) .EQ. 0) THEN
            IP(7)=IP(7)+1
            N2=N2/17
            GO TO 20
         ELSEIF (MOD(N2,19) .EQ. 0) THEN
            IP(8)=IP(8)+1
            N2=N2/19
            GO TO 20
         ELSEIF (MOD(N2,23) .EQ. 0) THEN
            IP(9)=IP(9)+1
            N2=N2/23
            GO TO 20        
         ELSE
            IP(10)=N2
         ENDIF
         write(*,*) 'High prime factor: FFT can be slowed down'
         if (ncomp.eq.1) then
            write(*,*) 'for the discretization along the x axis'
         elseif (ncomp.eq.2) then
            write(*,*) 'for the discretization along the y axis'
         else
            write(*,*) 'for the discretization along the z axis'
         endif
         write(*,FMT="('2^',I1'*3^',I1,'*5^',I1,'*7^',I1,'*11^',
     $   I1,'*13^',I1,'*17^',I1,'*19^',I1,'*23^',I1,'*N=',I3)")
     $  IP(1),IP(2),IP(3) ,IP(4) ,IP(5) ,IP(6),IP(7),IP(8) ,IP(9),IP(10)
         
      ENDIF

      RETURN
      END
