      SUBROUTINE FACTOR235(N,IP,test)
      implicit none
      integer N,n2,test,IP(3)
C
c     si test est egal 1 non multiple de 2 ou 3 ou 5
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      test=0
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND. MOD(N,5) .NE. 0)
     $     THEN
         test=1
         RETURN
      ENDIF

   10 IF (N2 .LE. 1) RETURN
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
         test=1
      ENDIF
      RETURN
      END
