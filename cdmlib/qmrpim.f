      SUBROUTINE PIMZQMR(XS,XI,XR,B,WRK,NORM,LDA,NDIM,NLAR,LAMBDA,KAPPA
     $     ,THETA,GAMMA,KSI,RHO,EPSILON,MU,TAU,DOTS1,DOTS2,DOTS3,DOTS4
     $     ,NOU,NT,ITNO,MAXIT ,TOLE ,TOL ,STATUS ,STEPERR)

      IMPLICIT NONE

*     .. Array Arguments ..
      DOUBLE COMPLEX B(LDA),WRK(LDA,NLAR),XI(LDA),XR(LDA),XS(LDA)
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX ABSGAMMA2,ABSTAU02,DEN,GAMMA,EPSILON,EPSILON0,
     +     GAMMA0,KAPPA,KAPPA0,KSI,KSI0,LAMBDA,LAMBDA0,MU,MU0,
     +     RHO,RHO0,TAU,TAU0,THETA,TMP1
      DOUBLE PRECISION NORM, TOL, TOLE
      INTEGER NDIM, NLAR, NOU, I, ITNO, LDA ,MAXIT, STATUS,STEPERR,NT
*     ..
*     .. Local Arrays ..
      DOUBLE COMPLEX DOTS1,DOTS2,DOTS3,DOTS4
*     ..

      if (nou.eq.0) goto 10
      if (nou.eq.1) goto 20
      if (nou.eq.2) goto 30
      if (nou.eq.3) goto 40
      if (nou.eq.4) goto 50
      if (nou.eq.5) goto 100

*     1. lambda=1, kappa=-1, theta=-1
 10   LAMBDA = (1.d0,0.d0)
      KAPPA = -(1.d0,0.d0)
      THETA = -(1.d0,0.d0)
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)           
      DO I=1,NDIM
         NORM=NORM+dreal(B(I)*DCONJG(B(I)))
         xs(i)=xi(i)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      
      NORM=dsqrt(NORM)

*     Loop
      STATUS = 0
      STEPERR = -1
      ITNO=0

*     2. wtilde=vtilde=r=b-Ax
*     r=b-Ax

c     A*x=wrk(3)
      nou=1
      nt=1
 
c     calcul de A*xi

      Return

 20   DOTS1 = 0.d0
      DOTS2 = 0.d0
      DOTS3 = 0.d0
c     3. p=q=d=s=0,
*     4. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
*     rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde, mu=0
*     Compute A^{T}wtilde
c     CALL TMATVEC(WRK(IWTILDE),WRK(IATWTILDE),IPAR)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DOTS1,DOTS2,DOTS3)      
      do i=1,ndim
         wrk(i,3)=xr(i)
         wrk(i,1)=B(i)-wrk(i,3)
         wrk(i,7)=wrk(i,1)
         wrk(i,8)=wrk(i,1)
         wrk(i,2)=0.d0
         wrk(i,4)=0.d0
         wrk(i,5)=0.d0
         wrk(i,6)=0.d0
         DOTS1 = DOTS1+wrk(i,7)*dconjg(wrk(i,7))
         DOTS2 = DOTS2+wrk(i,8)*dconjg(wrk(i,8))
         DOTS3 = DOTS3+wrk(i,7)*wrk(i,8)  
         XI(I)= WRK(I,8)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nt=2
      nou=2
      return


 30   DOTS4 = 0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DOTS4)          
      do i=1,ndim
         wrk(i,9)=xr(i)
         DOTS4 = DOTS4+wrk(i,7)*wrk(i,9)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL       

*     Accumulate simultaneously partial inner-products
c     CALL PDZSUM(4,DOTS)

      GAMMA = CDSQRT(DOTS1)
      KSI = CDSQRT(DOTS2)
      RHO = DOTS3
      EPSILON = DOTS4
      MU = 0.d0
*     5. tau=epsilon/rho
      IF (cdabs(RHO).EQ.0.d0) THEN
         ITNO = 0
         STATUS = -3
         STEPERR = 5
         GO TO 200
      END IF

      TAU = EPSILON/RHO
 100  ITNO=ITNO+1
      


      IF (cdabs(GAMMA).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 6
         GO TO 200
      END IF
      IF (cdabs(KSI).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 7
         GO TO 200
      END IF

*     6. p=1/gamma*vtilde-mu*p
*     7    . q=1/ksi*A^{T}wtilde-(gamma*mu)/ksi*q
*     8. vtilde=Ap-tau/gamma*vtilde
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)      
      do i=1,ndim
         wrk(i,2)=wrk(i,7)/gamma-mu*wrk(i,2)
         wrk(i,4)=(wrk(i,9)-gamma*mu*wrk(i,4))/ksi 
         XI(I)= WRK(I,2)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL
      
      nt=1
      nou=3
      return

 40   DOTS1 = 0.d0
      DOTS2 = 0.d0
      DOTS3 = 0.d0
*     9. wtilde=q-tau/ksi*wtilde
*     11. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
*     rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde
*     Compute A^{T}wtilde
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DOTS1,DOTS2,DOTS3)         
      do i=1,ndim
         wrk(i,3)=xr(i)
         wrk(i,7)=wrk(i,3)-tau/gamma*wrk(i,7)
         wrk(i,8)=wrk(i,4)-tau/ksi*wrk(i,8)     
         DOTS1 = DOTS1+wrk(i,7)*dconjg(wrk(i,7))
         DOTS2 = DOTS2+wrk(i,8)*dconjg(wrk(i,8))
         DOTS3 = DOTS3+wrk(i,7)*wrk(i,8)     
         XI(I)= WRK(I,8)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL
      
      nt=2
      nou=4
      return
 50   DOTS4 = 0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DOTS4)
      do i=1,ndim
         wrk(i,9)=xr(i)
         DOTS4 = DOTS4+wrk(i,7)*wrk(i,9)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
*     Accumulate simultaneously partial inner-products
c     CALL PDZSUM(4,DOTS)

      GAMMA0 = GAMMA
      GAMMA = CDSQRT(DOTS1)
      KSI0 = KSI
      KSI = CDSQRT(DOTS2)
      RHO0 = RHO
      RHO = DOTS3
      EPSILON0 = EPSILON
      EPSILON = DOTS4
*     12. mu=(gamma0*ksi0*rho)/(gamma*tau*rho0)

      DEN = GAMMA*TAU*RHO0
      IF (cdabs(DEN).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 12
         GO TO 200

      END IF
      MU0 = MU
      MU = (GAMMA0*KSI0*RHO)/DEN

*     13. tau=epsilon/rho-gamma*mu
      IF (cdabs(RHO).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 13
         GO TO 200

      END IF
      TAU0 = TAU
      TAU = EPSILON/RHO - GAMMA*MU
*     14. theta=(|tau0|^2*(1-lambda))/(lambda*|tau|^2+|gamma|^2)

      ABSTAU02 = CDABS(TAU0)**2.D0
      ABSGAMMA2 = CDABS(GAMMA)**2.D0
      DEN = LAMBDA*ABSTAU02 + ABSGAMMA2
      IF (cdabs(DEN).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 14
         GO TO 200

      END IF
      THETA = (ABSTAU02* ((1.d0,0.d0)-LAMBDA))/DEN

*     15. kappa=(-gamma0*CONJG(tau0)*kappa0)/(gamma0*|tau|^2+|gamma|^2)
      KAPPA0 = KAPPA
      KAPPA = - (GAMMA0*DCONJG(TAU0)*KAPPA0)/DEN

*     16. lambda=(lambda0*|tau0|^2)/(gamma0*|tau|^2+|gamma|^2)
      LAMBDA0 = LAMBDA
      LAMBDA = LAMBDA0*ABSTAU02/DEN

*     17. d=theta*d+kappa*p
*     18. s=theta*s+kappa*A*p
*     19. x=x+d
*     20. r=r-s
c     critere d'arret
      TMP1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:TMP1)
      do i=1,ndim
         wrk(i,5)=theta*wrk(i,5)+kappa*wrk(i,2)    
         wrk(i,6)=theta*wrk(i,6)+kappa*wrk(i,3) 
         xs(i)=xs(i)+wrk(i,5) 
         wrk(i,1)=wrk(i,1)-wrk(i,6)   
         TMP1=TMP1+WRK(I,1)*DCONJG(WRK(I,1))
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      
      TOLE=dsqrt(cdabs(TMP1))/NORM
      if (mod(ITNO,50).EQ.0) write(*,*) 'RESIDUE',TOLE,'iteration',ITNO
      IF (TOLE.LE.TOL) then
         NOU=5
         STATUS=1
         RETURN         
      ENDIF

      IF (ITNO.GT.MAXIT) THEN
         STATUS = -1
         ITNO = MAXIT
      END IF
      GOTO 100

 200  RETURN

      END

