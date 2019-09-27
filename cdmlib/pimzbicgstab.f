      SUBROUTINE PIMZBICGSTAB(X,Xi,XR,B,lda,nlar,ndim,nou,WRK,ALPHA
     $     ,OMEGA,RHO,NORM,TOLE,ITNO,MAXIT,STATUS,STEPERR)
      IMPLICIT NONE

      integer lda,ndim,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),X(lda),xr(lda),xi(lda)
      
*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,KAPPA,OMEGA,RHO,RHO0,XXI,ctmp1,ctmp2
      DOUBLE PRECISION NORM,TOLE
      INTEGER NOU,II,ITNO,MAXIT,STATUS,STEPERR
      
      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 100

*     Set indices for mapping local vectors into wrk
c     IR = 1
c     IRTILDE =2
c     IP = 3
c     IQ = 4
c     IS = 5
c     IT = 6
c     IV = 7
c     IW = 8
c     IZ = 9
c     IXOLD = 10

*     Loop
      STATUS = 0
      STEPERR = -1

*     1.    r=b-Ax
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)
      do ii=1,ndim
         wrk(ii,1)=b(ii)
         x(ii)=xi(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=1
      return

 10   NORM=0.d0
*     2. rtilde=r, 3. p=v=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)      
      do ii=1,ndim
         wrk(ii,8)=xr(ii)
         WRK(ii,1)=WRK(ii,1)-wrk(ii,8)
         NORM=NORM+dreal(WRK(ii,1)*dconjg(WRK(ii,1)))
         WRK(ii,2)=WRK(ii,1) 
         WRK(ii,3)=0.d0
         WRK(ii,7)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
*     4. rho=alpha=omega=1
      RHO = (1.d0,0.d0)
      ALPHA = (1.d0,0.d0)
      OMEGA = (1.d0,0.d0)


      ITNO=0
 100  ITNO=ITNO+1
*     5. rho=dot(rtilde,r)
      RHO0 = RHO
      RHO=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RHO)         
      do ii=1,ndim
c         RHO=RHO+dconjg(WRK(ii,2))*WRK(ii,1)
         RHO=RHO+dconjg(WRK(ii,1))*WRK(ii,1) 
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
*     6. beta=rho*alpha/(rho0*omega)
      KAPPA = RHO0*OMEGA
      IF (cdabs(KAPPA).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 6
         return
      END IF

      BETA = RHO*ALPHA/KAPPA

*     7. p=r+beta*(p-omega*v), 8. v=Q1AQ2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
c         WRK(ii,3)=-OMEGA*WRK(ii,7)
         WRK(ii,3)=WRK(ii,3)-OMEGA*WRK(ii,7)
         WRK(ii,8)=WRK(ii,3)
         WRK(ii,3)=WRK(ii,1)
         WRK(ii,3)=WRK(ii,3)+BETA*WRK(ii,8)
         xi(ii)=WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return

 20   XXI=0.d0
*     9. xi=dot(rtilde,v)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:XXI)              
      do ii=1,ndim
         WRK(ii,7)=xr(ii)    
         XXI=XXI+dconjg(WRK(ii,2))*WRK(ii,7)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
*     10. alpha=rho/xi
      IF (cdabs(XXI).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 10
         GO TO 9999
      END IF

      ALPHA = RHO/XXI
      KAPPA = 0.d0
*     11. s=r-alpha*v, 12. if ||s||<breaktol then soft-breakdown has occurred
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:KAPPA)
      do ii=1,ndim
         WRK(ii,5)=WRK(ii,1)    
         WRK(ii,5)=WRK(ii,5)-ALPHA*WRK(ii,7)
         KAPPA=KAPPA+dconjg(WRK(ii,5))*WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      KAPPA=cdsqrt(KAPPA)

      IF (cdABS(KAPPA).LT.1.d-15) THEN
         STATUS = -2
         STEPERR = 12
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)         
         do ii=1,ndim
            x(ii)=WRK(ii,3)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         GO TO 9999
      END IF

*     13. t=Q1AQ2s
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) 
      do ii=1,ndim
         xi(ii)=WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL       
      nou=3
      return

 30   ctmp1=0.d0
      ctmp2=0.d0
*     14. omega=dot(t,s)/dot(t,t)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1,ctmp2) 
      do ii=1,ndim
         WRK(ii,6)=xr(ii)  
         ctmp1=ctmp1+dconjg(WRK(ii,6))*WRK(ii,6)
         ctmp2=ctmp2+dconjg(WRK(ii,6))*WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      IF (cdabs(ctmp1).EQ.0.d0) THEN
         STATUS = -3
         STEPERR = 14
         GO TO 9999
      END IF

      OMEGA = ctmp2/ctmp1

*     15. x=x+alpha*p+omega*s, 16. r=s-omega*t
      ctmp1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1)       
      do ii=1,ndim
         x(ii)=x(ii)+ALPHA*WRK(ii,3)+OMEGA*WRK(ii,5)
         WRK(ii,1)=WRK(ii,5)
         WRK(ii,1)=WRK(ii,1)-OMEGA*WRK(ii,6)
         ctmp1=ctmp1+dconjg(WRK(ii,1))*WRK(ii,1)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      ctmp1=ctmp1/norm
*     17. check stopping criterion
      IF (ITNO.GT.MAXIT) THEN
         STATUS = -1
         ITNO = MAXIT
         goto 9999
      END IF
      if (cdabs(ctmp1).le.TOLE/10.d0) then
         STATUS=1
         nou=4
         return
      END IF

     

      goto 100


 9999 CONTINUE
      RETURN
      END
