      subroutine ZCG(XI,XR,B,NORM,WRK,ALPHA,RHO0,XXI,LDA,NDIM,NLAR,NOU
     $     ,NSTAT,NLOOP,MAXIT,TOLE,TOL)
      IMPLICIT NONE

      INTEGER I,NLOOP,MAXIT


      

C     array
      INTEGER LDA,NDIM,NLAR,NOU,NSTAT
      DOUBLE COMPLEX B(LDA),WRK(LDA,NLAR),XR(LDA),XI(LDA)

c     Local Scalars
      DOUBLE PRECISION TOL,NORM,TOLE
      DOUBLE COMPLEX ALPHA,BETA,DELTA,RHO,RHO0,XXI,TMP1

c     en entre nou=0,nloop=0

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40
c     NOU=4 si change tolerance en cours de route
      
*  1. r=(b-Ax) CALCUL DU RESIDU
      NOU=1  
*     Loop
      NSTAT=0
c     stock l'estimation initiale
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)      
      DO I=1,NDIM
         WRK(I,1)=XI(I)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL
     
C     CALCUL Ax AL'EXTERIEUR ET RENVOIE DANS XR
      RETURN
C     CALCUL b-Ax
 10   NORM=0.D0
      RHO0=0.d0
c     2. p=r
*     3. rdotr=dot(r,r)
*     4. w=Ap
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM,RHO0)     
      DO I=1,NDIM
C     RK-1
         WRK(I,2)=B(I)-XR(I)
C     INITIALISE XK
         NORM=NORM+dreal(B(I)*DCONJG(B(I)))
         WRK(I,3)=WRK(I,2)
         RHO0=RHO0+WRK(I,2)*WRK(I,2)
         XI(I)=WRK(I,3)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      
      NORM=dsqrt(NORM)
      NOU=2
    
      RETURN

 20   XXI=0.d0
*     5. xi=dot(p,w)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:XXI)   
      DO I=1,NDIM
         WRK(I,5)=XR(I)    
         XXI=XXI+WRK(I,3)*WRK(I,5)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      
      if (XXI.eq.0.d0) then
         write(*,*) 'problem with XXI',XXI
         stop
      endif

*     Loop
     
 50   NLOOP=NLOOP+1

      IF (NLOOP.GE.MAXIT) then
         NSTAT=-4
         NOU=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)  
         DO I=1,NDIM
            XI(I)=WRK(I,1)
         ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL         
         RETURN         
      ENDIF
*     6. alpha=rdotr/xi      
      ALPHA=RHO0/XXI
c      write(*,*) 'ALPHA',ALPHA
     

*     7. x=x+alpha*p
*     8. r=r-alpha*w      
*     9. check stopping criterion
      TMP1=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:TMP1)     
      DO I=1,NDIM
         WRK(I,1)= WRK(I,1)+ALPHA*WRK(I,3)
         WRK(I,2)= WRK(I,2)-ALPHA*WRK(I,5)    
         TMP1=TMP1+WRK(I,2)*DCONJG(WRK(I,2))
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL       
      TOLE=dsqrt(cdabs(TMP1))/NORM
      if (NLOOP.eq.50) then
         write(*,*) 'RESIDUE',TOLE,'number iteration',nloop
      endif
c      write(*,*) 'tole',cdabs(TMP1),NORM
      IF (TOLE.LE.TOL) then
         NSTAT=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)     
         DO I=1,NDIM
            XI(I)=WRK(I,1)
         ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL         
         NOU=4
         RETURN         
      ENDIF
      

*     10. s=Ar
 40   NOU=3
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)        
      DO I=1,NDIM
         XI(I)=WRK(I,2)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      
      RETURN
*     11. rdotr=dot(r,r)
 30   RHO=0.d0
*     12. delta=dot(r,s)
      DELTA=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)   REDUCTION(+:RHO,DELTA)         
      DO I=1,NDIM
         WRK(I,4)=XR(I)     
         RHO=RHO+WRK(I,2)*WRK(I,2)
         DELTA=DELTA+WRK(I,2)*WRK(I,4)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL      

*     13. beta=rdotr/rdotr0
      BETA=RHO/RHO0
c      write(*,*) 'BETA',BETA
      RHO0=RHO

*     14. p=r+beta*p
*     15. w=s+beta*w
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
!$OMP DO SCHEDULE(STATIC)    
      DO I=1,NDIM
         WRK(I,3)=WRK(I,2)+BETA*WRK(I,3)
         WRK(I,5)=WRK(I,4)+BETA*WRK(I,5)
      ENDDO
!$OMP ENDDO 
!$OMP END PARALLEL
      
*     16. xi=delta-beta^2*xi      
      XXI=DELTA-BETA*BETA*XXI

      GOTO 50


      END
