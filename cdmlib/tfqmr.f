      SUBROUTINE TFQMR(X,Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT ,TOL
     $     ,NORM,ALPHA,BETA,RHO,TAU,THETA,ETA,STATUS,STEPERR)
      IMPLICIT NONE

*           PIM -- The Parallel Iterative Methods package
*           ---------------------------------------------
*
*                      Rudnei Dias da Cunha
*     National Supercomputing Centre and Mathematics Institute
*         Universidade Federal do Rio Grande do Sul, Brasil
*
*                          Tim Hopkins
*     Computing Laboratory, University of Kent at Canterbury, U.K.
*
* ----------------------------------------------------------------------
*
    
*     .. Array Arguments ..
      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),X(lda),xi(lda),xr(lda)
 
*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,C,ETA,RHO,RHO0,TAU,THETA,THETAtmp
      DOUBLE PRECISION TOL,NORM

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:u
c     2:w
c     3:r
c     4:v
c     5:d
c     6:rtilde
c     7: u+1
c     8 Au
c
c     1
c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   REDUCTION(+:NORM)            
      do ii=1,ndim
         x(ii)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      NORM=dsqrt(NORM)
      nou=1
      return

c     initialise w0=u0=r0=b-Ax0
 10   nou=2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
         WRK(ii,3)=B(ii)-xr(ii)
         WRK(ii,2)=WRK(ii,3)
         WRK(ii,1)=WRK(ii,3)
         xi(ii)=WRK(ii,1)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      return
c     initialise v0 et d0
 20   RHO=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RHO)    
      do ii=1,ndim
         WRK(ii,4)=xr(ii)
         WRK(ii,5)=0.d0
         WRK(ii,8)=WRK(ii,4)
         WRK(ii,6)=WRK(ii,3)
         RHO=RHO+dconjg(WRK(ii,6))*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
c     2
c     initialise tau,theta,eta
      TAU=(1.d0,0.d0)*NORM
      THETA=0.d0
      ETA=0.d0

c     3
c     defini rtilde

c     4
c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1

c     5
c     regarde si boucle impair (even pair)
      if (mod(ITNO,2).eq.0) then
c     6
c     calcul alpha
         ALPHA=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ALPHA)       
         do ii=1,ndim
            ALPHA=ALPHA+dconjg(WRK(ii,6))*WRK(ii,4)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
         ALPHA=RHO/ALPHA
c     7
c     calul de u+1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)         
         do ii=1,ndim
            WRK(ii,7)=WRK(ii,1)-ALPHA*WRK(ii,4)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL       
      endif
c     8
c     fin condition
c     9
c     calcul w, calcul de d,   calcul de THETA
      THETAtmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:THETATMP)        
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,2)-ALPHA*WRK(ii,8)
         WRK(ii,5)=WRK(ii,1)+THETA*THETA/ALPHA*ETA*WRK(ii,5)
         THETAtmp=THETAtmp+dconjg(WRK(ii,2))*WRK(ii,2)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
c     11
c     et C
      THETA=cdsqrt(THETAtmp)/TAU
      C=(1.d0,0.d0)/cdsqrt(1.d0+THETA*THETA)
c     12
c     calcul de TAU et ETA
      TAU=TAU*THETA*C
      ETA=C*C*ALPHA
c     13
c     calcul de X
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)     
      do ii=1,ndim
         x(ii)=x(ii)+ETA*WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     

c     calcul de l'arret de la boucle
      if (dsqrt(dble(ITNO+1))*cdabs(TAU)/NORM.le.TOL) then
         STATUS=1       
         nou=4
         return
      endif

c     14
 40   if (mod(ITNO,2).eq.1) then
c     15
         RHO0=RHO
         RHO=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:RHO)     
         do ii=1,ndim
            RHO=RHO+dconjg(WRK(ii,6))*WRK(ii,2)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
         BETA=RHO/RHO0
c     16
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)        
         do ii=1,ndim
            WRK(ii,7)=WRK(ii,2)+BETA*WRK(ii,1)
            WRK(ii,9)=WRK(ii,8)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
      endif

c     calcul de Au
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)    
      do ii=1,ndim
         xi(ii)=WRK(ii,7)
         WRK(ii,1)=WRK(ii,7)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL       
      nou=3
      return
         
 30   nou=3
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   
      do ii=1,ndim
         WRK(ii,8)=xr(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      if (mod(ITNO,2).eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   
         do ii=1,ndim
            WRK(ii,4)=WRK(ii,8)+BETA*(WRK(ii,9)+BETA*WRK(ii,4))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
      endif
      
      if (ITNO.le.MAXIT) goto 100

      END
