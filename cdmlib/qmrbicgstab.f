      SUBROUTINE QMRBICGSTAB(X,Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT
     $     ,TOL,RESIDU,NORM,ALPHA,RHO,OMEGA,TAU,TAUT,THETA,THETAT,ETA
     $     ,ETAT ,METHODE,STATUS ,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),X(lda),xi(lda),xr(lda),ctmp
      integer METHODE

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,C,ETA,RHO,RHO0,TAU,THETA,THETAT,TAUT
     $     ,ETAT,OMEGA
      DOUBLE PRECISION TOL,NORM,RESIDU

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:r
c     2:p
c     3:v
c     4:d
c     5:rtilde0
c     6:s
c     7: xtilde
c     8: dtilde
c     9: t

c     1
c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)       
      do ii=1,ndim
         x(ii)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      NORM=dsqrt(NORM)
      nou=1
      return

c     initialise r0=b-Ax0,rOt=rO,p0=v0=d0
c     calcul residu initial
 10   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)     
      do ii=1,ndim
         WRK(ii,1)=B(ii)-xr(ii)
         WRK(ii,2)=0.d0
         WRK(ii,3)=0.d0
         WRK(ii,4)=0.d0
         WRK(ii,5)=WRK(ii,1)    
         ctmp=ctmp+WRK(ii,1)*dconjg(WRK(ii,1))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
      write(*,*) 'residue initial',cdabs(cdsqrt(ctmp))/norm

c     initialise rho,alpha,w=1,tau=norm,theta=0,eta=0
      RHO=(1.d0,0.d0)
      ALPHA=(1.d0,0.d0)
      OMEGA=(1.d0,0.d0)
      TAU=cdsqrt(ctmp)
      THETA=0.d0
      ETA=0.d0

      
c     commence la boucle
      ITNO=0
 100  ITNO=ITNO+1

c     calcul de rho et beta
      RHO0=RHO
      RHO=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RHO)       
      do ii=1,ndim
         RHO=RHO+dconjg(WRK(ii,5))*WRK(ii,1)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      BETA=RHO*ALPHA/RHO0/OMEGA
c      write(*,*) 'BETA',BETA

c     calcul de p
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)        
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,1)+BETA*(WRK(ii,2)-OMEGA*WRK(ii,3))
         xi(ii)=WRK(ii,2)
c         write(*,*) 'P',WRK(ii,2),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return

c     calcul de v
 20   ALPHA=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ALPHA)       
      do ii=1,ndim
         WRK(ii,3)=xr(ii) 
         ALPHA=ALPHA+dconjg(WRK(ii,5))*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      if (cdabs(ALPHA).eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return
      endif
      ALPHA=RHO/ALPHA
c      write(*,*) 'ALPHA',ALPHA
c     calcul de s
c     calcul thetat
      THETAT=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:THETAT)       
      do ii=1,ndim
         WRK(ii,6)=WRK(ii,1)-ALPHA*WRK(ii,3)
         THETAT=THETAT+dconjg(WRK(ii,6))*WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
      IF (cdabs(TAU).eq.0.d0) then
         STATUS=-1
         STEPERR=2
         return
      endif
      THETAT=cdsqrt(THETAT)/TAU
c      write(*,*) 'THETAT',THETAT
c     calcul  de c
      c=(1.d0,0.d0)/cdsqrt(1.d0+THETAT*THETAT)
c      write(*,*) 'C',C
c     calcul de taut
      TAUT=TAU*THETAT*c
c      write(*,*) 'TAUT',TAUT
c     calcul de etat
      ETAT=c*c*ALPHA
c      write(*,*) 'ETAT',ETAT
c     calcul de dt
c     calcul de xtilde
c     calcul de t
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
         WRK(ii,8)=WRK(ii,2)+THETA*THETA*ETA/ALPHA*WRK(ii,4)
         WRK(ii,7)=x(ii)+ETAT*WRK(ii,8)
         xi(ii)=WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=3
      return

 30   OMEGA=0.d0
      CTMP=0.d0
      if (METHODE.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:OMEGA,CTMP)               
         do ii=1,ndim
            WRK(ii,9)=xr(ii)
            OMEGA=OMEGA+dconjg(WRK(ii,9))*WRK(ii,9)
            CTMP=CTMP+dconjg(WRK(ii,6))*WRK(ii,9)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
      else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:OMEGA,CTMP)              
         do ii=1,ndim
            WRK(ii,9)=xr(ii)
            OMEGA=OMEGA+dconjg(WRK(ii,6))*WRK(ii,9)
            CTMP=CTMP+dconjg(WRK(ii,6))*WRK(ii,6)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
      endif
      
      if (cdabs(OMEGA).eq.0.d0) then
         STATUS=-1
         STEPERR=3
         return
      endif
      OMEGA=CTMP/OMEGA

c      write(*,*) 'OMEGA',OMEGA
c     calcul de rk
c     calcul de THETA
      THETA=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:THETA)         
      do ii=1,ndim
         WRK(ii,1)=WRK(ii,6)-OMEGA*WRK(ii,9)     
         THETA=THETA+dconjg(WRK(ii,1))*WRK(ii,1)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
      if (cdabs(TAUT).eq.0.d0) then
         STATUS=-1
         STEPERR=4
         return
      endif
c      write(*,*) 'THETAI',THETA,TAUT
      THETA=cdsqrt(THETA)/TAUT
c      write(*,*) 'THETA',THETA
c     caclul de c
      C=(1.d0,0.d0)/cdsqrt(1.d0+THETA*THETA)
c      write(*,*) 'C',C
c     calcul de TAU
      TAU=TAUT*THETA*c
c      write(*,*) 'TAU',TAU
c     calcul de ETA
c      write(*,*) 'ETAI',c*c,OMEGA
      ETA=c*c*OMEGA
c      write(*,*) 'ETA',ETA
c     calcul de d
      if (cdabs(OMEGA).eq.0.d0) then
         STATUS=-1
         STEPERR=5
         return
      endif
c     calcul de x
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
         WRK(ii,4)=WRK(ii,6)+THETAT*THETAT*ETAT/OMEGA*WRK(ii,8)
         x(ii)=WRK(ii,7)+ETA*WRK(ii,4)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c     calcul de l'arret de la boucle
      RESIDU=dsqrt(dble(ITNO+1))*cdabs(TAU)/NORM
      if (RESIDU.le.TOL) then
         STATUS=1
         do ii=1,ndim
            xi(ii)=x(ii)
         enddo
         nou=4
         return
 40   endif

      if (ITNO.le.MAXIT) goto 100

      END
