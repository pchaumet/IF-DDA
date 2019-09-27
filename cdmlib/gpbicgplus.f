      SUBROUTINE GPBICGplus(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT,TOL
     $     ,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda)

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,ETA,DZETA,R0RN,ctmp,ctmp1,ctmp2,ctmp3
     $     ,ctmp4,ctmp5
      DOUBLE PRECISION TOL,NORM,RESIDU
c      write(*,*) 'NORM',NORM
      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:r0
c     2:p
c     3:r
c     4:Ap
c     5:t
c     6:Ar
c     7:Au
c     8:u
c     9:Az
c     10:z
c     11:x
c     12:pas utilise

c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:NORM)        
      do ii=1,ndim
         WRK(ii,11)=xi(ii)
         NORM=NORM+cdabs(B(ii))**2.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      NORM=dsqrt(NORM)

      nou=1
      return

c     initialise r0,p=0,r,Ap-1=0,u-1=0,z-1=0,Au-1=0
      
 10   ctmp=0.d0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC)  REDUCTION(+:ctmp)       
      do ii=1,ndim
         WRK(ii,1)=B(ii)-xr(ii)
         WRK(ii,2)=0.d0
         WRK(ii,3)=WRK(ii,1)
         WRK(ii,4)=0.d0
         WRK(ii,7)=0.d0
         WRK(ii,8)=0.d0
         WRK(ii,10)=0.d0
         ctmp=ctmp+WRK(ii,1)*dconjg(WRK(ii,1))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      write(*,*) 'residu initial',cdabs(cdsqrt(ctmp))/norm

      R0RN=ctmp

c     initialise beta
      BETA=0.d0
      
c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1
c      write(*,*) 'ITNO',ITNO

c     calcul de p=r+beta*(p-u) et Ar
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) 
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,3)+BETA*(WRK(ii,2)-WRK(ii,8))
         xi(ii)=WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return

c     calcul de Au
 20   ctmp=0.d0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC)  REDUCTION(+:ctmp)       
      do ii=1,ndim
         WRK(ii,6)=xr(ii)
         xi(ii)=WRK(ii,8)
c     calcul de  Ap=Au+beta(Ap-Au)
         WRK(ii,4)=WRK(ii,6)+BETA*(WRK(ii,4)-WRK(ii,7))
c     calcul de alpha=r0r/r0Ap
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,4)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

      if (ctmp.eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return 
      endif

      ALPHA=R0RN/ctmp
c      write(*,*) 'ALPHA',Alpha

c     calcul des coeffs dzeta et eta

      if (ITNO.eq.0) then

        ETA=0.d0
        DZETA=0.d0
        ctmp=0.d0
c     DZETA=(Ar,Ar)/(Ar,r)
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC)  REDUCTION(+:DZETA,ctmp)            
        do ii=1,ndim
         DZETA=DZETA+dconjg(WRK(ii,6))*WRK(ii,3)
         ctmp=ctmp+dconjg(WRK(ii,6))*WRK(ii,6)
        enddo
!$OMP ENDDO 
!$OMP END PARALLEL         

        if (ctmp.eq.0.d0) then
           STATUS=-1
           STEPERR=2
           return 
        endif
        DZETA=DZETA/ctmp
c        write(*,*) 'dzeta',dzeta
      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1,ctmp2,ctmp3,ctmp4,ctmp5)          
         do ii=1,ndim
            ctmp1=ctmp1+dconjg(WRK(ii,9))*WRK(ii,9)
            ctmp2=ctmp2+dconjg(WRK(ii,6))*WRK(ii,6)
            ctmp3=ctmp3+dconjg(WRK(ii,6))*WRK(ii,3)
            ctmp4=ctmp4+dconjg(WRK(ii,9))*WRK(ii,3)
            ctmp5=ctmp5+dconjg(WRK(ii,9))*WRK(ii,6)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         ctmp=ctmp1*ctmp2-ctmp5*dconjg(ctmp5)
         
         if (ctmp.eq.0.d0) then
            STATUS=-1
            STEPERR=3
            return 
         endif

         DZETA=(ctmp1*ctmp3-ctmp4*dconjg(ctmp5))/ctmp
         ETA=(ctmp2*ctmp4-ctmp5*ctmp3)/ctmp
c         write(*,*) 'dzeta eta',dzeta,eta
      endif

!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)
c     calcul de u=dzeta*Ap+ETA*(Az+beta*u)
      do ii=1,ndim
         WRK(ii,8)=DZETA*WRK(ii,4)+ETA*(WRK(ii,9)+BETA*WRK(ii,8))
c     calcul de t=r-alpha*Ap
         WRK(ii,5)=WRK(ii,3)-ALPHA*WRK(ii,4)
c     calcul de z=dzeta*r+ETA*z-Alpha*u
         WRK(ii,10)=DZETA*WRK(ii,3)+ETA*WRK(ii,10)-ALPHA*WRK(ii,8)
c     calcul de Au
         xi(ii)=WRK(ii,8)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      nou=3
      return

 30   RESIDU=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)
      do ii=1,ndim
         WRK(ii,7)=xr(ii)
c     calcul de Az=dzeta*Ar+ETA*Az-Alpha*Au
         WRK(ii,9)=DZETA*WRK(ii,6)+ETA*WRK(ii,9)-ALPHA*WRK(ii,7)
c     calcul de x=x+alpha*p+z et r=t-Az
         WRK(ii,11)=WRK(ii,11)+ALPHA*WRK(ii,2)+WRK(ii,10)
         WRK(ii,3)=WRK(ii,5)-WRK(ii,9)
         RESIDU=RESIDU+cdabs(WRK(ii,3))**2.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL    

      RESIDU=dsqrt(RESIDU)/NORM
      if (mod(ITNO,50).EQ.0) write(*,*) 'RESIDU',RESIDU
      if (RESIDU.le.TOL) then
         STATUS=1
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)
         do ii=1,ndim
            xi(ii)=WRK(ii,11)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
         nou=4
         return
 40   endif


c     calcul de beta
      ctmp=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp) 
      do ii=1,ndim
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
      if (R0RN.eq.0.d0) then
         STATUS=-1
         STEPERR=4
         return 
      endif

      BETA=ALPHA*ctmp/DZETA/R0RN
c      write(*,*) 'beta',beta
      R0RN=ctmp



      if (ITNO.le.MAXIT) goto 100

      STATUS=1
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)
      do ii=1,ndim
         xi(ii)=WRK(ii,11)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
      END
