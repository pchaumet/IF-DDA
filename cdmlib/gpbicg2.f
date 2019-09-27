c$$$c****************************************************************
c$$$c     Iterative solver GPBICG2
c$$$c****************************************************************
c$$$c     Authors: P. C. Chaumet and A. Rahmani
c$$$c     Date: 04/02/2010
c$$$c     
c$$$c     Purpose: iterative solver for linear system Ax=b. There is no
c$$$c     condition on the matrix. Notice that the product A x is provided
c$$$c     by the user.
c$$$c     
c$$$c     Reference: if you use this routine in your research, please
c$$$c     reference, as appropriate: P. C. Chaumet and A. Rahmani, Efficient
c$$$c     discrete dipole approximation for magneto-dielectric scatterers
c$$$c     Opt. Lett. 34, 917 (2009). J. Tang, Y. Shen, Y. Zheng, and D. Qiu,
c$$$c     Coastal Eng. 51, 143 (2004).
c$$$c
c$$$c     license: GNU GPL
c$$$c     We cannot guarantee correctness of the programs...
c$$$c****************************************************************
c$$$c****************************************************************
c$$$c****************************************************************
c$$$c     Main program and example. The main program if provided for testing
c$$$c     purposes and should be commented out to use only the routine
c$$$c     GPBICG2
c$$$      implicit none
c$$$      integer LDA,NDIM,NLAR,NOU,STATUS,STEPERR,NLOOP,MAXIT,i,j,ncompte
c$$$      parameter(lda=200,nlar=12)
c$$$      double precision NORM,TOL,TOLE
c$$$      double complex XI(lda),XR(lda),B(lda),WRK(lda,nlar) ,mat(lda,lda)
c$$$     $     ,ALPHA,BETA,ETA,DZETA,R0RN,icomp
c$$$      
c$$$c     We want to solve Ax=b
c$$$c     mat:  matrix A
c$$$
c$$$
c$$$c     lda: Input: leading dimension array of the matrix
c$$$
c$$$c     ndim: Input: dimension  of the matrix: ndim.le.lda
c$$$
c$$$c     nlar: Input: size of the work vector: nlar.ge.12
c$$$
c$$$c     MAXIT: Input: Maximum of iteration for the iterative solver
c$$$
c$$$c     NLOOP: Output: number of iteration for the iterative solver Should
c$$$c     be initialize to 0.
c$$$
c$$$c     ncompte: number of Ax products.
c$$$
c$$$c     STATUS: Output: if STATUS.lt.0 a problem occured in GPBICG
c$$$c     STATUS=1 the requested tolerance has been reached or the maximum number of
c$$$c     iterations has been reached
c$$$
c$$$c     STEPERR: Output: if STEPERR.gt.0: indicates where the problem
c$$$c     occurs in GPBICG. STEPERR=0 the maximum number of iterations has
c$$$c     been reached.  STEPERR=-1 routine completed without error
c$$$
c$$$
c$$$c     tol: Input: tolerance requested by the user. At the end we have:
c$$$
c$$$c     r=||Ax-b||/||b||<tol
c$$$
c$$$c     b: Input: right-hand side of Ax=b 
c$$$
c$$$c     norm: Output: norm of b
c$$$
c$$$c     xi: Input: initial guess; output:solution of the linear equation
c$$$
c$$$c     xr: Input: xr = A xi
c$$$
c$$$c     NOU: Local integer used by GPBICG. Should be initialized to 0.
c$$$
c$$$c     ALPHA,BETA,ETA,DZETA,R0RN: Local complex needs for GPBICG 
c$$$
c$$$c     WRK: local array used by for GPBICG
c$$$      
c$$$c     Example
c$$$      ndim=100
c$$$      icomp=(0.d0,1.d0)
c$$$      do i=1,ndim
c$$$         do j=1,ndim
c$$$            mat(i,j)=1.d0
c$$$         enddo
c$$$         mat(i,i)=(1.d0,0.d0)*dble(i)*icomp
c$$$         b(i)=(1.d0,0.d0)*icomp
c$$$         xi(i)=(0.d0,0.d0)
c$$$      enddo
c$$$      
c$$$c     initialization
c$$$      tol=1.d-6
c$$$      nloop=0
c$$$      nou=0
c$$$      MAXIT=100
c$$$      ncompte=0
c$$$   
c$$$ 10   call GPBICG2(Xi,XR,B,lda,ndim,nlar,nou,WRK,NLOOP,MAXIT ,TOL ,NORM
c$$$     $     ,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
c$$$      
c$$$      if (STATUS.lt.0) then
c$$$         write(*,*) 'stop nstat',STATUS,STEPERR
c$$$         stop
c$$$      endif
c$$$      ncompte=ncompte+1
c$$$      do i=1,ndim
c$$$         xr(i)=0.d0
c$$$         do j=1,ndim
c$$$            xr(i)=xr(i)+mat(i,j)*xi(j)
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      if (STATUS.ne.1) goto 10
c$$$      
c$$$      DO I=1,NDIM
c$$$         WRITE(*,*) 'SOL',Xi(I)
c$$$      ENDDO
c$$$
c$$$      if (STEPERR.eq.0) then
c$$$         write(*,*) 'NLOOP has reached MAXIT',NLOOP,MAXIT
c$$$      endif
c$$$
c$$$c     Compute the relative error
c$$$      tole=0.d0
c$$$      do i=1,ndim
c$$$         tole=tole+cdabs(xr(i)-b(i))**2.d0
c$$$      enddo
c$$$      tole=dsqrt(tole)/norm
c$$$      write(*,*) 'relative error',tole,'<',tol
c$$$      write(*,*) 'number of iteration',NLOOP,MAXIT
c$$$      write(*,*) 'number of product Ax',ncompte
c$$$
c$$$      end
c****************************************************************
c****************************************************************
c****************************************************************

      SUBROUTINE GPBICG2(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT ,TOL
     $     ,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda)

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,ETA,DZETA,R0RN,ctmp,ctmp1,ctmp2,ctmp3
     $     ,ctmp4,ctmp5
      DOUBLE PRECISION TOL,NORM,RESIDU

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     Column index of the various variables stored in WRK array.
c     1:r0
c     2:p
c     3:r
c     4:y
c     5:t
c     6:Ap
c     7:At
c     8:u
c     9:w
c     10:z
c     11:x
c     12:t old

c     Compute norm and Ax0 (x0 initial guess)
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)   REDUCTION(+:NORM) 
      do ii=1,ndim
         WRK(ii,11)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      NORM=dsqrt(NORM)

      nou=1
      return

c     initialize r0=b-Ax0,rOt=rO,p0=v0=d0; Compute the initial residue
 10   R0RN=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   REDUCTION(+:R0RN)    
      do ii=1,ndim
         WRK(ii,1)=B(ii)-xr(ii)
         WRK(ii,2)=0.d0
         WRK(ii,3)=WRK(ii,1)
         WRK(ii,5)=0.d0
         WRK(ii,9)=0.d0
         WRK(ii,8)=0.d0
         WRK(ii,10)=0.d0
         R0RN=R0RN+WRK(ii,1)*dconjg(WRK(ii,1))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      write(*,*) 'Initial Residue for iterative method'
     $     ,cdabs(cdsqrt(R0RN))/norm

c     initialize rho,alpha,w=1,tau=norm,theta=0,eta=0
      BETA=0.d0
      
c     begin the iteration sequence
      ITNO=-1
 100  ITNO=ITNO+1

c     compute p=r+beta*(p-u)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) 
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,3)+BETA*(WRK(ii,2)-WRK(ii,8))
         xi(ii)=WRK(ii,2)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
      nou=2
      return

c     compute Ap,alpha=r0r/r0Ap
 20   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO  SCHEDULE(STATIC)  REDUCTION(+:ctmp)        
      do ii=1,ndim
         WRK(ii,6)=xr(ii)
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      if (cdabs(ctmp).eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return 
      endif

      ALPHA=R0RN/ctmp

c     compute y=t-r-alpha*w+alpha*Ap et de t=r-alpha AP
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)  
      do ii=1,ndim
         WRK(ii,4)=WRK(ii,5)-WRK(ii,3)-ALPHA*WRK(ii,9)+ALPHA*WRK(ii,6)
         WRK(ii,12)=WRK(ii,5)
         WRK(ii,5)=WRK(ii,3)-ALPHA*WRK(ii,6)
         xi(ii)=WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      nou=3
      return

c     compute At
 30   RESIDU=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)     
      do ii=1,ndim
         WRK(ii,7)=xr(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
c     compute dzeta and eta

      if (mod(ITNO,2).eq.0) then

        ETA=0.d0
        DZETA=0.d0
        ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DZETA,ctmp)        
        do ii=1,ndim
         DZETA=DZETA+dconjg(WRK(ii,7))*WRK(ii,5)
         ctmp=ctmp+dconjg(WRK(ii,7))*WRK(ii,7)
        enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
        if (cdabs(ctmp).eq.0.d0) then
           STATUS=-1
           STEPERR=2
           return 
        endif
        DZETA=DZETA/ctmp

      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1,ctmp2,ctmp3,ctmp4,ctmp5)        
         do ii=1,ndim
            ctmp1=ctmp1+dconjg(WRK(ii,7))*WRK(ii,7)
            ctmp2=ctmp2+dconjg(WRK(ii,4))*WRK(ii,4)
            ctmp3=ctmp3+dconjg(WRK(ii,7))*WRK(ii,4)
            ctmp4=ctmp4+dconjg(WRK(ii,7))*WRK(ii,5)
            ctmp5=ctmp5+dconjg(WRK(ii,4))*WRK(ii,5)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         ctmp=ctmp1*ctmp2-ctmp3*dconjg(ctmp3)
         
         if (cdabs(ctmp).eq.0.d0) then
            STATUS=-1
            STEPERR=3
            return 
         endif

         DZETA=(ctmp2*ctmp4-ctmp5*ctmp3)/ctmp
         ETA=(ctmp1*ctmp5-dconjg(ctmp3)*ctmp4)/ctmp

      endif

c     compute u,z, x,r
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)       
      do ii=1,ndim
         WRK(ii,8)=DZETA*WRK(ii,6)+ETA*(WRK(ii,12)-WRK(ii,3)+BETA*WRK(ii
     $        ,8))
         WRK(ii,10)=DZETA*WRK(ii,3)+ETA*WRK(ii,10)-ALPHA*WRK(ii,8)
         WRK(ii,11)=WRK(ii,11)+ALPHA*WRK(ii,2)+WRK(ii,10)
         WRK(ii,3)=WRK(ii,5)-ETA*WRK(ii,4)-DZETA*WRK(ii,7)
         RESIDU=RESIDU+dreal(WRK(ii,3)*dconjg(WRK(ii,3)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      
      RESIDU=dsqrt(RESIDU)/NORM

      if (RESIDU.le.TOL) then
         STATUS=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
         do ii=1,ndim
            xi(ii)=WRK(ii,11)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
         nou=4
         return
 40   endif

c     compute beta
      ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)       
      do ii=1,ndim
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      if (cdabs(R0RN).eq.0.d0) then
         STATUS=-1
         STEPERR=4
         return 
      endif

      BETA=ALPHA*ctmp/DZETA/R0RN
      R0RN=ctmp

c     compute w
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
         WRK(ii,9)=WRK(ii,7)+BETA*WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      if (ITNO.le.MAXIT) goto 100

      STATUS=1
      STEPERR=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)    
      do ii=1,ndim
         xi(ii)=WRK(ii,11)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      END

