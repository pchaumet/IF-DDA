c$$$      implicit none
c$$$      integer LDA,NDIM,NLAR,NOU,STATUS,STEPERR,NLOOP,MAXIT,i,j,IM
c$$$      parameter(lda=200,nlar=12)
c$$$      double precision NORM,TOL,TOLE
c$$$      double complex XI(lda),XR(lda),B(lda),WRK(lda,nlar) ,mat(lda,lda)
c$$$     $     ,diago(lda),ALPHA,BETA,ETA,DZETA,R0RN,icomp
c$$$      
c$$$  
c$$$      ndim=100
c$$$      icomp=(0.d0,1.d0)
c$$$      do i=1,ndim
c$$$         do j=1,ndim
c$$$            mat(i,j)=0.d0
c$$$         enddo
c$$$         mat(i,i)=(1.d0,0.d0)*dble(i)
c$$$         b(i)=(1.d0,0.d0)
c$$$         xi(i)=(0.d0,0.d0)
c$$$      enddo
c$$$      
c$$$      tol=1.d-6
c$$$      tole=tol
c$$$      nloop=0
c$$$      nou=0
c$$$      MAXIT=100
c$$$  
c$$$ 10   call GPBICGsafe(Xi,XR,B,lda,ndim,nlar,nou,WRK,NLOOP,MAXIT ,TOLE
c$$$     $     ,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
c$$$      write(*,*) '****************************************'
c$$$      write(*,*) 'donnees',nloop,nou,ALPHA,BETA,ETA,DZETA,R0RN
c$$$      
c$$$      write(*,*) 'TOL',tol,STATUS
c$$$      if (STATUS.lt.0) then
c$$$         write(*,*) 'stop nstat',STATUS,STEPERR
c$$$         stop
c$$$      endif
c$$$      
c$$$      do i=1,ndim
c$$$         xr(i)=0.d0
c$$$         do j=1,ndim
c$$$            xr(i)=xr(i)+mat(i,j)*xi(j)
c$$$         enddo
c$$$         write(*,*) 'xr xi',xr(i),xi(i)
c$$$      enddo
c$$$ 
c$$$      if (STATUS.ne.1) goto 10
c$$$  
c$$$      
c$$$      write(*,*) '********************NLOOP',NLOOP
c$$$      DO I=1,NDIM
c$$$         WRITE(*,*) 'SOL',Xi(I)
c$$$      ENDDO
c$$$
c$$$c     calcul erreur relative
c$$$      tole=0.d0
c$$$      do i=1,ndim
c$$$         tole=tole+cdabs(xr(i)-b(i))
c$$$      enddo
c$$$      tole=tole/norm
c$$$      write(*,*) 'tole',tole
c$$$
c$$$c      if (tole.ge.tol) goto 10
c$$$
c$$$  
c$$$      end
c*************************************************
      SUBROUTINE GPBICGsafe(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT
     $     ,TOL,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
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
c     3:u
c     4:r
c     5:Au
c     6:Ar
c     7:Ap
c     8:t
c     9:z
c     10:Az
c     11:x

c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)        
      do ii=1,ndim
         WRK(ii,11)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      NORM=dsqrt(NORM)

      nou=1
      return

c     initialise r0=b-Ax0,p=u=z=0 Ap=Au=Az=0 et rn=r0
 10   R0RN=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:R0RN)         
      do ii=1,ndim
         WRK(ii,1)=B(ii)-xr(ii)
         WRK(ii,2)=0.d0
         WRK(ii,3)=0.d0
         WRK(ii,9)=0.d0
         WRK(ii,7)=0.d0
         WRK(ii,5)=0.d0
         WRK(ii,10)=0.d0
         WRK(ii,4)=WRK(ii,1)
         R0RN=R0RN+WRK(ii,1)*dconjg(WRK(ii,1))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
      
c     calcul residu initial
   
      write(*,*) 'Initial Residue for iterative method'
     $     ,cdabs(cdsqrt(R0RN))/norm


c     initialise rho,alpha,w=1,tau=norm,theta=0,eta=0
      BETA=0.d0
      
c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1
c      write(*,*) 'ITNO',ITNO

c     calcul de p=r+beta*(p-u) et r(4) 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,4)+BETA*(WRK(ii,2)-WRK(ii,3))
         xi(ii)=WRK(ii,4)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL        
      nou=2
      return

c     calcul de Ar(6) et Ap (7)=Ar+beta*(Ap-Au) et calcul de alpha=r0r/r0Ap
 20   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)        
      do ii=1,ndim
         WRK(ii,6)=xr(ii)
         WRK(ii,7)=WRK(ii,6)+beta*(WRK(ii,7)-WRK(ii,5))      
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,7)
c         write(*,*) 'Ar',WRK(ii,6),'Ap',WRK(ii,7)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      if (ctmp.eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return 
      endif

      ALPHA=R0RN/ctmp
c      write(*,*) 'alpha',alpha
c     calcul des coeffs dzeta et eta

      if (ITNO.eq.0) then

        ETA=0.d0
        DZETA=0.d0
        ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp,DZETA)         
        do ii=1,ndim
         DZETA=DZETA+dconjg(WRK(ii,6))*WRK(ii,4)
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
c        write(*,*) 'DZETA',DZETA
      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1,ctmp2,ctmp3,ctmp4,ctmp5)           
         do ii=1,ndim
            ctmp1=ctmp1+dconjg(WRK(ii,10))*WRK(ii,10)
            ctmp2=ctmp2+dconjg(WRK(ii,6))*WRK(ii,4)
            ctmp3=ctmp3+dconjg(WRK(ii,10))*WRK(ii,4)
            ctmp4=ctmp4+dconjg(WRK(ii,6))*WRK(ii,10)
            ctmp5=ctmp5+dconjg(WRK(ii,6))*WRK(ii,6)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         ctmp=ctmp5*ctmp1-ctmp4*dconjg(ctmp4)
         
         if (ctmp.eq.0.d0) then
            STATUS=-1
            STEPERR=3
            return 
         endif

         DZETA=(ctmp1*ctmp2-ctmp3*ctmp4)/ctmp
         ETA=(ctmp5*ctmp3-dconjg(ctmp4)*ctmp2)/ctmp
c         write(*,*) 'DZETA',DZETA,ETA
      endif

c     calcul de u=dzeta*Ap+eta*(Az+beta*u)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)        
      do ii=1,ndim
         WRK(ii,3)=DZETA*WRK(ii,7)+ETA*(WRK(ii,10)+BETA*WRK(ii,3))
         xi(ii)=WRK(ii,3)
c         write(*,*) 'u',WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      nou=3
      return

c     Calcul de t,z,Az,x,r
 30   RESIDU=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)          
      do ii=1,ndim
c     Au (5)
         WRK(ii,5)=xr(ii)
c     t=r-alpha*Ap
         WRK(ii,8)= WRK(ii,4)-alpha*WRK(ii,7)
c         write(*,*) 't',WRK(ii,8)
c     z=dzeta r+eta  z -alpha u
         WRK(ii,9)= DZETA*WRK(ii,4)+ETA*WRK(ii,9)-alpha*WRK(ii,3)
c         write(*,*) 'z',WRK(ii,9)
c     Az=dzeta Ar+eta  Az -alpha Au
         WRK(ii,10)= DZETA*WRK(ii,6)+ETA*WRK(ii,10)-alpha*WRK(ii,5)
c         write(*,*) 'Az',WRK(ii,10)
c     x=x+alpha p+ z
         WRK(ii,11)=WRK(ii,11)+alpha*WRK(ii,2)+WRK(ii,9)
c         write(*,*) 'x',WRK(ii,11)
c     r=t-Az
         WRK(ii,4)=WRK(ii,8)-WRK(ii,10)
c         write(*,*) 'r',WRK(ii,4)
         RESIDU=RESIDU+dreal(WRK(ii,4)*dconjg(WRK(ii,4)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      RESIDU=dsqrt(RESIDU)/NORM
      if (mod(ITNO,50).EQ.0) write(*,*) 'RESIDUE',RESIDU,'iteration'
     $     ,itno
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
      endif


c     calcul de beta
 40   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)          
      do ii=1,ndim
         ctmp=ctmp+dconjg(WRK(ii,1))*WRK(ii,4)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      if (R0RN.eq.0.d0) then
         STATUS=-1
         STEPERR=4
         return 
      endif

      BETA=ALPHA*ctmp/DZETA/R0RN
c     write(*,*) 'beta',beta
      R0RN=ctmp

      if (ITNO.le.MAXIT) goto 100

      STATUS=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)             
      do ii=1,ndim
         xi(ii)=WRK(ii,11)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      END
