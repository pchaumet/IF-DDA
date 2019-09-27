c$$$      implicit none
c$$$      integer LDA,NDIM,NLAR,NOU,STATUS,STEPERR,NLOOP,MAXIT,i,j,IM
c$$$      parameter(lda=200,nlar=12)
c$$$      double precision NORM,TOL,TOLE
c$$$      double complex XI(lda),XR(lda),B(lda),WRK(lda,nlar) ,mat(lda,lda)
c$$$     $     ,diago(lda),ALPHA,BETA,ETA,DZETA,icomp,rho,rhoa
c$$$     $     ,FFloc(lda)
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
c$$$      tol=1.d-10
c$$$      tole=tol
c$$$      nloop=0
c$$$      nou=0
c$$$      MAXIT=100
c$$$      
c$$$  
c$$$ 10   call GPBICOR(Xi,XR,B,FFloc,lda,ndim,nlar,nou,WRK,nloop,MAXIT ,TOL
c$$$     $     ,NORM,ALPHA,BETA,ETA,DZETA,rho,rhoa,STATUS,STEPERR)
c$$$      
c$$$c      write(*,*) 'donnees',nloop,nou,ALPHA,BETA,ETA,DZETA,R0RN
c$$$      
c$$$c      write(*,*) 'TOL',tol,STATUS
c$$$      if (STATUS.lt.0) then
c$$$         write(*,*) 'stop nstat',STATUS,STEPERR
c$$$         stop
c$$$      endif
c$$$      
c$$$      
c$$$      do i=1,ndim
c$$$         xr(i)=0.d0
c$$$         do j=1,ndim
c$$$            xr(i)=xr(i)+mat(i,j)*xi(j)
c$$$         enddo
c$$$c         write(*,*) 'xr xi',xr(i),xi(i)
c$$$      enddo
c$$$      
c$$$     
c$$$      
c$$$      
c$$$      if (STATUS.ne.1) goto 10
c$$$  
c$$$      
c$$$      write(*,*) '********************',NLOOP,tole
c$$$      DO I=1,NDIM
c$$$         WRITE(*,*) 'SOL',FFloc(I)
c$$$      ENDDO
c$$$      do i=1,ndim
c$$$         xr(i)=0.d0
c$$$         do j=1,ndim
c$$$            xr(i)=xr(i)+mat(i,j)*FFloc(j)
c$$$         enddo
c$$$c         write(*,*) 'xr xi',xr(i),xi(i)
c$$$      enddo
c$$$c     calcul erreur relative
c$$$      tole=0.d0
c$$$      do i=1,ndim
c$$$         tole=tole+cdabs(xr(i)-b(i))
c$$$      enddo
c$$$      tole=tole/norm
c$$$      write(*,*) 'tole',tole
c$$$      stop
c$$$      if (tole.ge.tol) goto 10
c$$$
c$$$  
c$$$      end
c*************************************************
      SUBROUTINE GPBICOR(Xi,XR,B,FFloc,lda,ndim,nlar,nou,WRK,ITNO,MAXIT
     $     ,TOL,NORM,ALPHA,BETA,ETA,DZETA,rho,rhoa,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda),FFloc(lda)

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,ETA,DZETA,RHO,RHOA,ctmp,ctmp1,ctmp2
     $     ,ctmp3,ctmp4,ctmp5
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
c     3:rtilde
c     4:q
c     5:w ancien
c     6:ttilde ancien
c     7:t
c     8:u
c     9:utilde
c     10:r
c     11:z
c     12:t ancien
c     x=FFloc
c     qtilde et ttilde=xr
c     y= xi 
      
c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)            
      do ii=1,ndim
         FFloc(ii)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      NORM=dsqrt(NORM)

      nou=1
      return

c     initialise r0=b-Ax0, r=r0
 10   RESIDU=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)     
      do ii=1,ndim
         WRK(ii,1)=B(ii)-xr(ii)
         WRK(ii,10)= WRK(ii,1)
         RESIDU=RESIDU+dreal(WRK(ii,1)*dconjg(WRK(ii,1)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      write(*,*) 'Initial Residue for iterative method',dsqrt(RESIDU)
     $     /norm

      
c     commence la boucle
      ITNO=0
 100  ITNO=ITNO+1
c      write(*,*) 'ITNO',ITNO

c     calcul de rtilde=Ar
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   
      do ii=1,ndim        
         xi(ii)=WRK(ii,10)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return

c     calcul de rtilde=Ar et rho=(r0,rtilde)
 20   rho=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:rho)  
      do ii=1,ndim
         WRK(ii,3)=xr(ii)
         rho=rho+dconjg(WRK(ii,1))*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
c      write(*,*) 'rho',rho
      if (rho.eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return 
      endif


      if (ITNO.eq.1) then
c     initialise p=r0,q=rtilde,t=w=0
         beta=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)     
         do ii=1,ndim
            WRK(ii,2)= WRK(ii,1)
            WRK(ii,4)= WRK(ii,3)
            WRK(ii,12)= 0.d0
            WRK(ii,5)= 0.d0
            WRK(ii,8)= 0.d0
            WRK(ii,9)= 0.d0
            WRK(ii,6)= 0.d0            
            xi(ii)=WRK(ii,4)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      else
c     calcul beta,wa=ttilde+beta*q;p=r+beta*(p-u);q=rtilde+beta*(q-u)
         beta=rho/rhoa*alpha/dzeta
c         write(*,*) 'beta',beta
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)        
         do ii=1,ndim
            WRK(ii,5)=WRK(ii,6)+beta*WRK(ii,4)
            WRK(ii,2)=WRK(ii,10)+beta*(WRK(ii,2)-WRK(ii,8))
            WRK(ii,4)=WRK(ii,3)+beta*(WRK(ii,4)-WRK(ii,9))
            xi(ii)=WRK(ii,4)
c            write(*,*) 'w', WRK(ii,5),'p',WRK(ii,2),'q',WRK(ii,4),ii
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL        
      endif

      nou=3
      return

c     calcul a de A*q=xr

c     Calcul alpha=rho/(r0,qtilde)
 30   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)         
      do ii=1,ndim
         ctmp=ctmp+dconjg(WRK(ii,1))*xr(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      if (ctmp.eq.0.d0) then
         STATUS=-1
         STEPERR=2
         return 
      endif

      ALPHA=RHO/ctmp
c      write(*,*) 'alpha',alpha,rho,ctmp
c     y=y-r+alpha*(q-w);t=r-alpha*q;ttilde=rtilde-alpha*qtilde
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(STATIC)
      do ii=1,ndim
         xi(ii)=WRK(ii,12)-WRK(ii,10)+alpha*(WRK(ii,4)-WRK(ii,5))
         WRK(ii,7)=WRK(ii,10)-alpha*WRK(ii,4)
         xr(ii)=WRK(ii,3)-alpha*xr(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
      

c     calcul des coeffs dzeta et eta

      if (ITNO.eq.1) then

        ETA=0.d0
        DZETA=0.d0
        ctmp=0.d0
c     dzeta=(ttilde,t)/(ttilde,ttilde)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)     
        do ii=1,ndim
         DZETA=DZETA+dconjg(xr(ii))*WRK(ii,7)
         ctmp=ctmp+dconjg(xr(ii))*xr(ii)
        enddo
!$OMP ENDDO 
!$OMP END PARALLEL
c        write(*,*) 'dzeta',DZETA,ctmp
        if (ctmp.eq.0.d0) then
           STATUS=-1
           STEPERR=3
           return 
        endif
        DZETA=DZETA/ctmp
      
      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
c
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp1,ctmp2,ctmp3,ctmp4,ctmp5) 
         do ii=1,ndim
            ctmp1=ctmp1+dconjg(xr(ii))*xr(ii)
            ctmp2=ctmp2+dconjg(xi(ii))*WRK(ii,7)
            ctmp3=ctmp3+dconjg(xi(ii))*xr(ii)
            ctmp4=ctmp4+dconjg(xr(ii))*WRK(ii,7)
            ctmp5=ctmp5+dconjg(xi(ii))*xi(ii)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         ctmp=ctmp1*ctmp5-ctmp3*dconjg(ctmp3)
         if (ctmp.eq.0.d0) then
            STATUS=-1
            STEPERR=4
            return 
         endif
         ETA=(ctmp1*ctmp2-ctmp3*ctmp4)/ctmp
         DZETA=(ctmp5*ctmp4-ctmp2*dconjg(ctmp3))/ctmp


      endif
c      write(*,*) 'ETA',ETA,DZETA
c     calcul de u=dzeta*q+eta*(t-r+beta*u)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)     
      do ii=1,ndim
         WRK(ii,8)=dzeta*WRK(ii,4)+eta*(WRK(ii,12)-WRK(ii,10)+beta
     $        *WRK(ii,8))
         WRK(ii,9)=dzeta*(WRK(ii,3)-xr(ii))/alpha+eta*(WRK(ii,6)-WRK(ii
     $        ,3)+beta*WRK(ii,9))
         WRK(ii,11)=dzeta*WRK(ii,10)+eta*WRK(ii,11)-alpha*WRK(ii,8)
c         write(*,*) 'u',WRK(ii,8),'util',WRK(ii,9),'z',WRK(ii,11),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

c     calcul de x et r
c     passe ancien a nouveau            
      RESIDU=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+: RESIDU) 
      do ii=1,ndim
         WRK(ii,10)=WRK(ii,7)-eta*xi(ii)-dzeta*xr(ii)
         RESIDU=RESIDU+dreal(WRK(ii,10)*dconjg(WRK(ii,10)))
         FFloc(ii)=FFloc(ii)+alpha*WRK(ii,2)+WRK(ii,11)
         WRK(ii,6)=xr(ii)
         WRK(ii,12)=WRK(ii,7)
c         write(*,*) 'sol',FFloc(ii),'r',WRK(ii,10),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      RESIDU=dsqrt(RESIDU)/NORM

      if (mod(ITNO,50).EQ.0) write(*,*) 'RESIDUE',RESIDU,'iteration'
     $     ,itno
      if (RESIDU.le.TOL) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)
         do ii=1,ndim
            xi(ii)=FFloc(ii)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
         STATUS=1        
         nou=4
         return
      endif
 40   RHOA=RHO

      if (ITNO.le.MAXIT) goto 100

      STATUS=1

      END
