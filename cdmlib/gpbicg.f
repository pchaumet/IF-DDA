c$$$      implicit none
c$$$      integer LDA,NDIM,NLAR,NOU,STATUS,STEPERR,NLOOP,MAXIT,i,j,IM
c$$$      parameter(lda=200,nlar=12)
c$$$      double precision NORM,TOL,TOLE
c$$$      double complex XI(lda),XR(lda),B(lda),WRK(lda,nlar) ,mat(lda,lda)
c$$$     $     ,diago(lda),ALPHA,BETA,ETA,DZETA,R0RN,icomp
c$$$      
c$$$  
c$$$      ndim=5
c$$$      icomp=(0.d0,1.d0)
c$$$      do i=1,ndim
c$$$         do j=1,ndim
c$$$            mat(i,j)=0.d0
c$$$         enddo
c$$$         mat(i,i)=(1.d0,0.d0)*dble(i)*icomp
c$$$         b(i)=(1.d0,0.d0)*icomp
c$$$         diago(i)=1.d0
c$$$         xi(i)=(0.d0,0.d0)
c$$$      enddo
c$$$      
c$$$      tol=1.d-3
c$$$      tole=tol
c$$$      nloop=0
c$$$      nou=0
c$$$      MAXIT=100
c$$$      
c$$$      do i=1,ndim
c$$$         b(i)=b(i)*diago(i)
c$$$         xi(i)=xi(i)/diago(i)
c$$$      enddo
c$$$      
c$$$ 10   call GPBICG(Xi,XR,B,lda,ndim,nlar,nou,WRK,NLOOP,MAXIT ,TOLE ,NORM
c$$$     $     ,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
c$$$      
c$$$      write(*,*) 'donnees',nloop,nou,ALPHA,BETA,ETA,DZETA,R0RN
c$$$      
c$$$      write(*,*) 'TOL',tol,STATUS
c$$$      if (STATUS.lt.0) then
c$$$         write(*,*) 'stop nstat',STATUS,STEPERR
c$$$         stop
c$$$      endif
c$$$      
c$$$      do i=1,ndim
c$$$         xi(i)=xi(i)*diago(i)
c$$$      enddo
c$$$      
c$$$      do i=1,ndim
c$$$         xr(i)=0.d0
c$$$         do j=1,ndim
c$$$            xr(i)=xr(i)+mat(i,j)*xi(j)
c$$$         enddo
c$$$c         write(*,*) 'xr xi',xr(i),xi(i)
c$$$      enddo
c$$$      
c$$$      do i=1,ndim
c$$$         xr(i)=xr(i)*diago(i)
c$$$      enddo
c$$$      
c$$$      
c$$$      if (STATUS.ne.1) goto 10
c$$$  
c$$$      
c$$$      write(*,*) '********************',NLOOP
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
c$$$      if (tole.ge.tol) goto 10
c$$$
c$$$  
c$$$      end
c*************************************************
      SUBROUTINE GPBICG(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT ,TOL
     $     ,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda),zzero

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,ETA,DZETA,R0RN,ctmp,ctmp1,ctmp2,ctmp3
     $     ,ctmp4,ctmp5
      DOUBLE PRECISION TOL,NORM,RESIDU
      
      zzero=(0.d0,0.d0)

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:r0^*
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
c     12:t ancien

c     calcul norme et Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:NORM)      
      do ii=1,ndim
         WRK(ii,11)=xi(ii)
         NORM=NORM+dreal(B(ii)*dconjg(B(ii)))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      NORM=dsqrt(NORM)

      nou=1
      return

c     initialise r0=b-Ax0,rOt=rO,p0=v0=d0 et calcul residu
 10   R0RN=0.d0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC)  REDUCTION(+:R0RN)    
      do ii=1,ndim
         WRK(ii,2)=0.d0
         WRK(ii,3)=B(ii)-xr(ii)
         WRK(ii,1)=dconjg(WRK(ii,3))
         WRK(ii,5)=0.d0
         WRK(ii,9)=0.d0
         WRK(ii,8)=0.d0
         WRK(ii,10)=0.d0
         R0RN=R0RN+WRK(ii,1)*WRK(ii,3)
c         write(*,*) 'wrk',WRK(ii,3),B(ii),xr(ii),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      write(*,*) 'Initial Residue for iterative method'
     $     ,cdabs(cdsqrt(R0RN))/norm,'Norm',norm

c     initialise rho,alpha,w=1,tau=norm,theta=0,eta=0
      BETA=0.d0
      
c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1
c     write(*,*) 'ITNO',ITNO

c     calcul de p=r+beta*(p-u)
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) 
      do ii=1,ndim
         WRK(ii,2)=WRK(ii,3)+BETA*(WRK(ii,2)-WRK(ii,8))
         xi(ii)=WRK(ii,2)
c         write(*,*) 'wrk2',WRK(ii,2),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return

c     calcul de Ap et alpha=r0r/r0Ap
 20   ctmp=0.d0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC) REDUCTION(+:ctmp)        
      do ii=1,ndim
         WRK(ii,6)=xr(ii)
         ctmp=ctmp+WRK(ii,1)*WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      if (ctmp.eq.zzero) then
         STATUS=-1
         STEPERR=1
         return 
      endif

      ALPHA=R0RN/ctmp

c     calcul de y=t-r-alpha*w+alpha*Ap et de t=r-alpha AP
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC)   
      do ii=1,ndim
         WRK(ii,4)=WRK(ii,5)-WRK(ii,3)+ALPHA*(WRK(ii,6)-WRK(ii,9))
         WRK(ii,12)=WRK(ii,5)
         WRK(ii,5)=WRK(ii,3)-ALPHA*WRK(ii,6)
         xi(ii)=WRK(ii,5)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=3
      return

c     calcul de At
 30   RESIDU=0.d0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO SCHEDULE(STATIC) 
      do ii=1,ndim
         WRK(ii,7)=xr(ii)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c     calcul des coeffs dzeta et eta

      if (ITNO.eq.0) then

         ETA=0.d0
         DZETA=0.d0
         ctmp=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC) REDUCTION(+:DZETA,ctmp)         
         do ii=1,ndim
            DZETA=DZETA+dconjg(WRK(ii,7))*WRK(ii,5)
            ctmp=ctmp+dconjg(WRK(ii,7))*WRK(ii,7)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
         if (ctmp.eq.zzero) then
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
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
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
         
         if (ctmp.eq.zzero) then
            STATUS=-1
            STEPERR=3
            return 
         endif

         DZETA=(ctmp2*ctmp4-ctmp5*ctmp3)/ctmp
         ETA=(ctmp1*ctmp5-dconjg(ctmp3)*ctmp4)/ctmp

      endif

c     calcul de u, z, x et r
      
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
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
      if (mod(ITNO,20).EQ.0) write(*,*) 'RESIDUE',RESIDU,'iteration'
     $     ,ITNO
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
      endif


c     calcul de beta
 40   ctmp=0.d0
!$OMP PARALLEL   DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO  SCHEDULE(STATIC) REDUCTION(+:ctmp)  
      do ii=1,ndim
         ctmp=ctmp+WRK(ii,1)*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      if (R0RN.eq.zzero) then
         STATUS=-1
         STEPERR=4
         return 
      endif

      BETA=ALPHA*ctmp/DZETA/R0RN
      R0RN=ctmp


c     calcul de w
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO   SCHEDULE(STATIC)     
      do ii=1,ndim
         WRK(ii,9)=WRK(ii,7)+BETA*WRK(ii,6)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      if (ITNO.le.MAXIT) goto 100

      STATUS=1
      STEPERR=0
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ii)
!$OMP DO    SCHEDULE(STATIC)    
      do ii=1,ndim
         xi(ii)=WRK(ii,11)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      END
