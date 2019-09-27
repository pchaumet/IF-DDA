c$$$      implicit none
c$$$      integer LDA,NDIM,NLAR,NOU,STATUS,STEPERR,NLOOP,MAXIT,i,j,IM
c$$$      parameter(lda=200,nlar=10)
c$$$      double precision NORM,TOL,TOLE
c$$$      double complex XI(lda),XR(lda),B(lda),WRK(lda,nlar) ,mat(lda,lda)
c$$$     $     ,diago(lda),ALPHA,RHO,RHOA,BETA,icomp
c$$$      
c$$$  
c$$$      ndim=5
c$$$      icomp=(0.d0,1.d0)
c$$$      do i=1,ndim
c$$$         do j=1,ndim
c$$$            mat(i,j)=0.d0
c$$$         enddo
c$$$         mat(i,i)=(1.d0,0.d0)*dble(i)
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
c$$$ 10   call CORS(Xi,XR,B,lda,ndim,nlar,nou,WRK,NLOOP,MAXIT ,TOLE ,NORM
c$$$     $     ,ALPHA,RHO,RHOA,BETA,STATUS,STEPERR)
c$$$      
c$$$      write(*,*) 'donnees',nloop,nou,ALPHA,BETA,RHO,RHOA
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
      SUBROUTINE CORS(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT ,TOL
     $     ,NORM,ALPHA,RHO,RHOA,BETA,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda)
      DOUBLE COMPLEX ALPHA,BETA,RHO,RHOA,ctmp
      DOUBLE PRECISION TOL,NORM,RESIDU
c      write(*,*) 'NORM',NORM
      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:r0*=Ar0 ou r0
c     2:r
c     3: rtilde=Ar
c     4:e
c     5:d
c     6:q
c     7:h
c     8:f
c     9:qtilde=Aq
c     10:x

c     calcul norme et initialise xO et calcul Ax0
      NORM=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:NORM)         
      do ii=1,ndim
         WRK(ii,10)=xi(ii)
         NORM=NORM+cdabs(B(ii))**2.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      NORM=dsqrt(NORM)

      nou=1
      return

c     initialise r0=b-Ax0 et r0* possible ici
c     calcul residu initial pour information
 10   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)    
      do ii=1,ndim
         WRK(ii,2)=B(ii)-xr(ii) 
         ctmp=ctmp+WRK(ii,2)*dconjg(WRK(ii,2))
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
      write(*,*) 'Initial Residue for iterative method'
     $     ,cdabs(cdsqrt(ctmp))/norm,NORM

c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1

c     calcul de rtilde
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)      
      do ii=1,ndim
         xi(ii)=WRK(ii,2)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
      nou=2
      return

c     calcul de r0* et rtilde
 20   if (ITNO.eq.0) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) 
         do ii=1,ndim
            WRK(ii,1)=dconjg(xr(ii))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
      endif
c     calcul de rho
      rho=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:rho)    
      do ii=1,ndim
         WRK(ii,3)=xr(ii)
         rho=rho+WRK(ii,1)*WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
      
      if (rho.eq.0.d0) then
         STATUS=-1
         STEPERR=1
         return 
      endif
c     calcul de qtilde
      if (ITNO.EQ.0) then
c     e0=r0 et d0=rtilde et q0=rtilde
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)   
         do ii=1,ndim
            WRK(ii,4)= WRK(ii,2)
            WRK(ii,5)= WRK(ii,3)
            WRK(ii,6)= WRK(ii,3)
            xi(ii)=WRK(ii,6)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL       
      else
c     calcul de BETA
         BETA=RHO/RHOA
c     calcul de e, d et q
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)     
         do ii=1,ndim
            WRK(ii,4)= WRK(ii,2)+BETA*WRK(ii,7)
            WRK(ii,5)= WRK(ii,3)+BETA*WRK(ii,8)
            WRK(ii,6)= WRK(ii,5)+BETA*(WRK(ii,8)+BETA*WRK(ii,6))
            xi(ii)=WRK(ii,6)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
      endif

      nou=3
      return

c     calcul de Aqtilde
 30   ctmp=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ctmp)
      do ii=1,ndim
         WRK(ii,9)=xr(ii)
         ctmp=ctmp+WRK(ii,1)*WRK(ii,9)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
c     calcul de alpha
      if (ctmp.eq.0.d0) then
         STATUS=-1
         STEPERR=2
         return 
      endif
      ALPHA=RHO/ctmp
      RESIDU=0.d0
c     calcul de  h,f,x,r
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)   
      do ii=1,ndim
         WRK(ii,7)=WRK(ii,4)-ALPHA*WRK(ii,6)
         WRK(ii,8)=WRK(ii,5)-ALPHA*WRK(ii,9)
         WRK(ii,10)=WRK(ii,10)+ALPHA*(2.d0*WRK(ii,4)-ALPHA*WRK(ii,6))
         WRK(ii,2)=WRK(ii,2)-ALPHA*(2.d0*WRK(ii,5)-ALPHA*WRK(ii,9))
         RESIDU=RESIDU+cdabs(WRK(ii,2))**2.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL 

      RESIDU=dsqrt(RESIDU)/NORM

      if (mod(ITNO,50).EQ.0) write(*,*) 'RESIDU',RESIDU

      if (RESIDU.le.TOL) then
         STATUS=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO SCHEDULE(STATIC)
         do ii=1,ndim
            xi(ii)=WRK(ii,10)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
         nou=4
         return
 40   endif
      RHOA=RHO

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
