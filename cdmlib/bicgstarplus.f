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
c$$$         mat(i,i)=(1.d0,0.d0)*dble(i)*(icomp**i)
c$$$         b(i)=(1.d0,0.d0)*(icomp**i)
c$$$         diago(i)=1.d0
c$$$         xi(i)=(0.d0,0.d0)
c$$$      enddo
c$$$      
c$$$      tol=1.d-3
c$$$      tole=tol
c$$$      nloop=0
c$$$      nou=0
c$$$      MAXIT=10000
c$$$      
c$$$      do i=1,ndim
c$$$         b(i)=b(i)*diago(i)
c$$$         xi(i)=xi(i)/diago(i)
c$$$      enddo
c$$$      write(*,*) '*************************'
c$$$ 10   call GPBICGSTARPLUS(Xi,XR,B,lda,ndim,nlar,nou,WRK,NLOOP,MAXIT
c$$$     $     ,TOLE ,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
c$$$      
c$$$c      write(*,*) 'donnees',nloop,nou,ALPHA,BETA,ETA,DZETA,R0RN
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
c$$$  
c$$$      end
c*************************************************
      SUBROUTINE GPBICGSTARPLUS(Xi,XR,B,lda,ndim,nlar,nou,WRK,ITNO,MAXIT
     $     ,TOL,NORM,ALPHA,BETA,ETA,DZETA,R0RN,STATUS,STEPERR)
      IMPLICIT NONE

      INTEGER ii,nou,ndim,ITNO,LDA,MAXIT,STATUS,STEPERR,nlar
      DOUBLE COMPLEX B(lda),WRK(lda,nlar),xi(lda),xr(lda),zzero

*     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA,BETA,ETA,DZETA,R0RN,ctmp,ctmp1,ctmp2,ctmp3
     $     ,ctmp4,ctmp5,ctmp6,ctmp7,ctmp8
      DOUBLE PRECISION TOL,NORM,RESIDU
      
      zzero=(0.d0,0.d0)

      IF (NOU.EQ.1) GOTO 10
      IF (NOU.EQ.2) GOTO 20
      IF (NOU.EQ.3) GOTO 30
      IF (NOU.EQ.4) GOTO 40

      STATUS = 0
      STEPERR = -1

c     1:r0^*
c     2:Ar
c     3:r
c     4:c
c     5:y
c     6:v
c     7:Ap
c     8:w
c     9:t
c     10:Ac
c     11:x
c     12:Aw

c     calcul norme et Ax0
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

c     initialise r0=b-Ax0,rOt=rO,p0=v0=d0 et calcul residu
 10   R0RN=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:R0RN)    
      do ii=1,ndim        
         WRK(ii,3)=B(ii)-xr(ii)
         WRK(ii,1)=dconjg(WRK(ii,3))
         WRK(ii,8)=0.d0
         WRK(ii,4)=0.d0
         WRK(ii,9)=0.d0
         WRK(ii,5)=0.d0
         R0RN=R0RN+WRK(ii,1)*WRK(ii,3)
c         write(*,*) 'wrk',WRK(ii,3),B(ii),xr(ii),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      write(*,*) 'Initial Residue for iterative method'
     $     ,cdabs(cdsqrt(R0RN))/norm,'Norm',norm

      
c     commence la boucle
      ITNO=-1
 100  ITNO=ITNO+1
c     write(*,*) 'ITNO',ITNO

c     calcul de Ark
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC)
      do ii=1,ndim
         xi(ii)=WRK(ii,3)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL      
      nou=2
      return


c     calcul des coeffs dzeta et eta

 20   if (ITNO.eq.0) then

         ETA=0.d0
         BETA=0.d0
         ALPHA=0.d0
         DZETA=0.d0
     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:DZETA,ALPHA,ctmp1,ctmp2)         
         do ii=1,ndim
            WRK(ii,2)=xr(ii)
            ALPHA=ALPHA+WRK(ii,1)*WRK(ii,3)
            DZETA=DZETA+dconjg(WRK(ii,2))*WRK(ii,3)
            ctmp1=ctmp1+dconjg(WRK(ii,2))*WRK(ii,2)
            ctmp2=ctmp2+WRK(ii,1)*WRK(ii,2) 
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         R0RN=ALPHA
         if (ctmp1.eq.zzero) then
            STATUS=-1
            STEPERR=2
            return 
         endif
         DZETA=DZETA/ctmp1
         if (ctmp2.eq.zzero) then
            STATUS=-1
            STEPERR=3
            return 
         endif
         ALPHA=ALPHA/ctmp2
        
         
      else
         
         ctmp1=0.d0
         ctmp2=0.d0
         ctmp3=0.d0
         ctmp4=0.d0
         ctmp5=0.d0
         ctmp6=0.d0
         ctmp7=0.d0
         ctmp8=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC)
!$OMP& REDUCTION(+:ctmp1,ctmp2,ctmp3,ctmp4,ctmp5,ctmp6,ctmp7,ctmp8)
         do ii=1,ndim
            WRK(ii,2)=xr(ii)
            ctmp1=ctmp1+WRK(ii,1)*WRK(ii,3)
            ctmp2=ctmp2+WRK(ii,1)*WRK(ii,2)            
            ctmp3=ctmp3+dconjg(WRK(ii,5))*WRK(ii,5)
            ctmp4=ctmp4+dconjg(WRK(ii,2))*WRK(ii,3)            
            ctmp5=ctmp5+dconjg(WRK(ii,2))*WRK(ii,5)
            ctmp6=ctmp6+dconjg(WRK(ii,5))*WRK(ii,3)
            ctmp7=ctmp7+dconjg(WRK(ii,2))*WRK(ii,2)
            ctmp8=ctmp8+WRK(ii,1)*WRK(ii,12)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

         if (DZETA.eq.zzero) then
            STATUS=-1
            STEPERR=4
            return 
         endif
          if (R0RN.eq.zzero) then
            STATUS=-1
            STEPERR=5
            return 
         endif
         BETA=-ALPHA/DZETA*ctmp1/R0RN

         ctmp=ctmp2-BETA*ctmp8
         if (ctmp.eq.zzero) then
            STATUS=-1
            STEPERR=6
            return 
         endif
         ALPHA=ctmp1/ctmp

         ctmp=ctmp7*ctmp3-ctmp5*dconjg(ctmp5)
         if (ctmp.eq.zzero) then
            STATUS=-1
            STEPERR=7
            return 
         endif
         
         DZETA=(ctmp3*ctmp4-ctmp5*ctmp6)/ctmp
         ETA=(ctmp7*ctmp6-dconjg(ctmp5)*ctmp4)/ctmp
         R0RN=ctmp1
      endif
c      write(*,*) 'donnees',ALPHA,BETA,ETA,DZETA,R0RN
c     calcul de Ap, v, c
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC)       
      do ii=1,ndim
         WRK(ii,7)=WRK(ii,2)-BETA*WRK(ii,12)
         WRK(ii,6)=DZETA*WRK(ii,3)+ETA*WRK(ii,9)
         WRK(ii,4)=DZETA*WRK(ii,7)+ETA*(WRK(ii,5)-BETA*WRK(ii,4))
         xi(ii)=WRK(ii,4)
c         write(*,*) 'Ap',WRK(ii,7),'v',WRK(ii,6),'c',WRK(ii,4),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      nou=3
      return

 30   RESIDU=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii)  
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RESIDU)
      do ii=1,ndim
         WRK(ii,10)=xr(ii)
         WRK(ii,8)=WRK(ii,3)-BETA*WRK(ii,8)-WRK(ii,4)
         WRK(ii,12)=WRK(ii,7)-WRK(ii,10)
         WRK(ii,9)=WRK(ii,6)-ALPHA*WRK(ii,4)
         xi(ii)=DZETA*WRK(ii,2)+ETA*WRK(ii,5)
         WRK(ii,5)=xi(ii)-ALPHA*WRK(ii,10)
         WRK(ii,11)=WRK(ii,11)+WRK(ii,6)+ALPHA*WRK(ii,8)
         WRK(ii,3)=WRK(ii,3)-xi(ii)-ALPHA*WRK(ii,12)
         RESIDU=RESIDU+dreal(WRK(ii,3)*dconjg(WRK(ii,3)))
c         write(*,*) 'AC',WRK(ii,10),'W',WRK(ii,8),'AW',WRK(ii,12),'t'
c     $        ,WRK(ii,9),'y',WRK(ii,5),'x',WRK(ii,11),'r',WRK(ii,3),ii
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      
      RESIDU=dsqrt(RESIDU)/NORM
c      write(*,*) 'RESIDU',RESIDU
      if (mod(ITNO,20).EQ.0) write(*,*) 'RESIDUE',RESIDU,'iteration'
     $     ,ITNO
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


 40   if (ITNO.le.MAXIT) goto 100

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
