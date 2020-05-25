      subroutine inverserigopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $     ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $     ,ntotalm,ntotal,ldabi,nlar,nmax,nxm ,nym,nzm ,nx,ny ,nz ,nx2
     $     ,ny2,nxy2,nz2,nbsphere,nbsphere3,XI,XR,wrk,FF ,FF0 ,FFloc
     $     ,polarisa ,methodeit,tol,tol1,nloop,ncompte,planf,planb
     $     ,nstop,infostr)

      
      implicit none

      integer ntotalm,ntotal,nmax,nxm,nym,nzm,nx,ny,nz,nx2,ny2,nxy2,nz2
      integer ncompte,nt,ldabi, nlar,i,ii,jj,k,nbsphere,nbsphere3,nstop
      integer NLIM,ndim,nou,nstat,nloop,STEPERR
      DOUBLE PRECISION  NORM,TOL,tol1,tole  , tmp
      double complex ALPHA,BETA,GPETA,DZETA,R0RN,QMR1,QMR2,QMR3,QMR4
     $     ,QMR5,QMR6,QMR7,QMR8,QMR9,DOTS1,DOTS2,DOTS3,DOTS4


      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc
      double complex, dimension(8*nxm*nym*nzm) :: FFTTENSORxx,
     $     FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz, FFTTENSORzz
     $     ,vectx,vecty,vectz
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa
      character(12) methodeit
      character(64) infostr
      integer*8 planf,planb
      
      nlim=10000 
      nou=0
      ndim=nbsphere3

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  
      do i=1,nbsphere3
         xi(i)=FF0(i)
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL          

      if (methodeit(1:7).eq.'GPBICG1') then
 2002    call GPBICG(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim,TOL
     $        ,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR)      
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                   
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax,nxm,nym,nzm,nx,ny,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2002
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM
c     as we begin with ITNO=-1
         nloop=nloop+1
        
c         write(*,*) 'GPBICG1 tol1',tol1,tmp/NORM,ncompte
c         write(*,*) 'FF',FF
         
      elseif (methodeit(1:7).eq.'GPBICG2') then
 2009    call GPBICG2(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)            
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                     
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)            
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo      
!$OMP ENDDO 
!$OMP END PARALLEL
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)

         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2009
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'GPBICG2 tol1',tol1,tmp/NORM
         
      elseif (methodeit(1:10).eq.'GPBICGsafe') then
 2019    call GPBICGsafe(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL             
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo  
!$OMP ENDDO 
!$OMP END PARALLEL  
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)

         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2019
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'GPBICGsafe tol1',tol1,tmp/NORM
      elseif (methodeit(1:10).eq.'GPBICGplus') then
 2016    call GPBICGplus(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL             
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         endif
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)  
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo  
!$OMP ENDDO 
!$OMP END PARALLEL  
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2016
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)  
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'GPBICGsafe tol1',tol1,tmp/NORM
         
      elseif (methodeit(1:8).eq.'GPBICGAR') then
 2010    call GPBICGAR(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo  
!$OMP ENDDO 
!$OMP END PARALLEL              
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2010
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'GPBICGAR1 tol1',tol1,tmp/NORM
         
      elseif (methodeit(1:9).eq.'GPBICGAR2') then
 2011    call GPBICGAR2(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim
     $        ,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN,NSTAT,STEPERR) 
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)               
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo 
!$OMP ENDDO 
!$OMP END PARALLEL
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2011
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'GPBICGAR2 tol1',tol1,tmp/NORM
      elseif (methodeit(1:12).eq.'BICGSTARPLUS') then
 2015    call GPBICGSTARPLUS(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK ,NLOOP
     $        ,Nlim,TOL,NORM,ALPHA,BETA,GPETA,DZETA,R0RN ,NSTAT
     $        ,STEPERR)      
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                   
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2015
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM
c     as we begin with ITNO=-1
         nloop=nloop+1
c     write(*,*) 'GPBICG1 tol1',tol1,tmp/NORM
         
         
         
      elseif (methodeit(1:6).eq.'QMRCLA') then

 2003    call PIMZQMR(FFloc,XI,XR,FF0,WRK,NORM,LDABI,NDIM,NLAR,QMR1
     $        ,QMR2,QMR3,QMR4,QMR5,QMR6,QMR7,QMR8,QMR9,DOTS1,DOTS2,DOTS3
     $        ,DOTS4,NOU,NT ,nloop,NLIM,TOLE ,TOL ,NSTAT ,STEPERR)
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif

         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         ncompte=ncompte+1
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xi(i)=FFloc(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
            nt=1
         endif
c     write(*,*) 'ncompte',ncompte,nt,tole
         if (nt.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xr(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere
               k=3*(i-1)
               xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k
     $              +2)+polarisa(i,1,3)*xr(k+3)
               xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k
     $              +2)+polarisa(i,2,3)*xr(k+3)
               xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k
     $              +2)+polarisa(i,3,3)*xr(k+3)
            enddo   
!$OMP ENDDO 
!$OMP END PARALLEL  
            call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy
     $           ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx
     $           ,vecty,vectz,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny,nz
     $           ,nx2,ny2,nxy2,nz2,XI ,XR,planf,planb)

         elseif (nt.eq.2) then
c     calcul avec le transpose
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xr(i)=0.d0
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
            call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy
     $           ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx
     $           ,vecty,vectz,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny,nz
     $           ,nx2,ny2,nxy2,nz2,XI ,XR,planf,planb)
            
c     xr=-At*xi car A=At
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k,ii,jj)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere
               k=3*(i-1)
               do ii=1,3
                  wrk(k+ii,12)=0.d0
                  do jj=1,3
                     wrk(k+ii,12)=wrk(k+ii,12)+polarisa(i,jj,ii)*xr(k
     $                    +jj)
                  enddo
                  xr(k+ii)=wrk(k+ii,12)+xi(k+ii)
               enddo               
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                    
         endif

         if (nstat.ne.1) goto  2003
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM

c     write(*,*) 'QMR tol1',tol1
         
      elseif (methodeit(1:5).eq.'TFQMR') then
 2004    call TFQMR(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK,nloop
     $        ,NLIM,TOL,NORM,QMR1,QMR2,QMR3,QMR4,QMR5,QMR6,NSTAT
     $        ,STEPERR)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xi(i)=FFloc(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy ,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI ,XR,planf,planb)


         if (nstat.ne.1) goto  2004
         
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM
         nloop=nloop+1
c     write(*,*) 'TFQMR tol1',tol1,NORM,xr(1),FF(1),FF0(1)
         
      elseif (methodeit(1:5).eq.'CG') then

 2005    call ZCG(XI,XR,FF0,NORM,WRK,QMR1,QMR2,QMR3,LDABI,NDIM,NLAR
     $        ,NOU,NSTAT,NLOOP,NLIM,TOLE,TOL)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)       
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)               
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)       
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL           
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy ,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI ,XR,planf,planb)

         if (nstat.ne.1) goto  2005
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM

c     write(*,*) 'ZCG tol1',tol1
         
      elseif (methodeit(1:8).eq.'BICGSTAB') then

 2006    call PIMZBICGSTAB(FFLOC,Xi,XR,FF0,ldabi,nlar,ndim,nou,WRK
     $        ,QMR1,QMR2,QMR3,NORM,TOL,nloop,nlim,NSTAT,STEPERR)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xi(i)=FFLOC(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL               
         endif
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo 
!$OMP ENDDO 
!$OMP END PARALLEL             

         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy ,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax  ,nxm ,nym,nzm,nx ,ny,nz,nx2,ny2,nxy2
     $        ,nz2,XI ,XR,planf,planb)

         if (nstat.ne.1) goto  2006

         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM
         
c     write(*,*) 'PIMZBICGSTAB tol1',tol1
         
      elseif (methodeit(1:12).eq.'QMRBICGSTAB1') then
         nt=1
 2007    call QMRBICGSTAB(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK
     $        ,nloop,nloop,TOL,TOLE,NORM,QMR1,QMR2,QMR3,QMR4,QMR5
     $        ,QMR6,QMR7,QMR8,QMR9,NT,NSTAT ,STEPERR)

         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1

         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i)            
!$OMP DO SCHEDULE(STATIC) 
            do i=1,nbsphere3
               xi(i)=FFLOC(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo  
!$OMP ENDDO 
!$OMP END PARALLEL            
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy ,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax  ,nxm ,nym,nzm,nx ,ny,nz,nx2,ny2,nxy2
     $        ,nz2,XI ,XR,planf,planb)
         if (nstat.ne.1) goto  2007
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM

c     write(*,*) 'QMRBICGSTAB1 tol1',tol1
         
      elseif (methodeit(1:12).eq.'QMRBICGSTAB2') then

         nt=2
 2008    call QMRBICGSTAB(FFloc,Xi,XR,FF0,ldabi,ndim,nlar,nou,WRK
     $        ,nloop,nloop,TOL,TOLE,NORM,QMR1,QMR2,QMR3,QMR4,QMR5
     $        ,QMR6,QMR7,QMR8,QMR9,NT,NSTAT ,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               xi(i)=FFLOC(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                  
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k)    
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo 
!$OMP ENDDO 
!$OMP END PARALLEL             
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy ,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz ,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax  ,nxm ,nym,nzm,nx ,ny,nz,nx2,ny2,nxy2
     $        ,nz2,XI ,XR,planf,planb)

         if (nstat.ne.1) goto  2008
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         tol1=dsqrt(tol1)/NORM

c     write(*,*) 'QMRBICGSTAB2 tol1',tol1
      elseif (methodeit(1:7).eq.'GPBICOR') then
 2013    call GPBICOR(XI,XR,FF0,FFloc,ldabi,ndim,nlar,nou,WRK
     $        ,NLOOP,Nlim,TOL,NORM,ALPHA,BETA,GPETA,DZETA,QMR1
     $        ,QMR2,NSTAT,STEPERR)
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                   
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL  
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)

         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2013
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM

c     write(*,*) 'GPBICOR tol1',tol1

      elseif (methodeit(1:7).eq.'CORS') then
 2014    call CORS(XI,XR,FF0,ldabi,ndim,nlar,nou,WRK,NLOOP,Nlim,TOL
     $        ,NORM,ALPHA,GPETA,DZETA,BETA,NSTAT,STEPERR)
         
         if (nstat.lt.0) then
            nstop=1
            infostr='Problem to solve Ax=b!'
            write(*,*) 'Problem to solve Ax=b',nstat,STEPERR
            return
         endif
         ncompte=ncompte+1
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
         do i=1,nbsphere3
            xr(i)=xi(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         if (nstat.eq.1) then
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)                
            do i=1,nbsphere3
               FFloc(i)=xi(i)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL                   
         endif
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            xi(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            xi(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotal,nmax ,nxm ,nym,nzm,nx,ny ,nz,nx2 ,ny2,nxy2
     $        ,nz2,XI,XR,planf,planb)
         if (nstop == -1) then
            infostr = 'Calculation cancelled during iterative method'
            return
         endif

         if (nstat.ne.1) goto  2014
c     compute the Residue
         tol1=0.d0
!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(i) 
!$OMP DO SCHEDULE(STATIC)  REDUCTION(+:tol1)           
         do i=1,nbsphere3
            xr(i)=xr(i)-FF0(i)
            FF(i)=xi(i)
            tol1=tol1+dreal(xr(i)*dconjg(xr(i)))
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL          
         tol1=dsqrt(tol1)/NORM
c     as we begin with ITNO=-1
         nloop=nloop+1
c         write(*,*) 'CORS tol1',tol1,tmp/NORM
         
         
      else
         write(*,*) 'Iterative method not correct',methodeit
         nstop=1
         infostr='Iterative method not correct!'
         return         
      endif


      
      end
