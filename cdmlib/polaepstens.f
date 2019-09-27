      subroutine polaepstens(rayon1,eps,eps0,k0,dis,methode,inv
     $     ,polaeps)
      implicit none
      integer dis,ii,jj,kk,inv
      double precision rayon,volume,pi,k0,k03,eps0,rayon1,indice
      double complex eps(3,3),icomp,polaeps(3,3)
     $     ,polaeps0(3,3),mattemp(3,3),mattemp1(3,3)
      character*2 methode

      integer job,lda,info,nlda
      parameter (lda=3)
      integer ipvt(lda)
      double complex  work(lda),deterz(2)

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      indice=dsqrt(eps0)
      k03=k0*k0*k0*indice*indice*indice
      nlda=3
      rayon=rayon1

      if (dis.eq.0) then
         volume=rayon*rayon*rayon
      elseif (dis.eq.1) then
         rayon=rayon*((0.75d0/pi)**(1.d0/3.d0))
         volume=rayon*rayon*rayon
      else
         write(*,*) 'Bad value of  dis',dis
         stop
      endif

c     calcul des polas statiques
      do ii=1,3
         do jj=1,3
            mattemp(ii,jj)=eps(ii,jj)
            mattemp1(ii,jj)=eps(ii,jj)
         enddo
         mattemp(ii,ii)=mattemp(ii,ii)+2.d0*eps0
         mattemp1(ii,ii)=mattemp1(ii,ii)-eps0
      enddo

      job=11
      call zgefa(mattemp,lda,nlda,ipvt,info)
      if (INFO.ne.0) then
         write(*,*) 'problem in polatens1 in zgefa.f',info
         stop
      endif
      call zgedi(mattemp,lda,nlda,ipvt,deterz,work,job)      
      if (INFO.ne.0) then
         write(*,*) 'problem in polatens2 in zgedi.f',info
         stop
      endif
      do ii=1,3
         do jj=1,3
            polaeps0(ii,jj)=0.d0
            do kk=1,3
               polaeps0(ii,jj)=polaeps0(ii,jj)+mattemp1(ii,kk)
     $              *mattemp(kk,jj)
            enddo
            polaeps0(ii,jj)=polaeps0(ii,jj)*volume
         enddo
      enddo

      if (methode.eq.'CM') then
         do ii=1,3
            do jj=1,3
               polaeps(ii,jj)=polaeps0(ii,jj)*eps0
            enddo
         enddo
      elseif (methode.eq.'RR') then
c     calcul de I-(2/3)ik^3 alpha_0
         do ii=1,3
            do jj=1,3
               mattemp(ii,jj)=-2.d0/3.d0*icomp*k03*polaeps0(ii,jj)
            enddo
            mattemp(ii,ii)=mattemp(ii,ii)+1.d0
         enddo
c      calcul de (I-(2/3)ik^3 alpha_0)^{-1}
         job=11
         call zgefa(mattemp,lda,nlda,ipvt,info)
         if (INFO.ne.0) then
            write(*,*) 'problem in polatens5 in zgefa.f',info
            stop
         endif
         call zgedi(mattemp,lda,nlda,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            write(*,*) 'problem in polatens6 in zgedi.f',info
            stop
         endif
c     alpha= (I-(2/3)ik^3 alpha_0)^{-1}*alpha_0
         do ii=1,3
            do jj=1,3
               polaeps(ii,jj)=0.d0
               do kk=1,3
                  polaeps(ii,jj)=polaeps(ii,jj)+mattemp(ii,kk)
     $                 *polaeps0(kk,jj)*eps0
               enddo
            enddo
         enddo
      else
         write(*,*) 'incorrect method for the polarizability '
         write(*,*) 'method',methode
         stop
      endif

      if (inv.eq.-1) then
         job=11
         call zgefa(polaeps,lda,nlda,ipvt,info)
         if (INFO.ne.0) then
            write(*,*) 'problem in polatens9 in zgefa.f',info
            stop
         endif
         call zgedi(polaeps,lda,nlda,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            write(*,*) 'problem in polatens10 in zgedi.f',info
            stop
         endif
      endif



      end
