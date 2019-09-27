c     routine qui fait passer du champ local au champ macro et
c     reciproquement dans le cas avec RR
      subroutine local_macro(Eloc,Em,eps,aretecube,k0,nsens)
      implicit none
      integer ii,jj,nsens
      double precision aretecube,k0,pi,const
      double complex Eloc(3),Em(3),eps(3,3)
      double complex mattemp(3,3),mattemp1(3,3),icomp

      integer job,lda,info
      parameter (lda=3)
      integer ipvt(lda)
      double complex  work(lda),deterz(2)


      icomp=(0.d0,1.d0)
      pi=dacos(-1.d0)
      const=k0*k0*k0*aretecube*aretecube*aretecube/2.d0/pi


c     champ electrique
      do ii=1,3
         do jj=1,3
            mattemp(ii,jj)=eps(ii,jj)
            mattemp1(ii,jj)=eps(ii,jj)
         enddo
         mattemp(ii,ii)=mattemp(ii,ii)+2.d0
         mattemp1(ii,ii)=mattemp1(ii,ii)-1.d0
      enddo

      do ii=1,3
         do jj=1,3
            mattemp(ii,jj)=mattemp(ii,jj)-icomp*const*mattemp1(ii,jj)
         enddo
      enddo

      if (nsens.eq.1) then
c     passe du local au champ macro
         job=11
         call zgefa(mattemp,lda,lda,ipvt,info)
         if (INFO.ne.0) then
            write(*,*) 'problem in local-macro in zgefa.f',info
            stop
         endif
         call zgedi(mattemp,lda,lda,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            write(*,*) 'problem in local-macro in zgedi.f',info
            stop
         endif

         do ii=1,3
            Em(ii)=0.d0
            do jj=1,3
               Em(ii)=Em(ii)+mattemp(ii,jj)*Eloc(jj)
            enddo
            Em(ii)=Em(ii)*3.d0
         enddo
      elseif (nsens.eq.-1) then
c     passe du champ macro au local
         do ii=1,3
            Eloc(ii)=0.d0
            do jj=1,3
               Eloc(ii)=Eloc(ii)+mattemp(ii,jj)*Em(jj)
            enddo
            Eloc(ii)=Eloc(ii)/3.d0
         enddo
      else
         write(*,*) 'mauvaise valeur de nsens',nsens
         stop
      endif
      end
