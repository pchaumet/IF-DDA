      subroutine polamodifsphere(nbsphere,nmax,xs,ys,zs,k0,polarisa
     $     ,aretecube,eps,eps0)
      implicit none
      integer i,j,ii,jj,kk,nbsphere,nmax
      double precision xs(nmax),ys(nmax),zs(nmax),k0,aretecube,pi
     $     ,propaesplibrestat(3,3),eps0
      double complex polarisa(nmax,3,3),eps,mattmp(3,3),mattmp1(3,3)
     $     ,icomp

      integer job,lda,info
      parameter (lda=3)
      integer ipvt(lda)
      double complex  work(lda),deterz(2)

      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      
      do i=1,nbsphere
         do ii=1,3
            do jj=1,3
               mattmp(ii,jj)=0.d0
            enddo
         enddo
         
c     calcul de chaque pola en effectuant la somme

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,ii,jj,propaesplibrestat)   
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mattmp)
         do j=1,nbsphere               
            call propa_espace_libre_stat(xs(i),ys(i),zs(i),xs(j)
     $           ,ys(j),zs(j),propaesplibrestat)
            
            do ii=1,3
               do jj=1,3
                  mattmp(ii,jj)=mattmp(ii,jj)+propaesplibrestat(ii,jj)
               enddo
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  

         
c     calcul somme de T+(eps+2)/3 et inverse
         do ii=1,3
            do jj=1,3
c     write(*,*) 'mattmp somme',mattmp(ii,jj),ii,jj,i
               mattmp(ii,jj)=mattmp(ii,jj)*aretecube*aretecube
     $              *aretecube*(eps-eps0)/4.d0/pi
            enddo
            mattmp(ii,ii)=mattmp(ii,ii)+(eps+2.d0*eps0)/3.d0/eps0
         enddo
         job=11
         call zgefa(mattmp,lda,3,ipvt,info)
         if (INFO.ne.0) then
            write(*,*) 'probleme1 in zgefa.f',info,i
            stop
         endif
         call zgedi(mattmp,lda,3,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            write(*,*) 'probleme2 in zgedi.f',info,i
            stop
         endif          
c     calcul pola claussius-mossotti: alpha0:mattmp1
         do ii=1,3
            do jj=1,3
               mattmp1(ii,jj)=mattmp(ii,jj)*aretecube*aretecube
     $              *aretecube*(eps-eps0)/4.d0/pi
            enddo
         enddo
c     calcul de 1-2/3 ik_0^3 alpha_0
         do ii=1,3
            do jj=1,3
               mattmp(ii,jj)=-mattmp1(ii,jj)*2.d0/3.d0*icomp*k0*k0
     $              *k0
               if (ii.eq.jj) mattmp(ii,jj)=1.d0+mattmp(ii,jj)
            enddo
         enddo

c     calcul inverse de   1-2/3 ik_0^3 alpha_0
         job=11
         call zgefa(mattmp,lda,3,ipvt,info)
         if (INFO.ne.0) then
            write(*,*) 'probleme3 in zgefa.f',info,i
            stop
         endif
         call zgedi(mattmp,lda,3,ipvt,deterz,work,job)      
         if (INFO.ne.0) then
            write(*,*) 'probleme4 in zgedi.f',info,i
            stop
         endif
c     calcul  de la pola: alpha_0*(1-2/3 ik_0^3 alpha_0){-1}
         do ii=1,3
            do jj=1,3
               polarisa(i,ii,jj)=0.d0
               do kk=1,3
                  polarisa(i,ii,jj)=polarisa(i,ii,jj)+mattmp1(ii
     $                 ,kk)*mattmp(kk,jj)
               enddo
            enddo
         enddo
c     fin calcul d une pola          
      enddo


      end
c******************************************************************
c     *************************************************************
c     propagateur d espace libre
      subroutine propa_espace_libre_stat(x,y,z,x0,y0,z0,propaesplibre)
      implicit none
      integer i,j
      double precision x,y,z,x0,y0,z0,k0,Id,Rab,Rtenseur,Rvect
      double precision propaesplibre,const1
      dimension propaesplibre(3,3),Id(3,3),Rtenseur(3,3),Rvect(3)

      Rab=0.d0
      Rvect(1)=(x-x0)
      Rvect(2)=(y-y0)
      Rvect(3)=(z-z0)

      do i=1,3
         do j=1,3
            Id(i,j)=0.d0
            if (i.eq.j) Id(i,i)=1.d0
            Rtenseur(i,j)=Rvect(i)*Rvect(j)
         enddo
         Rab=Rab+Rvect(i)*Rvect(i)
      enddo
      Rab=dsqrt(Rab)

c     normalise pour avoir le vecteur unitaire
      do i=1,3
         do j=1,3
            Rtenseur(i,j)=Rtenseur(i,j)/(Rab*Rab)
         enddo
      enddo
    
      if (Rab.eq.0.d0) then
         do i=1,3
            do j=1,3
               propaesplibre(i,j)=0.d0
	   enddo
         enddo
      else
         const1=Rab**(-3.d0)
         do i=1,3
            do j=1,3
               propaesplibre(i,j)=(3.d0*Rtenseur(i,j)-Id(i,j))*const1
            enddo
         enddo
      endif
      end
