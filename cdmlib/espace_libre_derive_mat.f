c     propagateur d espace libre en utilisant la symetrie
      subroutine propesplibdermat(x,y,z,x0,y0,z0,k0,test,Stenseurd)
      implicit none
      integer test
      double precision x,y,z,x0,y0,z0,k0,a,rab,a2,rab2,rab3,rab4,rab5
      double complex const1,const2,Txx,Txy,Tzz,Txz,uncomp,icomp
      double complex dconst1,dconst2,c13c2,dc13c2,cdexpo,ddd
      double complex Txxdz,Txydz,Txzdz,Tzzdz,Txxda,Txyda,Txzda,Tzzda
      double complex Tzxdz,Tzxda,Tzx
      double complex Stenseurd(3,3,3)
      double precision aa,cphi,sphi,s2phi,c2phi,xxp,yyp


      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      xxp=x-x0
      yyp=y-y0

      a=dsqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))
      a2=a*a
      Rab=dsqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
      Rab2=Rab*Rab
      Rab3=Rab2*Rab
      Rab4=Rab3*Rab
      Rab5=Rab4*Rab 

      if (Rab.eq.0.d0) then
         Txx=0.d0
         Txy=0.d0
         Txz=0.d0
         Tzz=0.d0
         Txxdz=0.d0
         Txydz=0.d0
         Txzdz=0.d0
         Tzzdz=0.d0
         Txxda=0.d0
         Txyda=0.d0
         Txzda=0.d0
         Tzzda=0.d0
      else
         const1=uncomp/Rab3-icomp*k0/Rab2
         const2=uncomp*k0*k0/Rab
         c13c2=3.d0*const1-const2      
         dconst1=-(3.d0,0.d0)/Rab5+(0.d0,2.d0)*k0/Rab4
         dconst2=-uncomp*k0*k0/Rab3
         dc13c2=3.d0*dconst1-dconst2
c     derive de c13c2/Rab2
         ddd=(dc13c2*Rab2-c13c2*2.d0)/Rab4
         cdexpo=cdexp((0.d0,1.d0)*k0*Rab)
         
         Txx=((const2-const1)+a2*c13c2/Rab2/2.d0)*cdexpo
         
         Tzz=((const2-const1)+(z-z0)*(z-z0)*c13c2/Rab2)*cdexpo
         
         Txy=-a2*c13c2/Rab2/2.d0*cdexpo
         
         Txz=a*(z-z0)*c13c2/Rab2*cdexpo
         
         Txxdz=icomp*Txx*(z-z0)*k0/Rab+cdexpo*(
     *        (dconst2-dconst1)*(z-z0)+a2/2.d0*ddd*(z-z0))
         
         Tzzdz=icomp*Tzz*(z-z0)*k0/Rab+cdexpo*(
     *        (dconst2-dconst1)*(z-z0)+2.d0*(z-z0)*c13c2/Rab2+
     *        (z-z0)*(z-z0)*ddd*(z-z0))
         
         Txydz=icomp*Txy*(z-z0)*k0/Rab-a2*cdexpo
     *        /2.d0*ddd*(z-z0)

         Txzdz=icomp*Txz*(z-z0)*k0/Rab+a*cdexpo*(
     *        c13c2/Rab2+(z-z0)*ddd*(z-z0))
         
         Txxda=icomp*Txx*a*k0/Rab+cdexpo*(
     *        (dconst2-dconst1)*a+a*c13c2/Rab2+a2/2.d0*ddd*a)
         
         Tzzda=icomp*Tzz*a*k0/Rab+cdexpo*(
     *        (dconst2-dconst1)*a+(z-z0)*(z-z0)*ddd*a)
         
         Txyda=icomp*Txy*a*k0/Rab+cdexpo*(
     *     -a*c13c2/Rab2-a2/2.d0*ddd*a)
         
         Txzda=icomp*Txz*a*k0/Rab+cdexpo*(
     *        (z-z0)*c13c2/Rab2+a*(z-z0)*ddd*a)
         
      endif

      Tzx=Txz
      Tzxda=Txzda
      Tzxdz=Txzdz
      aa=a
      
c     reconstruit le tenseur derivee
c     derivee suivant l axe X
      if (test.eq.1) then
         if (a.le.1.d-30) then
            sphi=1.d0
            cphi=0.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,1)=
     *        sphi*Txxda-2.d0*cphi*s2phi*Txy/aa+c2phi*sphi*Txyda
         Stenseurd(1,2,1)=
     *        -2.d0*cphi*c2phi*Txy/aa-s2phi*sphi*Txyda
         Stenseurd(1,3,1)=cphi*cphi*Txz/aa+sphi*sphi*Txzda
         Stenseurd(2,1,1)=
     *        -2.d0*cphi*c2phi*Txy/aa-s2phi*sphi*Txyda  
         Stenseurd(2,2,1)=
     *        sphi*Txxda+2.d0*cphi*s2phi*Txy/aa-c2phi*sphi*Txyda
         Stenseurd(2,3,1)=-sphi*cphi*Txz/aa+cphi*sphi*Txzda
         Stenseurd(3,1,1)=cphi*cphi*Tzx/aa+sphi*sphi*Tzxda
         Stenseurd(3,2,1)=-sphi*cphi*Tzx/aa+cphi*sphi*Tzxda
         Stenseurd(3,3,1)=sphi*Tzzda
      elseif(test.eq.2) then
c     derivee suivant l axe Y
         aa=a
         if (a.le.1.d-30) then
            sphi=0.d0
            cphi=1.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,2)=
     *        cphi*Txxda+2.d0*sphi*s2phi*Txy/aa+c2phi*cphi*Txyda
         Stenseurd(1,2,2)=
     *        2.d0*sphi*c2phi*Txy/aa-s2phi*cphi*Txyda  
         Stenseurd(1,3,2)=-sphi*cphi*Txz/aa+sphi*cphi*Txzda
         Stenseurd(2,1,2)=
     *        2.d0*sphi*c2phi*Txy/aa-s2phi*cphi*Txyda
         Stenseurd(2,2,2)=
     *        cphi*Txxda-2.d0*sphi*s2phi*Txy/aa-c2phi*cphi*Txyda
         Stenseurd(2,3,2)=sphi*sphi*Txz/aa+cphi*cphi*Txzda
         Stenseurd(3,1,2)=-sphi*cphi*Tzx/aa+sphi*cphi*Tzxda
         Stenseurd(3,2,2)=sphi*sphi*Tzx/aa+cphi*cphi*Tzxda
         Stenseurd(3,3,2)=cphi*Tzzda
         
      elseif(test.eq.3) then
         
c     derivee suivant l axe Z
         aa=a
         if (a.le.1.d-30) then
            cphi=0.d0
            sphi=0.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi 
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,3)=Txxdz+c2phi*Txydz
         Stenseurd(1,2,3)=-s2phi*Txydz
         Stenseurd(1,3,3)=sphi*Txzdz
         Stenseurd(2,1,3)=-s2phi*Txydz
         Stenseurd(2,2,3)=Txxdz-c2phi*Txydz
         Stenseurd(2,3,3)=cphi*Txzdz
         Stenseurd(3,1,3)=sphi*Tzxdz
         Stenseurd(3,2,3)=cphi*Tzxdz
         Stenseurd(3,3,3)=Tzzdz
      elseif(test.eq.4) then
      
c     reconstruit le tenseur derivee
c     derivee suivant l axe X
         if (a.le.1.d-30) then
            sphi=1.d0
            cphi=0.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,1)=
     *        sphi*Txxda-2.d0*cphi*s2phi*Txy/aa+c2phi*sphi*Txyda
         Stenseurd(1,2,1)=
     *        -2.d0*cphi*c2phi*Txy/aa-s2phi*sphi*Txyda
         Stenseurd(1,3,1)=cphi*cphi*Txz/aa+sphi*sphi*Txzda
         Stenseurd(2,1,1)=
     *        -2.d0*cphi*c2phi*Txy/aa-s2phi*sphi*Txyda  
         Stenseurd(2,2,1)=
     *        sphi*Txxda+2.d0*cphi*s2phi*Txy/aa-c2phi*sphi*Txyda
         Stenseurd(2,3,1)=-sphi*cphi*Txz/aa+cphi*sphi*Txzda
         Stenseurd(3,1,1)=cphi*cphi*Tzx/aa+sphi*sphi*Tzxda
         Stenseurd(3,2,1)=-sphi*cphi*Tzx/aa+cphi*sphi*Tzxda
         Stenseurd(3,3,1)=sphi*Tzzda
         
c     derivee suivant l axe Y
         aa=a
         if (a.le.1.d-30) then
            sphi=0.d0
            cphi=1.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,2)=
     *        cphi*Txxda+2.d0*sphi*s2phi*Txy/aa+c2phi*cphi*Txyda
         Stenseurd(1,2,2)=
     *        2.d0*sphi*c2phi*Txy/aa-s2phi*cphi*Txyda  
         Stenseurd(1,3,2)=-sphi*cphi*Txz/aa+sphi*cphi*Txzda
         Stenseurd(2,1,2)=
     *        2.d0*sphi*c2phi*Txy/aa-s2phi*cphi*Txyda
         Stenseurd(2,2,2)=
     *        cphi*Txxda-2.d0*sphi*s2phi*Txy/aa-c2phi*cphi*Txyda
         Stenseurd(2,3,2)=sphi*sphi*Txz/aa+cphi*cphi*Txzda
         Stenseurd(3,1,2)=-sphi*cphi*Tzx/aa+sphi*cphi*Tzxda
         Stenseurd(3,2,2)=sphi*sphi*Tzx/aa+cphi*cphi*Tzxda
         Stenseurd(3,3,2)=cphi*Tzzda
         
c     derivee suivant l axe Z
         aa=a
         if (a.le.1.d-30) then
            cphi=0.d0
            sphi=0.d0
            aa=1.d300
         else
            sphi=xxp/a
            cphi=yyp/a
         endif
         s2phi=2.d0*sphi*cphi 
         c2phi=cphi*cphi-sphi*sphi
         Stenseurd(1,1,3)=Txxdz+c2phi*Txydz
         Stenseurd(1,2,3)=-s2phi*Txydz
         Stenseurd(1,3,3)=sphi*Txzdz
         Stenseurd(2,1,3)=-s2phi*Txydz
         Stenseurd(2,2,3)=Txxdz-c2phi*Txydz
         Stenseurd(2,3,3)=cphi*Txzdz
         Stenseurd(3,1,3)=sphi*Tzxdz
         Stenseurd(3,2,3)=cphi*Tzxdz
         Stenseurd(3,3,3)=Tzzdz        
      endif
      
      end
