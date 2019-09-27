      subroutine  coeffmie(rad,k0,eps,eps0,coeffmiee)
      implicit none
      integer i,j,k
      double precision rad,k0,eps0
      double complex eps,x,y,index,icomp
      double complex nume,deno,zzz,zzzd,zdsin,zdcos,coeffmiee

      icomp=(0.d0,1.d0)
      x=k0*rad*(1.d0,0.d0)*dsqrt(eps0)
      index=cdsqrt(eps)/dsqrt(eps0)

      if (dimag(eps).ne.0.d0) then
         if (dimag(index).lt.0.d0) index=-index
      endif   

      y=index*x
      coeffmiee=0.d0
      nume=index*(zdsin(y)/y-zdcos(y))
     *     *((zdcos(x)*x-zdsin(x))/x/x+zdsin(x))

      nume=nume-(zdsin(x)/x-zdcos(x))*
     *     ((zdcos(y)*y-zdsin(y))/y/y+zdsin(y))

      zzz=zdsin(x)/x-zdcos(x)-icomp*(zdcos(x)/x+zdsin(x))
      zzzd=(zdcos(x)*x-zdsin(x))/x/x+zdsin(x)-
     *     icomp*((-zdsin(x)*x-zdcos(x))/x/x+zdcos(x))

      deno=index*(zdsin(y)/y-zdcos(y))*zzzd-
     *     zzz*((zdcos(y)*y-zdsin(y))/y/y+zdsin(y))

      coeffmiee=nume/deno

c     transformation du mie coefficient en polarisabilite
      coeffmiee=coeffmiee*3.d0/2.d0*icomp*rad*rad*rad/x/x/x

      end
      double complex function zdcos(u)
      implicit none
      double complex icomp,u
      icomp=(0.d0,1.d0)
      zdcos=(cdexp(icomp*u)+cdexp(-icomp*u))/2.d0
      return
      end
      double complex function zdsin(u)
      implicit none
      double complex icomp,u
      icomp=(0.d0,1.d0)
      zdsin=(cdexp(icomp*u)-cdexp(-icomp*u))/2.d0/icomp
      return
      end
