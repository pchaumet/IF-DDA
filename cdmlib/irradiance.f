c     calcul E0 en fonction de l'irradiance pour une onde plane
c     I=(ExB)/(2 mu_0)=E^2/(2 c mu_0): eps_0 mu_0=1/c^2
c     I=eps_0 c^2 E^2/(2 c)=eps_0 c E^2/2
c     E^2=2 I/(eps_0 c)
c     F=4 pi eps_0 E^2= 4 pi eps_0 2 I/(eps_0 c)=8 pi I/c
      subroutine irradiance(P0,w0,E0,irra)
      implicit none

      double precision c,eps0,pi,P0,surf,w0,irra
      double complex E0

      pi=dacos(-1.d0)
      c=299792458.d0
      eps0=1.d0/(c*c*4.d0*pi*1.d-7)

      surf=w0*w0*pi
      Irra=P0/surf
      E0=dsqrt(2.d0*Irra/c/eps0)*(1.d0,0.d0)

      end
