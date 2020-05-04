      implicit none
      integer i,j,k,nx,ny,nz,ndipole,ntest
      double precision rayon,ray,x,y,z,xc,yc,zc,aretecube
      double complex eps,eps0


      nx=20
      ny=20
      nz=20
      rayon=100.d0
      aretecube=2.d0*rayon/dble(nx)
      xc=0.d0
      yc=0.d0
      zc=0.d0

      eps0=(1.d0,0.d0)
      eps=(2.25d0,0.d0)

      
      open(15,file='arbitraryobject')

      write(15,*)  nx,ny,nz
      write(15,*) aretecube
      ndipole=0
      do i=1,nz
         do j=1,ny
            do k=1,nx
               ndipole=ndipole+1
               x=-rayon+aretecube*(dble(k)-0.5d0)+xc
               y=-rayon+aretecube*(dble(j)-0.5d0)+yc
               z=-rayon+aretecube*(dble(i)-0.5d0)+zc
               write(15,*) x,y,z
            enddo
         enddo
      enddo
      write(*,*) ' ndipole',ndipole
      do i=1,nz
         do j=1,ny
            do k=1,nx
               x=-rayon+aretecube*(dble(k)-0.5d0)+xc
               y=-rayon+aretecube*(dble(j)-0.5d0)+yc
               z=-rayon+aretecube*(dble(i)-0.5d0)+zc
            
               ray=dsqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))
               ntest=0
               if (ray.le.rayon) then
                  write(15,*) eps
               else
                  write(15,*) eps0
               endif
            enddo
         enddo
      enddo
      
      end
      
