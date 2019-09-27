      subroutine derivechamp(nx,ny,nz,nmax,aretecube,FF,test,FFder)
      implicit none
      integer nx,ny,nz,nmax,test,i,j,k,nxy,ii
      double precision aretecube,dx,dy,dz
      double complex FF(3*nmax),FFder(3*nmax)
      double complex FFaxe(1000),FFaxeder(1000)

      nxy=nx*ny
      dx=aretecube
      dy=aretecube
      dz=aretecube

      
      if (test.eq.1) then
c     calcul derivee suivant x
         do i=1,nz
            do j=1,ny
c     calcul de Exdx
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(k)=FF(ii+1)
               enddo
               call deriveFF(FFaxe,FFaxeder,nx,dx)
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+1)=FFaxeder(k)
               enddo
c     calcul de Eydx
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(k)=FF(ii+2)
               enddo
               call deriveFF(FFaxe,FFaxeder,nx,dx)
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+2)=FFaxeder(k)
               enddo
c     calcul de Ezdx
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(k)=FF(ii+3)
               enddo
               call deriveFF(FFaxe,FFaxeder,nx,dx)
               do k=1,nx
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+3)=FFaxeder(k)
               enddo
            enddo
         enddo
      elseif (test.eq.2) then
c     calcul des derivees suivant y
         do i=1,nz
            do k=1,nx
               do j=1,ny          
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(j)=FF(ii+1)
               enddo
               call deriveFF(FFaxe,FFaxeder,ny,dy)
               do j=1,ny
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+1)=FFaxeder(j)
               enddo

               do j=1,ny          
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(j)=FF(ii+2)
               enddo
               call deriveFF(FFaxe,FFaxeder,ny,dy)
               do j=1,ny
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+2)=FFaxeder(j)
               enddo

               do j=1,ny          
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(j)=FF(ii+3)
               enddo
               call deriveFF(FFaxe,FFaxeder,ny,dy)
               do j=1,ny
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+3)=FFaxeder(j)
               enddo

            enddo
         enddo

c     calcul des derivees suivant z
      elseif (test.eq.3) then
         do k=1,nx
            do j=1,ny 
               do i=1,nz         
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(i)=FF(ii+1)
               enddo
               call deriveFF(FFaxe,FFaxeder,nz,dz)
               do i=1,nz
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+1)=FFaxeder(i)
               enddo

               do i=1,nz         
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(i)=FF(ii+2)
               enddo
               call deriveFF(FFaxe,FFaxeder,nz,dz)
               do i=1,nz
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+2)=FFaxeder(i)
               enddo


               do i=1,nz         
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFaxe(i)=FF(ii+3)
               enddo
               call deriveFF(FFaxe,FFaxeder,nz,dz)
               do i=1,nz
                  ii=3*(k+nx*(j-1)+nxy*(i-1)-1)
                  FFder(ii+3)=FFaxeder(i)
               enddo

            enddo
         enddo
      else
         write(*,*) 'pb with test for derivative'
      endif

      end
