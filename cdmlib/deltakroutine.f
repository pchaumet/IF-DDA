      subroutine deltakroutine(kxinc,kyinc,deltakx,deltaky,numaper,i,j)
      implicit none
      integer i,j
      double precision kx,ky,kxinc,kyinc,deltakx,deltaky,numaper

      i=idnint(kxinc/deltakx)
      j=idnint(kyinc/deltaky)
      kx=deltakx*dble(i)
      ky=deltaky*dble(j)
      if (kx*kx+ky*ky.le.numaper*numaper) then
         return
      else
         i=idint(kxinc/deltakx)
         j=idint(kyinc/deltaky)
      endif

      end
