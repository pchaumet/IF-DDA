      integer function comparaison(x,y,z,x0,y0,z0,lambda)
      implicit none
      double precision x,y,z,x0,y0,z0,lambda,erreur
      erreur=lambda/1000.d0
      comparaison=0
      if (x.ge.x0-erreur.and.x.le.x0+erreur.and.y.ge.y0
     $     -erreur.and.y.le.y0+erreur.and.z.ge.z0-erreur.and.z.le.z0
     $     +erreur) comparaison=1
      return
      end

      
