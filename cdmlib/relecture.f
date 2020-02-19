      subroutine relecture(aretecube,lambda,nrig,nfft2,beam,object,trope
     $     ,nnnr,tolinit,
c     cube, sphere (includes multiple)
     $     density, side, sidex, sidey, sidez, hauteur,
     $     numberobjet, rayonmulti, xgmulti, ygmulti, zgmulti,
     $     epsmulti, epsanimulti,lc,hc,ng,
c     ellipsoid+arbitrary
     $     demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,
     $     namefileobj,
c     planewavecircular.in / planewavelinear.in files
     $     theta, phi, pp, ss,thetam ,phim, ppm, ssm,E0m,nbinc, P0, w0,
     $     xgaus, ygaus, zgaus,namefileinc,numberobjetmax,polarizability
     $     ,nquad,filereread,nlecture,nx,ny,nz,nstop,infostr)
      implicit none
      integer ntmp,ierror,test,ii,jj,k,nlecture,nnnr,numberobjet
     $     ,numberobjetmax,nstop,long,long1,nquad,ng,nrig,nfft2,nx,ny,nz
      double precision tmp,lambda,tolinit, side, sidex, sidey, sidez,
     $     hauteur,demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj
     $     ,theta, phi, pp, ss, P0, w0, xgaus, ygaus, zgaus,xg,yg,zg
     $     ,rayon,xgmulti(numberobjetmax),ygmulti(numberobjetmax)
     $     ,zgmulti(numberobjetmax),rayonmulti(numberobjetmax),lc,hc
     $     ,density,aretecube
      integer nbinc
      double precision thetam(10), phim(10), ppm(10), ssm(10)
      double complex E0m(10)
      character (64) infostr,object,objecttmp,beam,beamtmp,namefileobj
     $     ,namefileinc,filereread,filereread1,filereread2
      character (3) trope,chatmp
      character (5) file1
      character (4) file2
      character (2) polarizability,polatmp
      double complex eps,epsani(3,3),epsmulti(numberobjetmax)
     $     ,epsanimulti(3,3,numberobjetmax),ctmp
      
      file1='.data'
      file2='.bin'
      xg = xgmulti(1)
      yg = ygmulti(1)
      zg = zgmulti(1)
      rayon = rayonmulti(1)
      eps = epsmulti(1)
      epsani = epsanimulti(:,:,1)
      long = len( trim(filereread  ) )
      long1 = len( trim( file1 ) )
      filereread1=filereread(1:long)//file1(1:long1)
      long1 = len( trim( file2 ) )
      filereread2=filereread(1:long)//file2(1:long1)

      if (long.eq.0) then
         infostr='No file given for read local field from file'
         nstop=1
         return
      endif
      
c     write(*,*) 'file',filereread1,filereread2
      open(1000,file=filereread1,status='old',iostat=ierror)
c     write(*,*) 'erreur',ierror
c     le fichier n'a pas encore ete cree et on le cree
      if (ierror.ne.0) then
         open(1000,file=filereread1,status='new')
         open(1001,file=filereread2,status='new',form='unformatted')
         write(1000,*) 'lambda  =',lambda
         write(1000,*) 'power   =',P0
         write(1000,*) 'waist   =',w0
         write(1000,*) 'beam    =',Beam
         write(1000,*) 'nrig    =',nrig
         write(1001) lambda
         write(1001) P0
         write(1001) w0
         write(1001) Beam
         write(1001) nrig
         
         if (Beam(1:11).eq.'pwavelinear') then
            write(1000,*) 'pola p  =',pp
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1001) pp
            write(1001) ss
            write(1001) theta
            write(1001) phi            
         elseif (Beam(1:13).eq.'pwavecircular') then
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1001) ss
            write(1001) theta
            write(1001) phi
         elseif (Beam(1:15).eq.'wavelinearmulti') then
            write(1000,*) 'number of incidence',nbinc
            do k=1,nbinc
               write(1000,*) 'pola p  =',ppm(k)
               write(1000,*) 'pola s  =',ssm(k)
               write(1000,*) 'theta   =',thetam(k)
               write(1000,*) 'phi     =',phim(k)
               write(1000,*) 'E0      =',E0m(k)
            enddo            
            write(1001) nbinc
            do k=1,nbinc
               write(1001) ppm(k)
               write(1001) ssm(k)
               write(1001) thetam(k)
               write(1001) phim(k)
               write(1001) E0m(k)
            enddo
         elseif (Beam(1:11).eq.'gwavelinear') then
            write(1000,*) 'pola p  =',pp
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1001) pp
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
         elseif (Beam(1:13).eq.'gwavecircular') then
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
         elseif (beam(1:15).eq.'gparawavelinear') then
            write(1000,*) 'pola p  =',pp
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1001) pp
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
         elseif (beam(1:17).eq.'gparawavecircular') then
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
         elseif (beam(1:14).eq.'gfftwavelinear') then
            write(1000,*) 'pola p  =',pp
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1000,*) 'nfft2   =',nfft2
            write(1001) pp
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
            write(1001) nfft2
         elseif (beam(1:16).eq.'gfftwavecircular') then
            write(1000,*) 'pola s  =',ss
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1000,*) 'nfft2   =',nfft2
            write(1001) ss
            write(1001) theta
            write(1001) phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
            write(1001) nfft2            
         elseif (beam(1:7).eq.'antenna') then         
            write(1000,*) 'pos x   =',xgaus
            write(1000,*) 'pos y   =',ygaus
            write(1000,*) 'pos z   =',zgaus
            write(1000,*) 'theta   =',theta
            write(1000,*) 'phi     =',phi
            write(1001) xgaus
            write(1001) ygaus
            write(1001) zgaus
            write(1001) theta
            write(1001) phi
         endif

         write(1000,*) 'discretization=',nnnr
         write(1000,*) 'object  =',Object
         write(1001) nnnr
         write(1001) Object
         if (object(1:6).eq.'sphere') then
            write(1000,*) 'rayon   =',rayon
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) rayon
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope

            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif

         elseif (object(1:13).eq.'randomsphere1') then
            write(1000,*) 'rayon   =',rayon
            write(1000,*) 'seed    =',ng
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'sidex   =',sidex
            write(1000,*) 'sidey   =',sidey
            write(1000,*) 'sidez   =',sidez
            write(1000,*) 'density =',density

            
            write(1001) rayon
            write(1001) ng
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) sidex
            write(1001) sidey
            write(1001) sidez
            write(1001) density

         elseif (object(1:13).eq.'randomsphere2') then
            write(1000,*) 'rayon   =',rayon
            write(1000,*) 'seed    =',ng
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'nx      =',sidex
            write(1000,*) 'ny      =',sidey
            write(1000,*) 'nz      =',sidez
            write(1000,*) 'meshsize=',aretecube
            write(1000,*) 'density =',density
            
            write(1001) rayon
            write(1001) ng
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) nx
            write(1001) ny
            write(1001) ny
            write(1001) aretecube
            write(1001) density
                        
         elseif (object(1:12).eq.'inhomosphere') then

            write(1000,*) 'rayon   =',rayon
            write(1000,*) 'lc      =',lc
            write(1000,*) 'st. de. =',hc
            write(1000,*) 'seed    =',ng
            write(1001) rayon
            write(1001) lc
            write(1001) hc
            write(1001) ng

         elseif (object(1:13).eq.'inhomocuboid1') then

            write(1000,*) 'sidex   =',sidex
            write(1000,*) 'sidey   =',sidey
            write(1000,*) 'sidez   =',sidez
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'lc      =',lc
            write(1000,*) 'st. de. =',hc
            write(1000,*) 'seed    =',ng

            write(1001) sidex
            write(1001) sidey
            write(1001) sidez
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) lc
            write(1001) hc
            write(1001) ng
         elseif (object(1:13).eq.'inhomocuboid2') then

            write(1000,*) 'nx      =',nx
            write(1000,*) 'ny      =',ny
            write(1000,*) 'nz      =',ny
            write(1000,*) 'meshsize=',aretecube
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'lc      =',lc
            write(1000,*) 'st. de. =',hc
            write(1000,*) 'seed    =',ng

            write(1001) nx
            write(1001) ny
            write(1001) ny
            write(1001) aretecube
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) lc
            write(1001) hc
            write(1001) ng
         elseif (object(1:4).eq.'cube') then
            write(1000,*) 'side    =',side
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) side
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope
            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif          
         elseif(object(1:7).eq.'cuboid1') then
            write(1000,*) 'sidex   =',sidex
            write(1000,*) 'sidey   =',sidey
            write(1000,*) 'sidez   =',sidez
            write(1000,*) 'psi     =',psiobj
            write(1000,*) 'theta   =',thetaobj
            write(1000,*) 'phi     =',phiobj
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) sidex
            write(1001) sidey
            write(1001) sidez
            write(1001) psiobj
            write(1001) thetaobj
            write(1001) phiobj
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope
            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif
         elseif(object(1:7).eq.'cuboid2') then
            write(1000,*) 'nx      =',nx
            write(1000,*) 'ny      =',ny
            write(1000,*) 'nz      =',ny
            write(1000,*) 'meshsize=',aretecube
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) nx
            write(1001) ny
            write(1001) ny
            write(1001) aretecube
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope
            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif  
            
         elseif(object(1:9).eq.'ellipsoid') then
            write(1000,*) 'demiaxex=',demiaxea
            write(1000,*) 'demiaxey=',demiaxeb
            write(1000,*) 'demiaxez=',demiaxec
            write(1000,*) 'psi     =',psiobj
            write(1000,*) 'theta   =',thetaobj
            write(1000,*) 'phi     =',phiobj
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) demiaxea
            write(1001) demiaxeb
            write(1001) demiaxec
            write(1001) psiobj
            write(1001) thetaobj
            write(1001) phiobj
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope
            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif  
         elseif(object(1:8).eq.'nspheres') then
            write(1000,*) 'num. obj=',numberobjet
            write(1000,*) 'iso-ani =',trope
            write(1001) numberobjet
            write(1001) trope
            if (trope.eq.'iso') then
               do k=1,numberobjet              
                  write(1000,*) 'pos x   =',xgmulti(k)
                  write(1000,*) 'pos y   =',ygmulti(k)
                  write(1000,*) 'pos z   =',zgmulti(k)
                  write(1000,*) 'rayon   =',rayonmulti(k)
                  write(1000,*) 'epsilon =',epsmulti(k)
                  write(1001) xgmulti(k)
                  write(1001) ygmulti(k)
                  write(1001) zgmulti(k)
                  write(1001) rayonmulti(k)
                  write(1001) epsmulti(k)
               enddo               
            else
               do k=1,numberobjet              
                  write(1000,*) 'pos x   =',xgmulti(k)
                  write(1000,*) 'pos y   =',ygmulti(k)
                  write(1000,*) 'pos z   =',zgmulti(k)
                  write(1000,*) 'rayon   =',rayonmulti(k)
                  write(1001) xgmulti(k)
                  write(1001) ygmulti(k)
                  write(1001) zgmulti(k)
                  write(1001) rayonmulti(k)
                  do ii=1,3
                     do jj=1,3
                        write(1000,*) 'epsilon =',epsanimulti(ii,jj,k)
                        write(1001) epsanimulti(ii,jj,k)
                     enddo
                  enddo
               enddo
            endif            
         elseif(object(1:8).eq.'cylinder') then
            write(1000,*) 'rayon   =',rayon
            write(1000,*) 'height  =',hauteur
            write(1000,*) 'psi     =',psiobj
            write(1000,*) 'theta   =',thetaobj
            write(1000,*) 'phi     =',phiobj
            write(1000,*) 'pos x   =',xg
            write(1000,*) 'pos y   =',yg
            write(1000,*) 'pos z   =',zg
            write(1000,*) 'iso-ani =',trope
            write(1001) rayon
            write(1001) hauteur
            write(1001) psiobj
            write(1001) thetaobj
            write(1001) phiobj
            write(1001) xg
            write(1001) yg
            write(1001) zg
            write(1001) trope
            if (trope.eq.'iso') then
               write(1000,*) 'epsilon =',eps
               write(1001) eps
            else
               do ii=1,3
                  do jj=1,3
                     write(1000,*) 'epsilon =',epsani(ii,jj)
                     write(1001) epsani(ii,jj)
                  enddo
               enddo
            endif  
         elseif(object(1:16).eq.'concentricsphere') then
            write(1000,*) 'num. obj=',numberobjet
            write(1000,*) 'iso-ani =',trope
            write(1001) numberobjet
            write(1001) trope
            if (trope.eq.'iso') then
               do k=1,numberobjet              
                  write(1000,*) 'pos x   =',xgmulti(k)
                  write(1000,*) 'pos y   =',ygmulti(k)
                  write(1000,*) 'pos z   =',zgmulti(k)
                  write(1000,*) 'rayon   =',rayonmulti(k)
                  write(1000,*) 'epsilon =',epsmulti(k)
                  write(1001) xgmulti(k)
                  write(1001) ygmulti(k)
                  write(1001) zgmulti(k)
                  write(1001) rayonmulti(k)
                  write(1001) epsmulti(k)
               enddo               
            else
               do k=1,numberobjet              
                  write(1000,*) 'pos x   =',xgmulti(k)
                  write(1000,*) 'pos y   =',ygmulti(k)
                  write(1000,*) 'pos z   =',zgmulti(k)
                  write(1000,*) 'rayon   =',rayonmulti(k)
                  write(1001) xgmulti(k)
                  write(1001) ygmulti(k)
                  write(1001) zgmulti(k)
                  write(1001) rayonmulti(k)
                  do ii=1,3
                     do jj=1,3
                        write(1000,*) 'epsilon =',epsanimulti(ii,jj,k)
                        write(1001) epsanimulti(ii,jj,k)
                     enddo
                  enddo
               enddo
            endif            
         elseif (object(1:9).eq.'arbitrary') then
            long = len( trim(namefileobj ) )
            write(1000,*) 'file obj =',namefileobj(1:long)
            write(1001) long
            write(1001) namefileobj(1:long)
         endif
c     type de polarisabilité
         write(1000,*) 'polarizability   =',polarizability
         write(1001) polarizability
c     type de la fonction Green
         write(1000,*) 'quadrature   =',nquad
         write(1001) nquad

         nlecture=0
         close(1000)
         close(1001)
         return
      else
c*******************************************************
c     ICI ON ECRIT LE FICHIER
c*******************************************************
         open(1000,file=filereread2,status='old',form='unformatted')
         read(1000) tmp
         call comparaisonreel(tmp,lambda,test)
         if (test.ne.0) goto 1000
         read(1000) tmp
         call comparaisonreel(tmp,P0,test)
         if (test.ne.0) goto 1000
         read(1000) tmp
         call comparaisonreel(tmp,w0,test)
         if (test.ne.0) goto 1000
         read(1000) Beamtmp
         if (.not.(lle(Beamtmp,Beam))) goto 1000
         read(1000) ntmp
         if (ntmp.ne.nrig) goto 1000
         
         if (Beam(1:11).eq.'pwavelinear') then
            read(1000) tmp
            call comparaisonreel(tmp,pp,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
         elseif (Beam(1:13).eq.'pwavecircular') then           
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
         elseif (Beam(1:15).eq.'wavelinearmulti') then           
            read(1000) ntmp
            if (ntmp.ne.nbinc) goto 1000            
            do k=1,nbinc
               read(1000) tmp
               call comparaisonreel(tmp,ppm(k),test)
               read(1000) tmp
               call comparaisonreel(tmp,ssm(k),test)
               if (test.ne.0) goto 1000
               read(1000) tmp
               call comparaisonreel(tmp,thetam(k),test)
               if (test.ne.0) goto 1000
               read(1000) tmp
               call comparaisonreel(tmp,phim(k),test)
               if (test.ne.0) goto 1000
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,E0m(k),test)
               if (test.ne.0) goto 1000
            enddo            
         elseif (Beam(1:11).eq.'gwavelinear') then
            read(1000) tmp
            call comparaisonreel(tmp,pp,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000           
         elseif (Beam(1:13).eq.'gwavecircular') then
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000
         elseif (Beam(1:15).eq.'gparawavelinear') then
            read(1000) tmp
            call comparaisonreel(tmp,pp,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000
         elseif (Beam(1:17).eq.'gparawavecircular') then
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000

         elseif (Beam(1:14).eq.'gfftwavelinear') then
            read(1000) tmp
            call comparaisonreel(tmp,pp,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nfft2) goto 1000
         elseif (Beam(1:16).eq.'gfftwavecircular') then
            read(1000) tmp
            call comparaisonreel(tmp,ss,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nfft2) goto 1000
         elseif (Beam(1:7).eq.'antena') then
            read(1000) tmp
            call comparaisonreel(tmp,xgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,ygaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zgaus,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,theta,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phi,test)
            if (test.ne.0) goto 1000        
         endif

         read(1000) ntmp
c     write(*,*) '8',ntmp,nnnr
         if (ntmp.ne.nnnr) goto 1000
         read(1000) Objecttmp
         if (Objecttmp(1:16).ne.Object(1:16)) goto 1000
         

         if (object(1:6).eq.'sphere') then
            read(1000) tmp
            call comparaisonreel(tmp,rayon,test)
c     write(*,*) '9',tmp,rayon,test
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
c     write(*,*) '10',tmp,xg,test
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
c     write(*,*) '11',tmp,yg,test
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
c     write(*,*) '11',tmp,zg,test
            if (test.ne.0) goto 1000
            read(1000) chatmp
c     write(*,*) '12',chatmp,trope
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
c     write(*,*) '13',ctmp,eps
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif

         elseif (object(1:13).eq.'randomsphere1') then
            
            read(1000) tmp
            call comparaisonreel(tmp,rayon,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ng) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) yg
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) zg
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidex,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidey,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidez,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,density,test)
            if (test.ne.0) goto 1000

         elseif (object(1:13).eq.'randomsphere2') then
            
            read(1000) tmp
            call comparaisonreel(tmp,rayon,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ng) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) yg
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) zg
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nx) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ny) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nz) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,aretecube,test)
            if (test.ne.0) goto 1000
            call comparaisonreel(tmp,density,test)
            if (test.ne.0) goto 1000
                        
         elseif (object(1:12).eq.'inhomosphere') then
            read(1000) tmp
            call comparaisonreel(tmp,rayon,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,lc,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,hc,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ng) goto 1000

         elseif (object(1:13).eq.'inhomocuboid1') then

            read(1000) tmp
            call comparaisonreel(tmp,sidex,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidey,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidez,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000          
            read(1000) tmp
            call comparaisonreel(tmp,lc,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,hc,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ng) goto 1000
            

         elseif (object(1:13).eq.'inhomocuboid2') then
            read(1000) ntmp
            if (ntmp.ne.nx) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ny) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nz) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,aretecube,test)
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000          
            read(1000) tmp
            call comparaisonreel(tmp,lc,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,hc,test)
            if (test.ne.0) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ng) goto 1000
            

         elseif (object(1:4).eq.'cube') then
            read(1000) tmp
            call comparaisonreel(tmp,side,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) yg
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) zg
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000            
            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif          
         elseif(object(1:7).eq.'cuboid1') then
            read(1000) tmp
            call comparaisonreel(tmp,sidex,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidey,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,sidez,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,psiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,thetaobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000

            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif


         elseif(object(1:7).eq.'cuboid2') then
            read(1000) ntmp
            if (ntmp.ne.nx) goto 1000
            read(1000) ntmp
            if (ntmp.ne.ny) goto 1000
            read(1000) ntmp
            if (ntmp.ne.nz) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,aretecube,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000

            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif     

            
         elseif(object(1:9).eq.'ellipsoid') then
            read(1000) tmp
            call comparaisonreel(tmp,demiaxea,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,demiaxeb,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,demiaxec,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,psiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,thetaobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000
            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif  
         elseif(object(1:8).eq.'nspheres') then
            read(1000) ntmp
            if (ntmp.ne.numberobjet) goto 1000
            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000
            if (trope.eq.'iso') then
               do k=1,numberobjet              
                  read(1000) tmp
                  call comparaisonreel(tmp,xgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,ygmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,zgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,rayonmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) ctmp
                  call comparaisoncomplexe(ctmp,eps,test)
                  if (test.ne.0) goto 1000
               enddo               
            else
               do k=1,numberobjet              
                  read(1000) tmp
                  call comparaisonreel(tmp,xgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,ygmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,zgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,rayonmulti(k),test)
                  if (test.ne.0) goto 1000
                  do ii=1,3
                     do jj=1,3
                        read(1000) ctmp
                        call comparaisoncomplexe(ctmp,epsanimulti(ii,jj
     $                       ,k),test)
                        if (test.ne.0) goto 1000
                     enddo
                  enddo
               enddo
            endif            
         elseif(object(1:8).eq.'cylinder') then
            read(1000) tmp
            call comparaisonreel(tmp,rayon,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,hauteur,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,psiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,thetaobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,phiobj,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,xg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp
            call comparaisonreel(tmp,yg,test)
            if (test.ne.0) goto 1000
            read(1000) tmp         
            call comparaisonreel(tmp,zg,test)
            if (test.ne.0) goto 1000
            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000

            if (trope.eq.'iso') then
               read(1000) ctmp
               call comparaisoncomplexe(ctmp,eps,test)
               if (test.ne.0) goto 1000
            else
               do ii=1,3
                  do jj=1,3
                     read(1000) ctmp
                     call comparaisoncomplexe(ctmp,epsani(ii,jj),test)
                     if (test.ne.0) goto 1000
                  enddo
               enddo
            endif  
         elseif(object(1:16).eq.'concentricsphere') then
            read(1000) ntmp
            if (ntmp.ne.numberobjet) goto 1000
            read(1000) chatmp
            if (chatmp.ne.trope) goto 1000
            if (trope.eq.'iso') then
               do k=1,numberobjet              
                  read(1000) tmp
                  call comparaisonreel(tmp,xgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,ygmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,zgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,rayonmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) ctmp
                  call comparaisoncomplexe(ctmp,eps,test)
                  if (test.ne.0) goto 1000
               enddo               
            else
               do k=1,numberobjet              
                  read(1000) tmp
                  call comparaisonreel(tmp,xgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,ygmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,zgmulti(k),test)
                  if (test.ne.0) goto 1000
                  read(1000) tmp
                  call comparaisonreel(tmp,rayonmulti(k),test)
                  if (test.ne.0) goto 1000
                  do ii=1,3
                     do jj=1,3
                        read(1000) ctmp
                        call comparaisoncomplexe(ctmp,epsanimulti(ii,jj
     $                       ,k),test)
                        if (test.ne.0) goto 1000
                     enddo
                  enddo
               enddo
            endif
         elseif (object(1:9).eq.'arbitrary') then
            long = len( trim(namefileobj ) )
            read(1000) ntmp
            if (ntmp.ne.long) goto 1000
            read(1000) objecttmp(1:ntmp)
            if (objecttmp(1:ntmp).ne.namefileobj(1:long)) goto 1000
         endif
         
         read(1000) polatmp
         if (polatmp.ne.polarizability) goto 1000
         read(1000) ntmp
         if (ntmp.ne.nquad) goto 1000

         close(1000)
         nlecture=1
         return
      endif

 1000 nstop=1
      infostr='Datafile do not match the data or file already exists'
      return

c     900  format(e22.17)
      end
