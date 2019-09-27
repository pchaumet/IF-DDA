      SUBROUTINE  CALLBHMIE(REFMED,CXEPS,RAD,WAVEL,MIECEXT,MIECABS
     $     ,MIECSCA,GSCA)
! Parameters:
      INTEGER  MXNANG
      PARAMETER(MXNANG=1000)
! Variables:
      INTEGER IREADEP,J,NAN,NANG,NANG0
      double precision AJ,ANG,DANG,GSCA,PI,POL
      double precision  MIECEXT,MIECABS,MIECSCA
      double precision QABS,QBACK,QEXT,QSCA,RAD,REFMED,S11,S12,S33,S34
     $     ,WAVEL,X
      double complex REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
!!$C***********************************************************************
!!$C Program to interactively call Bohren-Huffman Mie theory program
!!$C
!!$C CALLBHMIE will interactively prompt for:
!!$C 1. refractive index of surrounding medium
!!$C 2. either refractive index or dielectric constant of sphere
!!$C 3. radius of sphere
!!$C 4. wavelength (in vacuo)
!!$C 5. number of angles at which to calculate scattering intensities
!!$C
!!$C CALLBHMIE will return:
!!$C 1. Q_ext, Q_abs, Q_sca, g, Q_back
!!$C 2. If NANG>0, then will also return scattering matrix elements
!!$C    S_11, S_33, S_34, and POL
!!$C
!!$C Adapted by B.T.Draine, Princeton Univ. Obs.
!!$C***********************************************************************
      PI=dacos(-1.d0)
!!$      OPEN(UNIT=7,FILE='callbhmie.out',STATUS='UNKNOWN')
!!$      WRITE(*,*)' Enter (real) refractive index of surrounding medium'
!!$      READ(*,*)REFMED
!!$      WRITE(*,*)' Wish to enter refr.index or epsilon? (0 or 1)'
!!$      READ(*,*)IREADEP
!!$ 1000 IF(IREADEP.LE.0)THEN
!!$          WRITE(*,*)' Enter COMPLEX (KIND=KIND(0.d0)) :: refractive index of sphere'
!!$          READ(*,*)REFREL
!!$      ELSE
!!$          WRITE(*,*)' Enter COMPLEX (KIND=KIND(0.d0)) :: epsilon of sphere'
!!$          READ(*,*)CXEPS
         REFREL=cdSQRT(CXEPS)
!!$      ENDIF
         REFREL=REFREL/REFMED
!!$      WRITE(*,6012)REFREL
!!$      WRITE(7,6012)REFREL
!!$ 6012 FORMAT(' COMPLEX (KIND=KIND(0.d0)) :: refractive index=',1P2E10.3)
!!$      WRITE(*,*)' Enter radius'
!!$      READ(*,*)RAD
!!$      WRITE(*,*)' Enter wavelength'
!!$      READ(*,*)WAVEL
!!$      WRITE(*,*)' Enter NANG = number of angles between 0 and 90'
!!$      READ(*,*)NANG0

      NANG0=1

      IF(NANG0.GT.MXNANG)STOP'***Error: NANG > MXNANG'
      NANG=NANG0
      IF(NANG0.LT.2)NANG=2
      X=2.d0*PI*RAD*REFMED/WAVEL
!!$      WRITE(7,6013)RAD,WAVEL,X
!!$      WRITE(*,6013)RAD,WAVEL,X
!!$C**********
!!$C NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
!!$C Scattering matrix elements are calculated for 2*NANG-1 angles
!!$C including 0, 90, and 180 degrees.
!!$C**********
      IF(NANG.GT.1) DANG=0.5d0*PI/dble(NANG-1)
      CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
      QABS=QEXT-QSCA
      
      MIECEXT=pi*RAD*RAD*QEXT
      MIECABS=pi*RAD*RAD*QABS
      MIECSCA=pi*RAD*RAD*QSCA
      

!!$      WRITE(*,6065)QEXT,QSCA,QABS,GSCA,QBACK
!!$      WRITE(7,6065)QEXT,QSCA,QABS,GSCA,QBACK
!!$C POL=degree of polarization (for incident upolarized light)
      IF(NANG0.GT.1)THEN
          NAN=2*NANG-1
!!$          WRITE(*,6017)
!!$          WRITE(7,6017)
          DO J=1,NAN
             AJ=J
             S11=0.5d0*CdABS(S2(J))*CdABS(S2(J))
             S11=S11+0.5d0*CdABS(S1(J))*CdABS(S1(J))
             S12=0.5d0*CdABS(S2(J))*CdABS(S2(J))
             S12=S12-0.5d0*CdABS(S1(J))*CdABS(S1(J))
             POL=-S12/S11
             S33=dREAL(S2(J)*dCONJG(S1(J)))
             S34=dIMAG(S2(J)*dCONJG(S1(J)))
             ANG=DANG*(AJ-1.d0)*180.d0/PI
          END DO
      ENDIF

      END
c******************************************************************     
      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)

       IMPLICIT NONE

! Declare parameters:
! Note: important that MXNANG be consistent with dimension of S1 and S2
!       in calling routine!
      INTEGER  MXNANG,NMXX
      PARAMETER(MXNANG=1000,NMXX=15000)
! Arguments:
      INTEGER  NANG
      double precision GSCA,QBACK,QEXT,QSCA,X
      double COMPLEX  REFREL
      double COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)
! Local variables:
      INTEGER  J,JJ,N,NSTOP,NMX,NN
      double precision APSI,APSI1,CHI,CHI0,CHI1,DANG,FN,P,PII,RN,THETA
     $     ,XSTOP,YMOD
      double precision AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG)
     $     ,TAU(MXNANG)
      double precision PSI0,PSI1,PSI,DN,DX
      double COMPLEX  AN,AN1,BN,BN1,XI,XI1,Y
      double COMPLEX  D(NMXX)
!***********************************************************************
! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
!    to calculate scattering and absorption by a homogenous isotropic
!    sphere.
! Given:
!    X = 2*pi*a/lambda
!    REFREL = (COMPLEX (KIND=KIND(0.d0)) :: refr. index of sphere)/(REAL (KIND=KIND(0.d0)) :: index of medium)
!    NANG = number of angles between 0 and 90 degrees
!           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
!           if called with NANG<2, will set NANG=2 and will compute
!           scattering for theta=0,90,180.
! Returns:
!    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
!                                scatt. E perp. to scatt. plane)
!    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
!                                scatt. E parr. to scatt. plane)
!    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
!    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
!    QBACK = (dC_sca/domega)/pi*a**2
!          = backscattering efficiency
!    GSCA = <cos(theta)> for scattering
!
! Original program taken from Bohren and Huffman (1983), Appendix A
! Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
! in order to compute <cos(theta)>
! 91/05/07 (BTD): Modified to allow NANG=1
! 91/08/15 (BTD): Corrected error (failure to initialize P)
! 91/08/15 (BTD): Modified to enhance vectorizability.
! 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
! 91/08/15 (BTD): Changed definition of QBACK.
! 92/01/08 (BTD): Note that this version has been superceded by
!                 fully double precision version = bhmie.f which,
!                 unfortunately, is not standard f77.
!                 However, retain this in case standard f77 version
!                 is required for porting to some other system.
!***********************************************************************
!*** Safety checks
      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2
!*** Obtain pi:
      PII=dACOS(-1.d0)
      DX=X
      Y=X*REFREL
      YMOD=cdABS(Y)
!
!*** Series expansion terminated after NSTOP terms
!    Logarithmi! derivatives calculated from NMX on down
      XSTOP=X+4.d0*X**0.3333d0+2.d0
!*** Original code:
!      NMX=AMAX1(XSTOP,YMOD)+15
!      NSTOP=XSTOP
!*** Experimental code:
      NMX=1.d0*DMAX1(XSTOP,YMOD)+15
      NSTOP=1.d0*XSTOP
!
      IF(NMX.GT.NMXX)THEN
          WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
          STOP
      ENDIF
!*** Require NANG.GE.1 in order to calculate scattering intensities
      DANG=0.d0
      IF(NANG.GT.1)DANG=0.5d0*PII/dble(NANG-1)
      DO J=1,NANG
         THETA=dble(J-1)*DANG
         AMU(J)=dCOS(THETA)
      END DO

      DO J=1,NANG
         PI0(J)=0.d0
         PI1(J)=1.d0
      END DO

      NN=2*NANG-1
      DO J=1,NN
         S1(J)=(0.d0,0.d0)
         S2(J)=(0.d0,0.d0)
      END DO
!
!*** Logarithmi! derivative D(J) calculated by downward recurrence
!    beginning with initial value (0.,0.) at J=NMX
!
      D(NMX)=(0.d0,0.d0)
      NN=NMX-1
      DO N=1,NN
         RN=NMX-N+1
         D(NMX-N)=(RN/Y)-(1.d0/(D(NMX-N+1)+RN/Y))
      END DO
!
!*** Riccati-Bessel functions with REAL argument X
!    calculated by upward recurrence
!
      PSI0=dCOS(DX)
      PSI1=dSIN(DX)
      CHI0=-dSIN(X)
      CHI1=dCOS(X)
! APSI0 never used, so this line removed from program:
!      APSI0=PSI0
      APSI1=PSI1
! XI0 never used, so this line removed from program:
!      XI0=CMPLX(APSI0,-CHI0)
      XI1=dCMPLX(APSI1,-CHI1)
      QSCA=0.d0
      GSCA=0.d0
      P=-1.d0
      BIG_LOOP : DO N=1,NSTOP

         DN=N
         RN=N
         FN=(2.d0*RN+1.d0)/(RN*(RN+1.d0))
         PSI=(2.d0*DN-1.d0)*PSI1/DX-PSI0
         APSI=PSI
         CHI=(2.d0*RN-1.d0)*CHI1/X-CHI0
          XI=dCMPLX(APSI,-CHI)
          !
          !*** Store previous values of AN and BN for use
          !    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
             AN1=AN
             BN1=BN
          ENDIF
          !
          !*** Compute AN and BN:
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
          !
          !*** Augment sums for Qsca and g=<cos(theta)>
          QSCA=QSCA+(2.d0*RN+1.d0)*(CdABS(AN)**2+CdABS(BN)**2)
          GSCA=GSCA+((2.d0*RN+1.d0)/(RN*(RN+1.d0)))*(dREAL(AN)*dREAL(BN)
     $         +dIMAG(AN)*dIMAG(BN))
          IF(N.GT.1)THEN
             GSCA=GSCA+((RN-1.d0)*(RN+1.d0)/RN)*(dREAL(AN1)*dREAL(AN)
     $            +dIMAG(AN1)*dIMAG(AN)+dREAL(BN1)*dREAL(BN)+dIMAG(BN1)
     $            *dIMAG(BN))
          ENDIF
          !
          !*** Now calculate scattering intensity pattern
          !    First do angles from 0 to 90
          DO J=1,NANG
             JJ=2*NANG-J
             PI(J)=PI1(J)
             TAU(J)=RN*AMU(J)*PI(J)-(RN+1.d0)*PI0(J)
             S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
             S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
          END DO
          !
          !*** Now do angles greater than 90 using PI and TAU from
          !    angles less than 90.
          !    P=1 for N=1,3,...; P=-1 for N=2,4,...
          P=-P
          DO J=1,NANG-1
             JJ=2*NANG-J
             S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
             S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
          END DO
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=dCMPLX(APSI1,-CHI1)
          !
          !*** Compute pi_n for next value of n
          !    For each angle J, compute pi_n+1
          !    from PI = pi_n , PI0 = pi_n-1
          DO J=1,NANG
             PI1(J)=((2.d0*RN+1.d0)*AMU(J)*PI(J)-(RN+1.d0)*PI0(J))/RN
             PI0(J)=PI(J)
          END DO

      END DO BIG_LOOP

!
!*** Have summed sufficient terms.
!    Now compute QSCA,QEXT,QBACK,and GSCA
      GSCA=2.d0*GSCA/QSCA
      QSCA=(2.d0/(X*X))*QSCA
      QEXT=(4.d0/(X*X))*dREAL(S1(1))
      QBACK=CdABS(S1(2*NANG-1))*CdABS(S1(2*NANG-1))/(PII*X*X)
      RETURN

      END



