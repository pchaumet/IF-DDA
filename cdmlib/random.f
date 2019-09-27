C************************************************************
C************************************************************
C                                                            
C                                                            
C       0 , 1   - U N I F O R M  DISTRIBUTION                
C                                                                 
C                                                                 
C*******************************************************
C                                                                 
C     FOR DETAILS SEE:                                            
C                                                                   
C               AHRENS, J.H., DIETER, U. AND GRUBE, A.             
C               PSEUDO-RANDOM NUMBERS:  A NEW PROPOSAL        
C                     FOR THE CHOICE OF MULTIPLICATORS           
C               COMPUTING, 6 (1970), 121 - 138                  
C                                                               
C*************************************************************
C                                                                  
      double precision FUNCTION SUNIF(IR)
      DOUBLE PRECISION R,FACTOR,TWO28
      SAVE R
C
C     FACTOR - INTEGER OF THE FORM 8*K+5 AS CLOSE AS POSSIBLE
C              TO  2**26 * (SQRT(5)-1)/2     (GOLDEN SECTION)
C     TWO28  = 2**28  (I.E. 28 SIGNIFICANT BITS FOR DEVIATES)
C
      DATA FACTOR /41475557.0D0/, TWO28 /268435456.0D0/
C
C     RETURNS SAMPLE U FROM THE  0,1 -UNIFORM DISTRIBUTION
C     BY A MULTIPLICATIVE CONGRUENTIAL GENERATOR OF THE FORM
C        R := R * FACTOR (MOD 1) .
C     IN THE FIRST CALL R IS INITIALIZED TO
C        R := IR / 2**28 ,
C     WHERE IR MUST BE OF THE FORM  IR = 4*K+1.
C     THEN R ASSUMES ALL VALUES  0 < (4*K+1)/2**28 < 1 DURING
C     A FULL PERIOD 2**26 OF SUNIF.
C     THE PARAMETER IR IS USED ONLY IN THE FIRST CALL FOR
C     INITIALIZATION OF SUNIF. THEREAFTER (WHEN NEGATIVE)
C     IR BECOMES A DUMMY VARIABLE.
C
      IF (IR .GE. 0) GO TO 1
C
C     STANDARD CASE:  SAMPLING
C
      R=DMOD(R*FACTOR,1.0D0)
      SUNIF=SNGL(R)
      RETURN
C
C     FIRST CALL: INITIALIZATION
C
1     R=DBLE(FLOAT(IR))/TWO28
      R=DMOD(R*FACTOR,1.0D0)
      SUNIF=SNGL(R)
      IR=-1
      RETURN
      END
c*********************************************************
c*********************************************************
      SUBROUTINE rannor (x,m,s)

c     Generates a pseudo-random normal deviate (x) with mean (m) and
c     standard deviation (s) using a sum of uniform deviates
c     tactic. Calls system subroutine RANDOM, a psuedo-random uniform
c     U{0,1} generator with argument (uniform) that contains the uniform
c     variate. RANDOM must be initialized with a proper seed prior to
c     the first call of rannor.

c       input: m,s - mean and standard deviation of normal 
c       output: x - pseudo random normal deviate

      implicit none
      double precision x,m,s,ran,SUNIF
      INTEGER i,IR

      ir=-1

      x = 0.d0
      DO i = 1,12                   
         ran=SUNIF(IR)
         x = ran + x 
      ENDDO
      x = s * ( x - 6.d0 ) + m
      END 
c********************************************************
