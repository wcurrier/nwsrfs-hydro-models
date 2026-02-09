C     ==============================================================================
C     DUAMEL - Unit Hydrograph Routing with Restart Capability
C     ==============================================================================
C     
C     MODIFICATIONS (January 2026):
C     - Added warm start capability via QPREV/USE_QPREV
C     - Added QPREV_OUT for restart state saving
C     - Added QPREV_TS for optional time series output (controlled by SAVE_TS)
C     - Added SAVE_TS flag to prevent memory issues when time series not needed
C
C     ==============================================================================
C
      SUBROUTINE DUAMEL(Q,UN1,UT,DT,N,MM,K,NTAU,QB,
     &                  QPREV,USE_QPREV,QPREV_OUT,QPREV_TS,SAVE_TS)
C
C     ARGUMENT DOCUMENTATION:
C     =================================================================
C     INPUTS:
C       Q         : TCI (total channel inflow) vector - unrouted flow
C       UN1       : Unit hydrograph shape parameter (gamma distribution)
C       UT        : Unit hydrograph scale parameter (gamma distribution) 
C       DT        : Timestep of the UH function (days or fractions thereof)
C       N         : sim_length + MM (includes UH tail)
C       MM        : Unit hydrograph length
C       K         : UH computation flag (0 = skip UH build)
C       NTAU      : Lag parameter
C       QPREV     : Previous M flow values from prior simulation (warm start)
C       USE_QPREV : Flag to enable warm start (1=warm, 0=cold)
C       SAVE_TS   : Flag to save time series (1=yes, 0=no)
C
C     OUTPUTS:
C       QB        : Routed flow vector
C       QPREV_OUT : Last M values of Q for next restart
C       QPREV_TS  : qprev time series (only used if SAVE_TS=1)
C     =================================================================

      IMPLICIT REAL (A-H,O-Z)
      INTEGER A,B,I,J,M,IOR
      INTEGER N,MM,K,NTAU,USE_QPREV,M_SAVE,SAVE_TS
      INTEGER NQ,QIDX
      
      REAL Q(*), QB(*), QPREV(*), QPREV_OUT(*), QPREV_TS(*)
      REAL U1(MM)
      REAL QVAL
C
      M  = MM
      M_SAVE = MM
      NQ = N - MM
C
C----- Build unit hydrograph (unchanged from original)
C
      IF (UN1 .LT. 0.0) THEN
        U1(1) = 1.0
        M = 1
        GO TO 60
      ELSE
        IF (K .EQ. 0) GOTO 60
      END IF

      SP = 0.0
      TOC = GF(UN1)
      TOC = LOG(TOC*UT)

      DO 10 I=1,M
        TOP = I*DT/UT
        TOR = (UN1-1.0)*LOG(TOP) - TOP - TOC
        U1(I) = 0.0
        IF(TOR.GT.-8.0) THEN
          U1(I) = EXP(TOR)
        ELSE
          IF (I .GT. 1) THEN
            M = I
            GO TO 20
          END IF
        END IF
        SP = SP + U1(I)
   10 CONTINUE
   20 CONTINUE

      IF (SP .EQ. 0.0) SP = 1.0E-5
      SP = 1.0/SP
      DO 30 I=1,M
        U1(I) = U1(I)*SP
   30 CONTINUE

   60 CONTINUE
C
C----- CONVOLUTION WITH RESTART CAPABILITY
C
      IOC = N + NTAU

      DO 100 I=1,IOC
        QB(I) = 0.0
        
        DO 90 J=1,M
          IOR = J
          QIDX = I - J + 1
          
C         Determine which Q value to use based on QIDX position
          IF (QIDX .GE. 1 .AND. QIDX .LE. NQ) THEN
C           QIDX in range [1, NQ] - use current Q array
            QVAL = Q(QIDX)
          ELSE IF (QIDX .LE. 0 .AND. USE_QPREV .EQ. 1) THEN
C           QIDX in range [-M+1, 0] and warm start enabled - use QPREV
            QIDX = M_SAVE + QIDX
            IF (QIDX .GE. 1 .AND. QIDX .LE. M_SAVE) THEN
              QVAL = QPREV(QIDX)
            ELSE
              QVAL = 0.0
            END IF
          ELSE
C           Cold start or out of range - use zero
            QVAL = 0.0
          END IF
          
          QB(I) = QB(I) + QVAL*U1(IOR)
   90   CONTINUE
   
C       Save qprev state at this timestep (ONLY if SAVE_TS=1)
        IF (SAVE_TS .EQ. 1 .AND. I .GE. 1 .AND. I .LE. NQ) THEN
          DO 95 J=1,M_SAVE
            QIDX = I - J + 1
            IF (QIDX .GE. 1 .AND. QIDX .LE. NQ) THEN
              QPREV_TS((I-1)*M_SAVE + J) = Q(QIDX)
            ELSE IF (QIDX .LE. 0 .AND. USE_QPREV .EQ. 1) THEN
              QIDX = M_SAVE + QIDX
              IF (QIDX .GE. 1 .AND. QIDX .LE. M_SAVE) THEN
                QPREV_TS((I-1)*M_SAVE + J) = QPREV(QIDX)
              ELSE
                QPREV_TS((I-1)*M_SAVE + J) = 0.0
              END IF
            ELSE
              QPREV_TS((I-1)*M_SAVE + J) = 0.0
            END IF
   95     CONTINUE
        END IF
        
  100 CONTINUE
C
C----- SAVE FINAL M VALUES FOR RESTART
C
      DO 200 I=1,M_SAVE
        IF (NQ - M_SAVE + I .GT. 0) THEN
          QPREV_OUT(I) = Q(NQ - M_SAVE + I)
        ELSE IF (USE_QPREV .EQ. 1) THEN
          QIDX = I - (M_SAVE - NQ)
          IF (QIDX .GE. 1 .AND. QIDX .LE. M_SAVE) THEN
            QPREV_OUT(I) = QPREV(QIDX)
          ELSE
            QPREV_OUT(I) = 0.0
          END IF
        ELSE
          QPREV_OUT(I) = 0.0
        END IF
  200 CONTINUE
C
      RETURN
      END
      
C
C=================================================================
C
      FUNCTION GF(Y)
      REAL Y, X, H, GF
      GF=0.0
      H=1
      X=Y
 38   IF(X.LE.0.)GO TO 39
      IF(X.EQ.2.)GO TO 42
      IF(X.GT.2.)GO TO 40
      H=H/X
      X=X+1
      GO TO 38
  40  IF(X.LE.3.)GO TO 44
      X=X-1
      H=H*X
      GO TO 38
  44  X=X-2
      H=(((((((.0016063118*X+0.0051589951)*X+0.0044511400)*X+.0721101567
     *)*X+.0821117404)*X+.4117741955)*X+.4227874605)*X+.9999999758)*H
 42   GF=H
  39  RETURN
      END
