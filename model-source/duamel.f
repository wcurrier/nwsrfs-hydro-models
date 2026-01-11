C     Modified DUAMEL to support restart states
C     New arguments:
C       QPREV   - Previous Q values (M values) for warm start
C       USE_QPREV - Flag: 1 = use QPREV for warm start, 0 = cold start
C       QPREV_OUT - Output: last M values of Q for next restart
C
C     Modified DUAMEL - corrected restart logic
C     Modified DUAMEL - cleaner restart logic
C Modified DUAMEL - treating input as [QPREV|Q] concatenation
      SUBROUTINE DUAMEL(Q,UN1,UT,DT,N,MM,K,NTAU,QB,
     &                  QPREV,USE_QPREV,QPREV_OUT)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER A,B,I,J,M,IOR
      INTEGER N,MM,K,NTAU,USE_QPREV,M_SAVE
      INTEGER NQ,QIDX
      
      REAL Q(*), QB(*), QPREV(*), QPREV_OUT(*)
      REAL U1(MM)
      REAL QVAL
C
      M  = MM
      M_SAVE = MM
      NQ = N - MM
C
C----- Build unit hydrograph (unchanged code)
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
C----- Convolution: conceptually convolving with [QPREV|Q]
C     where QPREV has indices -M+1 to 0, Q has indices 1 to NQ
C
      IOC = N + NTAU

      DO 100 I=1,IOC
        QB(I) = 0.0
        
        DO 90 J=1,M
          IOR = J
          QIDX = I - J + 1
          
C         Determine which Q value to use
          IF (QIDX .GE. 1 .AND. QIDX .LE. NQ) THEN
C           Use current Q array
            QVAL = Q(QIDX)
          ELSE IF (QIDX .LE. 0 .AND. USE_QPREV .EQ. 1) THEN
C           Use QPREV array (QIDX ranges from 0 down to -M+1)
            QIDX = M_SAVE + QIDX
            IF (QIDX .GE. 1 .AND. QIDX .LE. M_SAVE) THEN
              QVAL = QPREV(QIDX)
            ELSE
              QVAL = 0.0
            END IF
          ELSE
            QVAL = 0.0
          END IF
          
          QB(I) = QB(I) + QVAL*U1(IOR)
   90   CONTINUE
  100 CONTINUE

C
C----- Save last M values of Q for restart
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
