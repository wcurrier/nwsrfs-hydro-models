C     ==============================================================================
C     Documentation NOTES from AWW
C     ==============================================================================
C      SUBROUTINE DUAMEL(Q,UN1,UT,DT,N,MM,K,NTAU,QB)
C      AWW adding documentation for input variables
C      Q      : TCI (total channel inflow) vector ... .ie unrouted
C      U1     : unit hydrograph vector
C      UN1    : unit hydrograph shape parameter (for gamma dist)
C      UT     : unit hydrograph scale parameter (for gamma dist)
C      DT     : timestep of the UH function (in days or fractions thereof)
C      N      : sim_length + uh_length ...ie length of U1
C      M      : max UH length?
C      QB     : routed flow vector
C      K      : 
C      NTAU   : 

C      AWW adding documentation for local variables
C      SP     : 
      
C     ==============================================================================
C     WRC MODIFICATION NOTES - RESTART CAPABILITY ADDITIONS (January 2026)
C     ==============================================================================
C     The following modifications enable unit hydrograph routing restart capability
C     for operational forecasting workflows. Changes allow:
C     - Warm start from previous routing states (last M flow values)
C     - Continuation of routing operations across simulation breaks
C     - Support for forecast ensemble generation from same initial conditions
C     
C     Key additions:
C     1. New input arguments: QPREV (previous M flow values), USE_QPREV (flag)
C     2. New output argument: QPREV_OUT (final M values for next restart)
C     3. Modified convolution logic to incorporate previous flow history
C     4. State capture logic for continuous restart capability
C     
C     Conceptual approach:
C     - Convolution operates on concatenated array [QPREV|Q]
C     - QPREV indices: -M+1 to 0 (previous simulation)
C     - Q indices: 1 to NQ (current simulation)
C     - This maintains routing continuity across simulation boundaries

C     New arguments:
C       QPREV     - Previous Q values (M values) for warm start
C       USE_QPREV - Flag: 1 = use QPREV for warm start, 0 = cold start
C       QPREV_OUT - Output: last M values of Q for next restart
C
C     =================================================================
C     WRC: RESTART INDEXING CONCEPT
C     =================================================================
C     Conceptual combined array for convolution: [QPREV|Q]
C
C     QPREV indices:  -M+1, -M+2, ..., -1,  0
C     Q indices:         1,    2, ..., NQ-1, NQ
C                        ^                   ^
C                        |                   |
C                   Start of Q          End of Q
C                                       (saved to QPREV_OUT)
C
C     During convolution at timestep I:
C       - Loop over unit hydrograph positions J=1 to M
C       - Compute QIDX = I - J + 1
C       - If QIDX >= 1: use Q(QIDX) from current simulation
C       - If QIDX <= 0: use QPREV(M + QIDX) from previous simulation
C       - This creates seamless routing across simulation boundaries
C     =================================================================
C     WRC: USAGE EXAMPLES
C     =================================================================
C     Cold Start (original behavior):
C       CALL DUAMEL(Q, UN1, UT, DT, N, MM, K, NTAU, QB,
C    &              QPREV_DUMMY, 0, QPREV_OUT)
C       - USE_QPREV = 0 ignores QPREV_DUMMY
C       - QPREV_OUT captures final M values for potential restart
C
C     Warm Start (restart from previous run):
C       CALL DUAMEL(Q, UN1, UT, DT, N, MM, K, NTAU, QB,
C    &              QPREV_FROM_PREVIOUS, 1, QPREV_OUT)
C       - USE_QPREV = 1 uses QPREV_FROM_PREVIOUS for routing history
C       - QPREV_OUT captures new final M values for next restart
C
C     =================================================================
C     WRC: TECHNICAL NOTES
C     =================================================================
C     1. Array sizing: M_SAVE stores original MM to ensure QPREV_OUT
C        always has correct size even if effective M shrinks due to
C        early UH decay
C
C     2. Index mapping: QPREV(M + QIDX) where QIDX ∈ [-M+1, 0]
C        ensures QPREV(1) corresponds to oldest value (QIDX = -M+1)
C        and QPREV(M) corresponds to newest value (QIDX = 0)
C
C     3. Edge case handling: When NQ < M_SAVE (very short simulation),
C        QPREV_OUT blends available Q values with carried-over QPREV
C        values to maintain full M-length state
C
C     4. Backward compatibility: Setting USE_QPREV=0 exactly replicates
C        original DUAMEL behavior (QVAL=0 for QIDX<=0)
C===========================================================
C
C     THIS SUBROUTINE PERFORM UNIT HYDROGRAPH ROUTING
C
      SUBROUTINE DUAMEL(Q,UN1,UT,DT,N,MM,K,NTAU,QB,
     &                  QPREV,USE_QPREV,QPREV_OUT)
     
C
C     WRC: ARGUMENT DOCUMENTATION (Updated for restart capability)
C     =================================================================
C     INPUTS:
C       Q      : TCI (total channel inflow) vector - unrouted flow (NQ values)
C       UN1    : Unit hydrograph shape parameter (gamma distribution)
C       UT     : Unit hydrograph scale parameter (gamma distribution) 
C       DT     : Timestep of the UH function (days or fractions thereof)
C       N      : sim_length + MM (includes UH tail for routing completion)
C       MM     : Unit hydrograph length (M values needed for convolution)
C       K      : UH computation flag (0 = skip UH build)
C       NTAU   : Lag parameter for routing
C       QPREV  : WRC: Previous M flow values from prior simulation (warm start)
C       USE_QPREV : WRC: Flag to enable warm start (1=warm, 0=cold)
C
C     OUTPUTS:
C       QB     : Routed flow vector (N + NTAU values)
C       QPREV_OUT : WRC: Last M values of Q for next restart
C
C     LOCAL VARIABLES:
C       U1     : Unit hydrograph ordinates (M values)
C       M      : Actual UH length used (may be < MM if UH decays early)
C       M_SAVE : WRC: Original MM value for array sizing
C       NQ     : WRC: Actual Q array length (N - MM)
C       QIDX   : WRC: Index into conceptual [QPREV|Q] array
C       QVAL   : WRC: Selected Q value during convolution
C     =================================================================

      IMPLICIT REAL (A-H,O-Z)
      INTEGER A,B,I,J,M,IOR
      INTEGER N,MM,K,NTAU,USE_QPREV,M_SAVE   ! WRC: Added USE_QPREV flag and M_SAVE
      INTEGER NQ,QIDX ! WRC: Added for restart indexing logic
      
      REAL Q(*), QB(*), QPREV(*), QPREV_OUT(*) ! WRC: Added QPREV and QPREV_OUT arrays
      REAL U1(MM)
      REAL QVAL ! WRC: Added to hold selected Q value during convolution
C
      M  = MM
      M_SAVE = MM ! WRC: Save original M for restart state sizing
      NQ = N - MM ! WRC: Actual length of Q array (N includes UH tail)
C
C----- WRC: Build unit hydrograph (unchanged from original)
C----- This section constructs the gamma-distribution unit hydrograph
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
C----- WRC: ===== MODIFIED CONVOLUTION WITH RESTART CAPABILITY (BEGIN) =====
C----- WRC: Convolution now operates on conceptual array [QPREV|Q]
C-----      QPREV covers indices -M+1 to 0 (from previous simulation)
C-----      Q covers indices 1 to NQ (current simulation)
C-----      This maintains routing continuity when restarting simulations
C
      IOC = N + NTAU

      DO 100 I=1,IOC
        QB(I) = 0.0
        
        DO 90 J=1,M
          IOR = J ! WRC: Unit hydrograph index
          QIDX = I - J + 1 ! WRC: Index into [QPREV|Q] concatenated array
          
C         WRC: Determine which Q value to use based on QIDX position
          IF (QIDX .GE. 1 .AND. QIDX .LE. NQ) THEN
C           WRC: QIDX in range [1, NQ] - use current Q array
            QVAL = Q(QIDX)
          ELSE IF (QIDX .LE. 0 .AND. USE_QPREV .EQ. 1) THEN
C           WRC: QIDX in range [-M+1, 0] and warm start enabled - use QPREV
C           WRC: Map QIDX to QPREV array: QIDX=0→M, QIDX=-1→M-1, etc.
            QIDX = M_SAVE + QIDX
            IF (QIDX .GE. 1 .AND. QIDX .LE. M_SAVE) THEN
              QVAL = QPREV(QIDX)
            ELSE
              QVAL = 0.0  ! WRC: Out of bounds safety
            END IF
          ELSE
C           WRC: Cold start or out of range - use zero (original behavior)
            QVAL = 0.0
          END IF
          
          QB(I) = QB(I) + QVAL*U1(IOR)  ! WRC: Accumulate convolution
   90   CONTINUE
  100 CONTINUE
C
C----- WRC: ===== MODIFIED CONVOLUTION WITH RESTART CAPABILITY (END) =====

C
C----- WRC: ===== SAVE FINAL M VALUES FOR RESTART (BEGIN) =====
C----- WRC: Capture the last M flow values for use in next simulation
C----- WRC: This enables continuous routing across simulation boundaries
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
C----- WRC: ===== SAVE FINAL M VALUES FOR RESTART (END) =====
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
