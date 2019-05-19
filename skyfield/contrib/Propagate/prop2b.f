C$Procedure      PROP2B ( Propagate a two-body solution )
 
      SUBROUTINE PROP2B ( GM, PVINIT, DT, PVPROP )
 
C$ Abstract
C
C     Given a central mass and the state of massless body at time t_0,
C     this routine determines the state as predicted by a two-body
C     force model at time t_0 + DT.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     CONIC
C     EPHEMERIS
C     UTILITY
C
C$ Declarations
 
      DOUBLE PRECISION      GM
      DOUBLE PRECISION      PVINIT ( 6 )
      DOUBLE PRECISION      DT
      DOUBLE PRECISION      PVPROP ( 6 )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     GM         I   Gravity of the central mass.
C     PVINIT     I   Initial state from which to propagate a state.
C     DT         I   Time offset from initial state to propagate to.
C     PVPROP     O   The propagated state.
C
C$ Detailed_Input
C
C     GM         is the gravitational constant G times the mass M of the
C                central body.
C
C     PVINIT     is the state at some specified time relative to the
C                central mass.  The mass of the object is assumed to
C                be negligible when compared to the central mass.
C
C     DT         is a offset in time from the time of the initial
C                state to which the two-body state should be
C                propagated. (The units of time and distance must be
C                the same in GM, PVINIT, and DT).
C
C$ Detailed_Output
C
C     PVPROP     is the two-body propagation of the initial state
C                DT units of time past the epoch of the initial state.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If GM is not positive, the error SPICE(NONPOSITIVEMASS) will
C        be signalled.
C
C     2) If the position of the initial state is the zero vector, the
C        error SPICE(ZEROPOSITION) will be signalled.
C
C     3) If the velocity of the initial state is the zero vector, the
C        error SPICE(ZEROVELOCITY) will be signalled.
C
C     4) If the cross product of the position and velocity of PVINIT
C        has squared length of zero, the error SPICE(NONCONICMOTION)
C        will be signalled.
C
C     5) The value of DT must be "reasonable".  In other words, DT
C        should be less than 10**20 seconds for realistic solar system
C        orbits specified in the MKS system.  (The actual bounds
C        on DT are much greater but require substantial computation.)
C        The "reasonableness" of DT is checked at run-time.  If DT is
C        so large that there is a danger of floating point overflow
C        during computation, the error SPICE(DTOUTOFRANGE) is
C        signalled and a message is generated describing the problem.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine uses a universal variables formulation for the
C     two-body motion of an object in orbit about a central mass. It
C     propagates an initial state to an epoch offset from the
C     epoch of the initial state by time DT.
C
C     This routine does not suffer from the finite precision
C     problems of the machine that are inherent to classical
C     formulations based on the solutions to Kepler's equation:
C
C           n( t - T ) = E - e Sin(E)         elliptic case
C           n( t - T ) = e sinh(F) - F        hyperbolic case
C
C     The derivation used to determine the propagated state is a
C     slight variation of the derivation in Danby's book
C     `Fundamentals of Celestial Mechanics' [1] .
C
C$ Examples
C
C     When the eccentricity of an orbit is near 1, and the epoch
C     of classical elements is near the epoch of periapse, classical
C     formulations that propagate a state from elements tend to
C     lack robustness due to the finite precision of floating point
C     machines. In those situations it is better to use a universal
C     variables formulation to propagate the state.
C
C     By using this routine, you need not go from a state to elements
C     and back to a state. Instead, you can get the state from an
C     initial state.
C
C     If PV is your initial state and you want the state 3600
C     seconds later, the following call will suffice.
C
C          Look up GM somewhere
C
C          DT = 3600.0D0
C
C          CALL PROP2B ( GM, PV, DT, PVDT )
C
C     After the call, PVDT will contain the state of the
C     object 3600 seconds after the time it had state PV.
C
C$ Restrictions
C
C     Users should be sure that GM, PVINIT and DT are all in the
C     same system of units ( for example MKS ).
C
C$ Literature_References
C
C     [1] `Fundamentals of Celestial Mechanics', Second Edition
C         by J.M.A. Danby;  Willman-Bell, Inc., P.O. Box 35025
C         Richmond Virginia;  pp 168-180
C
C$ Author_and_Institution
C
C     W.L. Taber     (JPL)
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.1.0, 31-AUG-2005 (NJB)
C
C        Updated to remove non-standard use of duplicate arguments
C        in VSCL call.
C
C-    SPICELIB Version 2.0.1, 22-AUG-2001 (EDW)
C
C        Corrected ENDIF to END IF.
C     
C-    Spicelib Version 2.0.0 16-May-1995  (WLT)
C
C        The initial guess at a solution to Kepler's equation was
C        modified slightly and a loop counter was added to the
C        bisection loop together with logic that will force termination
C        of the bisection loop.
C
C-    Spicelib Version 1.0.0, 10-Mar-1992 (WLT)
C
C
C-&
 
C$ Index_Entries
C
C      Propagate state vector using two-body force model
C
C-&

C$ Revisions
C
C-    SPICELIB Version 2.1.0, 31-AUG-2005 (NJB)
C
C        Updated to remove non-standard use of duplicate arguments
C        in VSCL call.
C
C-& 
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      BRCKTD
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      VDOT
      DOUBLE PRECISION      VNORM
 
      INTEGER               BRCKTI
 
      LOGICAL               RETURN
      LOGICAL               VZERO
 
C
C     Local Parameters
C
      INTEGER               BUFSIZ
      PARAMETER           ( BUFSIZ = 3 )
 
      INTEGER               MAXBIT
      PARAMETER           ( MAXBIT = 64 )

C
C     Local variables
C
C
C     The following quantities are needed in the solution of Kepler's
C     equation and in the propagation of the input state.  They are
C     described as they are introduced in the code below.
C
      DOUBLE PRECISION      B
      DOUBLE PRECISION      B2RV
      DOUBLE PRECISION      BQ
      DOUBLE PRECISION      BR
      DOUBLE PRECISION      BR0
      DOUBLE PRECISION      C0
      DOUBLE PRECISION      C1
      DOUBLE PRECISION      C2
      DOUBLE PRECISION      C3
      DOUBLE PRECISION      E
      DOUBLE PRECISION      EQVEC  ( 3 )
      DOUBLE PRECISION      F
      DOUBLE PRECISION      FX2
      DOUBLE PRECISION      H2
      DOUBLE PRECISION      HVEC   ( 3 )
      DOUBLE PRECISION      KFUN
      DOUBLE PRECISION      KFUNL
      DOUBLE PRECISION      KFUNU
      DOUBLE PRECISION      LOGBND
      DOUBLE PRECISION      LOWER
      DOUBLE PRECISION      OLDX
      DOUBLE PRECISION      PC
      DOUBLE PRECISION      PCDOT
      DOUBLE PRECISION      POS    ( 3 )
      DOUBLE PRECISION      Q
      DOUBLE PRECISION      QOVR0
      DOUBLE PRECISION      R0
      DOUBLE PRECISION      RV
      DOUBLE PRECISION      TMPVEC ( 3 )
      DOUBLE PRECISION      UPPER
      DOUBLE PRECISION      VC
      DOUBLE PRECISION      VCDOT
      DOUBLE PRECISION      VEL    ( 3 )
      DOUBLE PRECISION      X
      DOUBLE PRECISION      X2
      DOUBLE PRECISION      X3
 
C
C     The variables below store intermediate results that can be
C     reused if PVINIT is supplied more than once to this routine.
C     In this way, the number of redundant computations can be reduced.
C
      DOUBLE PRECISION      SAVEPV ( 6, BUFSIZ )
      DOUBLE PRECISION      SAVEGM (    BUFSIZ )
      DOUBLE PRECISION      SBOUND (    BUFSIZ )
      DOUBLE PRECISION      SB2RV  (    BUFSIZ )
      DOUBLE PRECISION      SBQ    (    BUFSIZ )
      DOUBLE PRECISION      SBR0   (    BUFSIZ )
      DOUBLE PRECISION      SF     (    BUFSIZ )
      DOUBLE PRECISION      SQOVR0 (    BUFSIZ )
 
      INTEGER               NEWEST (    BUFSIZ )
      INTEGER               NSAVED
 
C
C     Variables used to bracket X in our solution of Kepler's equation.
C
      DOUBLE PRECISION      BOUND
      DOUBLE PRECISION      FIXED
      DOUBLE PRECISION      LOGDPM
      DOUBLE PRECISION      LOGF
      DOUBLE PRECISION      LOGMXC
      DOUBLE PRECISION      MAXC
      DOUBLE PRECISION      ROOTF
 
      INTEGER               BUMPED
      INTEGER               I
      INTEGER               K
      INTEGER               MOSTC
      INTEGER               LCOUNT
  
      LOGICAL               NEW
 
C
C     Save everything.
C
      SAVE
 
C
C     Initial values
C
      DATA                  NSAVED  / 0       /
      DATA                  NEWEST  / 1, 2, 3 /
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      ELSE
         CALL CHKIN ( 'PROP2B' )
      END IF
 
C
C     Life will be easier if we use POS and VEL to hold the state.
C
      POS(1) = PVINIT(1)
      POS(2) = PVINIT(2)
      POS(3) = PVINIT(3)
 
      VEL(1) = PVINIT(4)
      VEL(2) = PVINIT(5)
      VEL(3) = PVINIT(6)
 
C
C     If we propagate many states from the same initial state,
C     most of the variables used to propagate the state will
C     not change in value.
C
C     To save time needed to compute these variables, we recompute
C     variables that depend upon the initial state only when the
C     initial state is not one of those already buffered by this
C     routine.
C
C     Determine whether or not this GM and state are the same as the
C     one of those already buffered.  Note that we look through the
C     saved states and GM from the most recently input values of PVINIT
C     and GM to the oldest saved state and GM.
C
C     NEWEST(1)  contains the most recently input initial conditions
C     NEWEST(2)  contains the next most recently input intial conditions
C                etc.
C
C     Also note that when this routine starts up there will be no
C     buffered states or GMs.  Every time we encounter a new state, we
C     will increment the number of saved states NSAVED until we have
C     BUFSIZ states buffered.  From that point on, when a new state is
C     encountered we will overwrite the oldest buffered state.
C
      I   =  0
      NEW = .TRUE.
 
      DO WHILE ( ( I .LT. NSAVED ) .AND. NEW )
 
         I   = I + 1
         K   = NEWEST(I)
 
         NEW =     ( PVINIT(1) .NE. SAVEPV(1,K) )
     .        .OR. ( PVINIT(2) .NE. SAVEPV(2,K) )
     .        .OR. ( PVINIT(3) .NE. SAVEPV(3,K) )
     .        .OR. ( PVINIT(4) .NE. SAVEPV(4,K) )
     .        .OR. ( PVINIT(5) .NE. SAVEPV(5,K) )
     .        .OR. ( PVINIT(6) .NE. SAVEPV(6,K) )
     .        .OR. ( GM        .NE. SAVEGM(  K) )
 
      END DO
 
 
      IF ( .NOT. NEW ) THEN
C
C        We update the order vector NEWEST so that the state being
C        used this time becomes the "youngest" state.
C
         K      = I
         BUMPED = NEWEST(K)
 
         DO I = K, 2, -1
            NEWEST(I) = NEWEST(I-1)
         END DO
 
         NEWEST(1) = BUMPED
         K         = BUMPED
 
C
C        Now look up all of the other saved quantities.
C
         B2RV   = SB2RV (K)
         BOUND  = SBOUND(K)
         BQ     = SBQ   (K)
         BR0    = SBR0  (K)
         F      = SF    (K)
         QOVR0  = SQOVR0(K)
 
      ELSE
 
C
C        We have a new state, new GM or both.  First let's make sure
C        there is nothing obviously wrong with them.  (We buffer
C        only states, GMs and intermediate values that are "good.")
C        First check for nonpositive mass.
C
         IF ( GM .LE. 0.0D0 ) THEN
            CALL SIGERR ( 'SPICE(NONPOSITIVEMASS)'   )
            CALL CHKOUT ( 'PROP2B'                   )
            RETURN
         END IF
 
C
C        Next for a zero position vector
C
         IF ( VZERO(POS) ) THEN
            CALL SIGERR ( 'SPICE(ZEROPOSITION)'   )
            CALL CHKOUT ( 'PROP2B'                )
            RETURN
         END IF
 
C
C        Finally for a zero velocity vector
C
         IF ( VZERO(VEL) ) THEN
            CALL SIGERR ( 'SPICE(ZEROVELOCITY)'   )
            CALL CHKOUT ( 'PROP2B'                )
            RETURN
         END IF
 
C
C        Obvious problems have been checked. Here are the relevant
C        equations. Let ...
C
C           GM        be the gravitational attraction of the central
C                     mass.
C
C           POS and   be the initial position and velocity respectively
C           VEL       of the orbiting object.
C
C           R0       be the magnitude of the position vector POS
C
C           RV       be the value of the dot product  POS * VEL
C
         R0 = VNORM ( POS      )
         RV = VDOT  ( POS, VEL )
 
C
C        Let HVEC be the specific angular momentum vector and let Q be
C        the distance at periapse.
C
C                   1)    HVEC  =   POS  x  VEL
C
C                                       2
C                   2)    H2    = |HVEC|  =  GM*(1+E)*Q
C
C
         CALL VCRSS ( POS,  VEL,     HVEC )
         H2 = VDOT  ( HVEC, HVEC          )
 
 
C
C        Let's make sure we are not in the pathological case of
C        rectilinear motion.
C
         IF ( H2 .EQ. 0 ) THEN
            CALL SIGERR ( 'SPICE(NONCONICMOTION)'   )
            CALL CHKOUT ( 'PROP2B'                  )
            RETURN
         END IF
 
C
C        Let E be the eccentricity of the orbit.
C
C        Let QVEC be the unit vector that points toward perihelion, and
C        let EQVEC be QVEC scaled by E.
C
C                                   VEL X HVEC      POS
C                    1)  E*QVEC  =  ----------  -   ---
C                                       GM           R0
C
C
C                                         VEL X HVEC      POS
C                    2)  E       = NORM ( ----------  -   --- )
C                                            GM            R0
C
C
         CALL VCRSS ( VEL,       HVEC,                      TMPVEC  )
         CALL VLCOM ( 1.0D0/GM,  TMPVEC,  -1.0D0/R0,  POS,  EQVEC   )
 
         E  = VNORM ( EQVEC )
 
C
C        Solve the equation H2 = GM*Q*(1+E) for Q.
C
         Q  = H2 / ( GM * (1+E) )
 
C
C        From the discussion of the universal variables formulation in
C        Danby's book on pages 174 and 175 (see the reference listed
C        above) you can show that by making the substitutions
C
C              F  =  1 - E
C
C        and
C
C                       _____
C                      /  Q
C              S  =   / -----    X   = B * X
C                   \/   GM
C
C
C        that DT satisfies the universal variables Kepler's equation:
C
C                                   2     2     2        2
C              DT =  B*R0*X*C_1( F*X ) + B *RV*X C_2( F*X )
C
C                                               3        2
C                                      +   B*Q*X C_3( F*X )
C
C                 =  KFUN( X )
C
C        (where C_k is used to denote the Stumpff functions. This is
C        the universal variables formulation of Kepler's equation.
C        KFUN is our abbreviation for "Kepler function.")
C
C        (One might wonder, "Why make such a change of variables?"
C        By making this substitution early in the derivation supplied
C        in Danby's book, you can always deal with physically
C        meaningful quantities --- the pure numeric value of F and the
C        distance of periapse.  Thus one does not need to be concerned
C        about infinite or negative semi-major axes or with discussing
C        how to interpret these somewhat artificial artifacts of the
C        classical derivations for two body motion.)
C
C        Given the unique X for which this Kepler's equation is
C        satisfied, we can compute the state of the orbiting object
C        at a time DT past the epoch of the state POS and VEL.
C        Evidently we will need the constants:
C
         F      = 1.0D0 - E
         B      = DSQRT ( Q / GM )
         BR0    = B*R0
         B2RV   = B*B*RV
         BQ     = B*Q
 
C
C        The state corresponding to the value of X that solves this
C        equation is given by
C
C              PC * POS + VC * VEL              ( position )
C
C        and
C
C              PCDOT * POS + VCDOT * VEL        ( velocity )
C
C        where
C                                            2        2
C           ( 1 )    PC    =  1  -  ( Q/R0 )X C_2( F*X )
C
C                                            3        2
C           ( 2 )    VC    =  DT -  ( B*Q  )X C_3( F*X )
C
C
C                                       Q               2
C           ( 3 )    PCDOT =     -  ( ------ ) X C_1( F*X )
C                                     B*R*R0
C
C                                      B*Q     2        2
C           ( 4 )    VCDOT =  1  -  (  ---  ) X C_2( F*X )
C                                      B*R
C
C        Here R denotes the distance from the center of CP*POS + CV*VEL
C        It turns out that R can be computed as:
C
C                                        2     2             2
C           ( 5 )   B*R    = B*R0 C_0(F*X ) + B *RV X C_1(F*X )
C
C                                                 2       2
C                                        +   B*Q X C_2(F*X )
C
C
C        Therefore we will also need the constant
C
         QOVR0   = Q / R0
 
C
C        We will have to find the unique value of X such that
C
C             DT = KFUN ( X )
C
C        where KFUN stands for the "Kepler function" defined by the
C        equation below:
C
C                                   2
C        KFUN(X) =   B*R0*X * C_1(FX )
C
C                   2     2        2
C                + B *RV*X * C_2(FX )
C
C                         3        2
C                +   B*Q*X * C_3(FX )
C
C
C        (There is a unique solution to this equation. KFUN(X) is
C        unbounded above and below and is an increasing function
C        over all real X for all non-rectilinear orbits. To see this
C        we note that the variable X is a function of DT and is given
C        by the integral from 0 to DT of the differential:
C
C                   dt
C                 ------
C                 B*R(t)
C
C        where R(t) is the range of the body as a function of time.
C        Therefore X is an increasing function of DT, and DT must
C        also be an increasing function of X.
C
C        Thus, there is a unique value of X  that solves this
C        equation).
C
C        If F is less than zero, we can have the computation of C0,...
C        overflow.  This is because for X < 0
C
C
C               C_0(X) = COSH( DSQRT(-X) )
C
C               C_1(X) = SINH( DSQRT(-X) )
C                        -----------------
C                              DSQRT(-X)
C
C
C
C        and from the recursion relationship we know that
C
C
C               C_2(X) =  ( 1/0! - C_0(X) ) / X
C
C               C_3(X) =  ( 1/1! - C_1(X) ) / X
C
C
C                         1 - COSH( DSQRT(-X) )
C               C_2(X) = ------------------------
C                                  X
C
C                         1  - SINH( DSQRT(-X) ) / DSQRT(-X)
C               C_3(X) = -----------------------------------
C                                    X
C
C        Clearly for negative values of F*X*X having large magnitude,
C        it is easy to get an overflow.
C
C        In the case when F is less than 0 we choose X so that we can
C        compute all of the following:
C
C               | COEF_0 * X**0 * C_0(FX**2) |
C
C               | COEF_1 * X**1 * C_1(FX**2) |
C
C               | COEF_2 * X**2 * C_2(FX**2) |
C
C               | COEF_3 * X**3 * C_3(FX**2) |
C
C
C         where COEF_n are coefficients that will be used in forming
C         linear combinations of X**n C_n(FX**2) terms.
C
C         The variable portion of the last 3 terms above can be
C         rewritten as:
C
C
C                                   SINH ( DSQRT(-F)*|X| )
C        | X**1 * C_1(FX**2) |  =   ----------------------
C                                          DSQRT(-F)
C
C
C
C                                   1 - COSH( DSQRT(-F)*|X| )
C        | X**2 * C_2(FX**2) |  =  ----------------------------
C                                             -F
C
C
C                                  DSQRT(-F)*|X|   - SINH(DSQRT(-F)*|X|)
C        | X**3 * C_3(FX**2) |  =  -------------------------------------
C                                              F*DSQRT(-F)
C
C
C        For large |X| the absolute values of these expressions are well
C        approximated by
C
C                                         0.0
C               COSH( DSQRT(-F)|X| ) * |F|
C
C                                         -0.5
C               SINH( DSQRT(-F)|X| ) * |F|
C
C                                         -1.0
C               COSH( DSQRT(-F)|X| ) * |F|
C
C                                         -1.5
C               SINH( DSQRT(-F)|X| ) * |F|
C
C
C        For large |X| the logarithms of these expressions are well
C        approximated by:
C
C
C               DSQRT(-F)|X| - LOG(2) - 0.0*LOG(-F)
C
C               DSQRT(-F)|X| - LOG(2) - 0.5*LOG(-F)
C
C               DSQRT(-F)|X| - LOG(2) - 1.0*LOG(-F)
C
C               DSQRT(-F)|X| - LOG(2) - 1.5*LOG(-F)
C
C        respectively.
C
C
C        To ensure that we can form a linear combination of these terms
C        we will require that:
C
C
C           |COEF_N*X**N * C_N(FX**2)| < DPMAX / 4
C
C
C
C        for N=0,1,2,3.  This is equivalent to
C
C              LOG ( X**N * C_N(FX**2) )   <      LOG ( DPMAX )
C            + LOG (|COEF_N|)                   - 2 LOG ( 2     )
C
C
C
C        or
C
C              LOG ( X**N * C_N(FX**2) )   <      LOG ( DPMAX    )
C                                             -   LOG ( |COEF_N| )
C                                             - 2*LOG ( 2        ).
C
C
C        Replacing the left hand side with the magnitude expressions
C        computed above we have:
C
C            DSQRT(-F)|X| - LOG(2) - N*0.5*LOG( -F )  <   LOG ( DPMAX  )
C                                                      -  LOG (|COEF_N|)
C                                                      -2*LOG ( 2      )
C
C         So that:
C
C
C            |X|  <    {   LOG ( DPMAX  )
C                        - LOG (|COEF_N|)
C                        - LOG (  2     )
C                        + LOG ( -F     )*N*0.5 } / DSQRT(-F)
C
C         Let MAXC be the maximum of 1.0D0 and the various coefficients
C         of the Stumpff functions.  We can then set our absolute value
C         bound on X to be:
C
C
C             MIN        LOG(DPMAX/2) - LOG(MAXC) + (n/2)LOG(-F)
C            n = 0,3  {  -----------------------------------------  }
C                               DSQRT(-F)
C
C        (Actually we know that the minimum must occur for n = 0 or
C        for n = 3).
C
C
         MAXC   = MAX   ( 1.0D0, DABS(BR0),
     .                           DABS(B2RV),
     .                           DABS(BQ),
     .                           DABS(QOVR0/BQ) )
 
         IF ( F .LT. 0 ) THEN
 
            LOGMXC = LOG   ( MAXC            )
            LOGDPM = LOG   ( DPMAX() / 2.0D0 )
 
            FIXED  = LOGDPM  - LOGMXC
 
            ROOTF  = DSQRT ( - F             )
            LOGF   = LOG   ( - F             )
 
            BOUND  = MIN (  (FIXED              )/ROOTF,
     .                      (FIXED + 1.5D0*LOGF )/ROOTF )
 
C
C           Note that in the above, we can always perform the division
C           by ROOTF.  To see this we note that -F is at least the
C           machine precision (we got it by subtracting E from 1.)
C           Thus its square root is a reasonably large number (if F is
C           10**-N then ROOTF is 10**(-N/2) )  The value of FIXED is
C           about 3*M where M is the largest exponent such that 2**M
C           is representable on the host machine.  Thus BOUND is at
C           worst M*10**(N/2)  This will always be computable.
C
        ELSE
 
C
C
C           In the case when F is non-negative we must be sure we
C           can compute all of the following.
C
C               | COEF_0 * X**0 * C_0(FX**2) | < | COEF_0          |
C
C               | COEF_1 * X**1 * C_1(FX**2) | < | COEF_1*|X|      |
C
C               | COEF_2 * X**2 * C_2(FX**2) | < | COEF_2*X**2 / 2 |
C
C               | COEF_3 * X**3 * C_3(FX**2) | < | COEF_3*X**3 / 6 |
C
C           If we assume that COEF_0 is computable, all of these are
C           bounded above by:
C
C                       | MAX(COEF_1,...COEF_3) * X**3 / 6 |
C
C           We want to make sure we can add these terms so we need to
C           make sure that
C
C              | MAX(COEF_1,...,COEF_3) * X**3 / 6 | < DPMAX() / 4.
C
C           Thus we need:
C
C              |X**3| <          1.5*DPMAX / MAX(COEF_1,...,COEF_3)
C              |X|    <  DCBRT ( 1.5*DPMAX / MAX(COEF_1,...,COEF_3) )
C
C           (We'll use logarithms to compute the upper bound for |X|.)
C
            LOGBND =      ( DLOG(1.5D0) + DLOG(DPMAX()) - DLOG(MAXC) )
     .               /      3.0D0
 
            BOUND  = DEXP ( LOGBND )
 
         END IF
C
C        All the obvious problems have been checked, move everybody
C        on the list down and put the new guy on top of the list.
C
 
         NSAVED = BRCKTI ( NSAVED+1, 1, BUFSIZ )
         BUMPED = NEWEST ( NSAVED              )
 
         DO I = NSAVED, 2, -1
            NEWEST(I) = NEWEST(I-1)
         END DO
 
         NEWEST(1) = BUMPED
         K         = BUMPED
 
         SAVEPV(1,K) = PVINIT(1)
         SAVEPV(2,K) = PVINIT(2)
         SAVEPV(3,K) = PVINIT(3)
         SAVEPV(4,K) = PVINIT(4)
         SAVEPV(5,K) = PVINIT(5)
         SAVEPV(6,K) = PVINIT(6)
         SAVEGM(  K) = GM
 
C
C        Finally we save the results of all of the above
C        computations so that we won't have to do them again,
C        if this initial state and GM are entered again.
C
         SB2RV (K) = B2RV
         SBOUND(K) = BOUND
         SBQ   (K) = BQ
         SBR0  (K) = BR0
         SF    (K) = F
         SQOVR0(K) = QOVR0
 
      END IF
C
C
C     We are now ready to find the unique value of X such that
C
C             DT = KFUN ( X )
C
C     First we must bracket the root. The basic idea is this:
C
C     1) KFUN(0) = 0 so we will let one endpoint of our initial
C        guess of a bracketing interval be 0.
C
C     2) We get our initial guess at the other endpoint of the
C        bracketing interval by recalling that
C
C                   dt
C         dX  =   ------
C                 B*R(t)
C
C        From this observation it follows that
C
C                   DT
C          X  <  -------
C                   B*Q
C
C        Thus the solution to
C
C             DT = KFUN ( X )
C
C        Satisifies
C
C                     DT
C         0 < X  <  -------
C                    B*Q
C
C
C        We now have a guess at a bracketing interval. In the case
C        DT is positive it looks like
C
C                0        X
C         -------[--------]-----------------------------
C
C        This is ok mathematically, but due to rounding etc it is
C        conceivable that we might not have bracketed the root.
C        We check and if not we will double the
C        endpoint farthest from zero and call this X, and make
C        the other endpoint the old value of X.
C
C
C                0
C         -------+--------[--------]--------------------
C
C
C        We continue this process ...
C
C                0
C         -------+-----------------[-----------------]--
C
C        ...until the root is bracketed. (One shift is certain
C        to do the job).
C
C        If we perform this interval shift, we will have to take
C        care that X does not run out of the domain for which
C        we can safely compute KFUN.  Thus we will make sure that
C        the endpoints of these shifted intervals always stay safely
C        inside the domain for which KFUN can be computed.
C
      X    = DT / BQ
      X    = BRCKTD ( X, -BOUND, BOUND )
      FX2  = F*X*X
 
      CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
      KFUN =   X*(BR0*C1 + X*( B2RV*C2 + X*(BQ*C3) ) )
 
      IF ( DT .LT. 0 ) THEN
 
         UPPER = 0
         LOWER = X
 
         DO WHILE ( KFUN .GT. DT )
 
            UPPER = LOWER
            LOWER = LOWER * 2.0D0
            OLDX  = X
            X     = BRCKTD ( LOWER, -BOUND, BOUND )
 
C
C           Make sure we are making progress. (In other words make sure
C           we don't run into the boundary of values that X can assume.
C           If we do run into the boundary, X will be unchanged and
C           there's nothing further we can do.  We'll have to call it
C           quits and tell the user what happened.)
C
            IF ( X .EQ. OLDX ) THEN
 
               FX2 = F*BOUND*BOUND
 
               CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
               KFUNL = -BOUND*(BR0*C1 - BOUND*(B2RV*C2 - BOUND*BQ*C3))
               KFUNU =  BOUND*(BR0*C1 + BOUND*(B2RV*C2 + BOUND*BQ*C3))
 
               CALL SETMSG ( 'The input delta time (DT) has a '       //
     .                       'value of #.  This is beyond the '       //
     .                       'range of DT for which we can '          //
     .                       'reliably propagate states. The limits ' //
     .                       'for this GM and initial state are '     //
     .                       'from # to #. '                          )
               CALL ERRDP  ( '#', DT                                  )
               CALL ERRDP  ( '#', KFUNL                               )
               CALL ERRDP  ( '#', KFUNU                               )
               CALL SIGERR ( 'SPICE(DTOUTOFRANGE)'                    )
               CALL CHKOUT ( 'PROP2B'                                 )
               RETURN
            END IF
 
            FX2   = F*X*X
 
            CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
            KFUN = X*(BR0*C1 + X*( B2RV*C2 + X*(BQ*C3) ) )
 
         END DO
 
 
      ELSE IF ( DT .GT. 0 ) THEN
 
         LOWER = 0
         UPPER = X
 
         DO WHILE ( KFUN .LT. DT )
 
            LOWER = UPPER
            UPPER = UPPER * 2.0D0
            OLDX  = X
            X     = BRCKTD ( UPPER, -BOUND, BOUND )
C
C           Make sure we are making progress.
C
            IF ( X .EQ. OLDX ) THEN
 
               FX2 = F*BOUND*BOUND
 
               CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
               KFUNL = -BOUND*(BR0*C1 - BOUND*(B2RV*C2 - BOUND*BQ*C3))
               KFUNU =  BOUND*(BR0*C1 + BOUND*(B2RV*C2 + BOUND*BQ*C3))
 
               CALL SETMSG ( 'The input delta time (DT) has a '       //
     .                       'value of #.  This is beyond the '       //
     .                       'range of DT for which we can '          //
     .                       'reliably propagate states. The limits ' //
     .                       'for this GM and initial state are '     //
     .                       'from # to #. '                          )
               CALL ERRDP  ( '#', DT                                  )
               CALL ERRDP  ( '#', KFUNL                               )
               CALL ERRDP  ( '#', KFUNU                               )
               CALL SIGERR ( 'SPICE(DTOUTOFRANGE)'                    )
               CALL CHKOUT ( 'PROP2B'                                 )
               RETURN
            END IF
 
            FX2  = F*X*X
 
            CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
            KFUN = X*(BR0*C1 + X*( B2RV*C2 + X*BQ*C3 ) )
 
         END DO
 
      ELSE
 
         CALL VEQUG  ( PVINIT, 6, PVPROP )
         CALL CHKOUT ( 'PROP2B'          )
         RETURN
 
      END IF
 
C
C     Ok. We've bracketed the root.  Now for lack of anything more
C     clever, we just bisect to find the solution.
C
C     We add a loop counter so that we can ensure termination of the
C     loop below.
C
C     On some systems the computed midpoint is stored in an extended
C     precision register.  Thus the midpoint is always different from
C     UPPER and LOWER.  Yet when the new value of LOWER and UPPER
C     are assigned UPPER and LOWER do not change and hence the
C     loop fails to terminate.  With the loop counter we force
C     termination of the loop.
C
      X    = MIN ( UPPER, MAX ( LOWER, (LOWER + UPPER)/2.0D0) )
      FX2  = F*X*X
 
      CALL STMP03 ( FX2, C0, C1, C2, C3 )
      LCOUNT = 0
      MOSTC  = 1000
 
      DO WHILE (         X      .GT. LOWER
     .           .AND.   X      .LT. UPPER
     .           .AND.   LCOUNT .LT. MOSTC )
 
         KFUN = X*(BR0*C1 + X*( B2RV*C2 + X*BQ*C3 ) )
 
         IF ( KFUN .GT. DT ) THEN
            UPPER = X
         ELSE IF ( KFUN .LT. DT ) THEN
            LOWER = X
         ELSE
            UPPER = X
            LOWER = X
         END IF
C
C        As soon as the bracketting values move away from
C        zero we can modify the count limit.
C
         IF ( MOSTC .GT. MAXBIT ) THEN
 
            IF ( UPPER .NE. 0.0D0 .AND. LOWER .NE. 0.0D0 ) THEN
               MOSTC   = MAXBIT
               LCOUNT = 0
            END IF
 
         END IF
 
         X    = MIN ( UPPER, MAX ( LOWER, (LOWER + UPPER)/2.0D0 ) )
         FX2  = F*X*X
 
         CALL STMP03 ( FX2, C0, C1, C2, C3 )
 
         LCOUNT = LCOUNT + 1
 
      END DO
 
 
C
C     With X in hand we simply compute BR, PC, VC, PCDOT and VCDOT
C     described in equations (1) --- (5) above. (Note, by our choice
C     of BOUND above, one can show that none of the computations
C     below can cause an overflow).
C
      X2    = X *X
      X3    = X2*X
      BR    = BR0*C0 + X*( B2RV*C1  + X*(BQ*C2) )
 
 
      PC    = 1.0D0     -      QOVR0        * X2 * C2
      VC    = DT        -      BQ           * X3 * C3
      PCDOT =           -     (QOVR0 / BR ) * X  * C1
      VCDOT = 1.0D0     -     (BQ    / BR ) * X2 * C2
 
 
C
C     ... and compute the linear combinations needed to get PVPROP
C
      CALL VLCOM ( PC,    POS, VC,    VEL,       PVPROP(1) )
      CALL VLCOM ( PCDOT, POS, VCDOT, VEL,       PVPROP(4) )
 
      CALL CHKOUT ( 'PROP2B' )
      RETURN
 
      END
 
