#include "symbol.inc"
!*******************************************************************
!
! this module contains all subroutines required to calculate
! the core level eigenvalues (initial state core levels)
!
!*******************************************************************

! "new" is the only version recommened
! in this case the local part of the PAW potentials is modified
! to correspond to the desired valency
! in the all version a local potential is added by the routine
! CL_TEST_POT

#define new


!#define RCORE_MAIN 0.5

MODULE cl
    
    USE prec
    IMPLICIT NONE
    INTEGER, PRIVATE,SAVE :: ICORELEVEL  ! 0=approximate initial state cl. shift
                                         ! 1=calculate accurate initial state shift
                                         ! 2=final state core level shift
    INTEGER, PRIVATE,SAVE :: NT_CL       ! type of ion for which final state shifts
                                         ! are calculated
    INTEGER, PRIVATE,SAVE :: N_CL        ! main quantum number
    INTEGER, PRIVATE,SAVE :: L_CL        ! l quantum number
    REAL(q), PRIVATE,SAVE :: Z_CL        ! occupancy
    INTEGER, PRIVATE,SAVE :: CL_MAXNL=20 ! maximal dimension for CL_SHIFT

    REAL(q),ALLOCATABLE,SAVE :: CL_SHIFT(:,:)
                                    !  calculated core level shifts
    COMPLEX(q),ALLOCATABLE,SAVE :: CVTEST(:)
    RGRID,ALLOCATABLE,SAVE :: CHTEST(:)

    CONTAINS

!*******************************************************************
!
! the routine init_cl_shift reads the line ICORELEVEL from
! INCAR and allocates the CL_SHIFT array
!
!*******************************************************************

    SUBROUTINE INIT_CL_SHIFT(IU5, IU0, NIONS, NTYP )
      IMPLICIT NONE 

      INTEGER :: IU5,IU0 ! input unit
      INTEGER :: NIONS   ! number of ions
      INTEGER :: NTYP    ! number of types

      LOGICAL :: LOPEN,LDUM
      INTEGER :: IDUM, N, IERR
      REAL(q) :: RDUM
      COMPLEX(q)  :: CDUM
      CHARACTER*1 :: CHARAC

      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')
! core level shifts
      ICORELEVEL=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'ICORELEVEL','=','#',';','I', &
     &            ICORELEVEL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ICORELEVEL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ICORELEVEL','I',ICORELEVEL,RDUM,CDUM,LDUM,CHARAC,N)
! index of ion
      NT_CL=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'CLNT','=','#',';','I', &
     &            NT_CL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''CLNT'' from file INCAR.'
         GOTO 150
      ENDIF
      NT_CL=MIN(MAX(NT_CL,1),NTYP)
      CALL XML_INCAR('CLNT','I',NT_CL,RDUM,CDUM,LDUM,CHARAC,N)

! main quantum number (n)
      N_CL=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'CLN','=','#',';','I', &
     &            N_CL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''CLN'' from file INCAR.'
         GOTO 150
      ENDIF
      N_CL=MIN(MAX(N_CL,0),7)
      CALL XML_INCAR('CLN','I',N_CL,RDUM,CDUM,LDUM,CHARAC,N)

! l quantum number
      L_CL=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'CLL','=','#',';','I', &
     &            L_CL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''CLL'' from file INCAR.'
         GOTO 150
      ENDIF
      L_CL=MAX(MIN(L_CL,N_CL-1),0)
      CALL XML_INCAR('INCAR','I',L_CL,RDUM,CDUM,LDUM,CHARAC,N)
!
      Z_CL=0.0
      CALL RDATAB(LOPEN,'INCAR',IU5,'CLZ','=','#',';','F', &
     &            IDUM,Z_CL,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''CLZ'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('CLZ','F',IDUM,Z_CL,CDUM,LDUM,CHARAC,N)

      ALLOCATE(CL_SHIFT(CL_MAXNL,NIONS))
      RETURN

  150 CONTINUE
      IF (IU0>=0) &
      WRITE(IU0,151) IERR,N
  151 FORMAT(' Error code was IERR=',I1,' ... . Found N=',I5,' data.')
  152 CONTINUE
      STOP

   END SUBROUTINE INIT_CL_SHIFT

!*******************************************************************
!
! write the tags associated with the core level shifts if required
!
!*******************************************************************

    SUBROUTINE XML_WRITE_CL_SHIFT
      IMPLICIT NONE 

      LOGICAL :: LDUM
      INTEGER :: IDUM
      REAL(q) :: RDUM
      COMPLEX(q)  :: CDUM
      CHARACTER*1 :: CHARAC

      CALL XML_INCAR('ICORELEVEL','I',ICORELEVEL,RDUM,CDUM,LDUM,CHARAC,1)
      IF (ICORELEVEL==0) RETURN

      CALL XML_INCAR('CLNT','I',NT_CL,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('CLN','I',N_CL,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('INCAR','I',L_CL,RDUM,CDUM,LDUM,CHARAC,1)
      CALL XML_INCAR('CLZ','F',IDUM,Z_CL,CDUM,LDUM,CHARAC,1)

   END SUBROUTINE XML_WRITE_CL_SHIFT

!*******************************************************************
!
! write information about core level shifts
!
!*******************************************************************

   SUBROUTINE WRITE_CL_SHIFT(IU6)
     
    IMPLICIT NONE
    INTEGER IU6

    IF (ICORELEVEL==0 .OR. IU6<0 ) RETURN
    
    WRITE(IU6,100) ICORELEVEL
    IF (ICORELEVEL >1) THEN
       WRITE(IU6,110) NT_CL, N_CL, L_CL, Z_CL
    ENDIF
100 FORMAT(' Core level calculations are selected ICORELEVEL = ',I2)
110 FORMAT('  CLNT =       ',I5,'       CLN = ',I5,/&
           '  CLL  =       ',I5,'       CLZ = ',F5.3)

   END SUBROUTINE WRITE_CL_SHIFT

!*******************************************************************
!
! function to querry whether accurate core levels are required
!
!*******************************************************************

    FUNCTION  ACCURATE_CORE_LEVEL_SHIFTS()
      IMPLICIT NONE
      LOGICAL ACCURATE_CORE_LEVEL_SHIFTS

      IF (ICORELEVEL>0) THEN
         ACCURATE_CORE_LEVEL_SHIFTS=.TRUE.
      ELSE
         ACCURATE_CORE_LEVEL_SHIFTS=.FALSE.
      ENDIF

    END FUNCTION ACCURATE_CORE_LEVEL_SHIFTS

!*******************************************************************
!
! function to querry whether accurate core levels are required
!
!*******************************************************************

    FUNCTION  FINALST_CORE_LEVEL_SHIFTS()
      IMPLICIT NONE
      LOGICAL FINALST_CORE_LEVEL_SHIFTS

      IF (ICORELEVEL>1) THEN
         FINALST_CORE_LEVEL_SHIFTS=.TRUE.
      ELSE
         FINALST_CORE_LEVEL_SHIFTS=.FALSE.
      ENDIF

    END FUNCTION FINALST_CORE_LEVEL_SHIFTS

!*******************************************************************
!
! this subroutine sets up a local pseudo test charge in the PAW
! sphere with norm 1
!
!*******************************************************************

    SUBROUTINE SET_CL_DRHOCORE_PS(DRHOCORE, NT, R, AUG)
      USE radial
      USE constant
      IMPLICIT NONE
      REAL(q) DRHOCORE(:)    ! charge distribution 
      INTEGER     :: NT      ! index of ion currently treated
      TYPE(rgrid) :: R       ! radial grid descriptor
      REAL(q) :: AUG(:,0:)   ! 1-normalized L-dep compensation charge
      
      ! shortcut if nothing needs to be done
      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS() .OR. NT /= NT_CL) RETURN
      
      ! set augmentation charge
      ! we need to multiply by Y(00) since 
      ! AUG is one normalised 
      DRHOCORE(1:R%NMAX)=AUG(:,0)*Z_CL*(1._q/(2*SQRT(PI)))

    END SUBROUTINE SET_CL_DRHOCORE_PS

!*******************************************************************
!
! this subroutine sets up a local AE test charge in the PAW
! sphere with norm 1
!
!*******************************************************************

    SUBROUTINE SET_CL_DRHOCORE_AE(DRHOCORE, NT,  RHOC, POTVAL , R, ZCORE, ZVAL) 
      USE constant
      USE radial
      IMPLICIT NONE
      REAL(q) DRHOCORE(:)    ! charge distribution 
      INTEGER     :: NT      ! index of ion currently treated
      REAL(q) RHOC(:)        ! electronic core charge distribution rho(r)
      REAL(q) POTVAL(:)      ! valence contribution of the potential for the reference
                             ! state (V_atom)
      TYPE(rgrid) :: R       ! radial grid descriptor
      REAL(q) ZCORE          ! electronic core charge
      REAL(q) ZVAL           ! number of valence electrons
! local variables
      REAL(q) APOT(R%NMAX)   ! total potential
      REAL(q) SCALE, DHARTREE
      LOGICAL :: LINIT=.FALSE.
      REAL(q),ALLOCATABLE,SAVE :: DRHO(:,:),EIG(:)
      INTEGER :: Z           ! atomic number

      ! shortcut if nothing needs to be done
      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS() .OR. NT /= NT_CL) RETURN

      !
      ! special case N_CL==0, just return without modifications of DRHOCORE
      !
      IF (N_CL==0) RETURN

      IF (.NOT.LINIT) THEN
         LINIT=.TRUE.
         ! set hartree potential
         ! thats just Y_00
         SCALE=1/(2*SQRT(PI))

         !     first evaluate the atomic reference inside the PAW
         !     sphere
         !     calculate the Hatree potential of the core electrons
         CALL RAD_POT_HAR(0, R, APOT, RHOC, DHARTREE)
!         WRITE(0,'(8F14.7)') APOT(R%NMAX-7:R%NMAX)* R%R(R%NMAX-7:R%NMAX)/RYTOEV/AUTOA*SCALE
         !     add the potential of the nucleus (delta like charge at origin)
         APOT=APOT*SCALE - FELECT/R%R*(ZCORE+ZVAL)
         !     add reference valence electron potential POTVAL
         !     (this one contains essentially valence only contributions)
         APOT= APOT-POTVAL
         Z = int(ZCORE + ZVAL)

         ALLOCATE(DRHO(1:R%NMAX,1),EIG(1))
         CALL CORE_WAVE_FKT(DRHO, EIG, N_CL, L_CL, APOT, R ,Z, 1)
  
         DRHO=DRHO*SCALE

!         WRITE(94,*) 'eigenvalue',N_CL,L_CL,EIG
!         WRITE(94,'(F14.7)') DRHO(:,1)
      ENDIF

      DRHOCORE(1:R%NMAX)=DRHO(:,1)*Z_CL
         

    END SUBROUTINE SET_CL_DRHOCORE_AE

!*******************************************************************
!
! this subroutine adds the Hartee potential of the test charge
! to the current local potential in the PAW sphere
! to speed up things the type index is passed down
! and checked against  NT_CL
! the resulting potential is more attractive (negative)
!*******************************************************************

    SUBROUTINE ADD_CL_HARTREE_POT(DRHOCORE, NT, ISPIN, POT, R) 
      USE radial
      USE constant
      IMPLICIT NONE
      REAL(q) DRHOCORE(:)    ! charge distribution 
      INTEGER     :: NT      ! index of type currently treated
      TYPE(rgrid) :: R       ! radial grid descriptor
      REAL(q) :: POT(:,:,:)  ! potential
      INTEGER :: ISPIN       ! spin polarised
    ! local
      REAL(q) :: DHARTREE
      REAL(q) APOT(R%NMAX)   ! total potential

!      WRITE(0,'("B",8F14.7)') POT(R%NMAX-7:R%NMAX,1,1)* R%R(R%NMAX-7:R%NMAX)/RYTOEV/AUTOA/(2*SQRT(PI))
      ! shortcut if nothing needs to be done
      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS() .OR. NT /= NT_CL) RETURN
      
      ! set hartree potential
      APOT=0
      CALL RAD_POT_HAR(0, R, APOT, DRHOCORE, DHARTREE)
!      WRITE(0,'("A",8F14.7)') APOT(R%NMAX-7:R%NMAX)* R%R(R%NMAX-7:R%NMAX)/RYTOEV/AUTOA/(2*SQRT(PI))

      POT(1:R%NMAX, 1,1) =POT(1:R%NMAX,1,1)-APOT

      IF (ISPIN==2) POT(1:R%NMAX,1,2)=POT(1:R%NMAX,1,2)-APOT
!      WRITE(0,'("C",8F14.7)') POT(R%NMAX-7:R%NMAX,1,1)* R%R(R%NMAX-7:R%NMAX)/RYTOEV/AUTOA/(2*SQRT(PI))

    END SUBROUTINE ADD_CL_HARTREE_POT


!*******************************************************************
!
! this subroutine calculates the AE wavefunctions for a 
! particular PAW sphere
! as arguments it takes the AE potential in the PAW sphere
! (derived from  RHOC, ZCORE, ZVAL, POTATOM and possibly DPOT)
! and returns the core wavefunctions W
!
!*******************************************************************

    SUBROUTINE SET_CORE_WF( RHOC, POTVAL , R, ZCORE, ZVAL, W, N, L, EIG, DPOT )

      USE prec
      USE constant
      USE radial
      IMPLICIT NONE

      REAL(q) RHOC(:)        ! electronic core charge distribution rho(r)
      REAL(q) POTVAL(:)      ! valence contribution of the potential for the reference
                             ! state (V_atom)
      TYPE(rgrid) :: R       ! radial grid descriptor
      REAL(q) ZCORE          ! electronic core charge
      REAL(q) ZVAL           ! number of valence electrons
      REAL(q) W(:,:)         ! on return: AE wavefunctions phi(r,(n,l)) 
              ! the main and l quantum numbers are storen in n and l      
      INTEGER:: N(:)         ! on return: main quantum number       
      INTEGER:: L(:)         ! on return: l quantum number       
      REAL(q):: EIG(:)       ! on return: eigenenergies
      REAL(q), OPTIONAL :: DPOT(:) 
                             ! difference potential V_actual - V_atom
! local variables
      INTEGER N1,L1          ! main and l quantum number
      REAL(q) APOT(R%NMAX)   ! total potential
      INTEGER :: Z           ! atomic number
      INTEGER :: INDEX	     ! combined main and l quantum number index
      INTEGER :: MAXNL	     ! maximal value for index
      LOGICAL :: NEXTNL
      REAL(q) SUM,SCALE

!     thats just Y_00
      SCALE=1/(2*SQRT(PI))

!     first evaluate the atomic reference inside the PAW
!     sphere
!     calculate the Hatree potential of the core electrons
      CALL RAD_POT_HAR(0, R, APOT, RHOC, SUM)
!     add the potential of the nucleus (delta like charge at origin)
      APOT=APOT*SCALE - FELECT/R%R*(ZCORE+ZVAL)
!      WRITE(94,"(2F14.7)") (R%R(N1), APOT(N1)*R%R(N1)/RYTOEV/AUTOA,N1=1,R%NMAX)
!     add reference valence electron potential POTVAL
!     (this one contains essentially valence only contributions)
      APOT= APOT-POTVAL
      IF (PRESENT(DPOT)) THEN
         APOT=APOT+DPOT(1:R%NMAX)*SCALE
      ENDIF
!      WRITE(94,"(2F14.7)") (R%R(N1), APOT(N1)*R%R(N1)/RYTOEV/AUTOA,N1=1,R%NMAX)

      Z = int(ZCORE + ZVAL)
 
      CALL CALCULATE_MAX_N_L(ZCORE,MAXNL)

      CALL INIT_N_L(INDEX,N1,L1)

      DO
         NEXTNL = NEXT_N_L(MAXNL, INDEX, N1, L1)
         IF (.NOT. NEXTNL) EXIT
         CALL CORE_WAVE_FKT(W, EIG, N1, L1, APOT, R ,Z,INDEX)
         N(INDEX)=N1
         L(INDEX)=L1
      ENDDO

    END SUBROUTINE


!*******************************************************************
!
! calculates the storage requirement for  storing the 
! core wavefunctions of a particular species with a core charge Z
! i.e. how many wavefunctions need to be stored
! 
!*******************************************************************

    SUBROUTINE CALCULATE_MAX_N_L( Z, MAXNL )
      USE prec
      IMPLICIT NONE
      REAL(q) Z              ! number of core electrons
      INTEGER MAXNL          ! on return: storage requirement

      SELECT CASE (INT(Z))
        CASE (:0)
           MAXNL=0
        CASE (1:2)
           MAXNL=1
        CASE (3:4)
           MAXNL=2
        CASE (5:10)
	   MAXNL=3
	CASE (11:12)
	   MAXNL=4
	CASE (13:18)
	   MAXNL=5
	CASE (19:20)
	   MAXNL=6
	CASE (21:30)
	   MAXNL=7
	CASE (31:36)
	  MAXNL=8
	CASE (37:38)
	  MAXNL=9
	CASE (39:48)
	  MAXNL=10
	CASE (49:54)
	  MAXNL=11
	CASE (55:56)
	  MAXNL=12
	CASE (57:70)
	  MAXNL=13
	CASE (71:80)
	  MAXNL=14
	CASE (81:86)
	  MAXNL=15
	CASE (87:88)
	  MAXNL=16
	CASE (89:102)
	  MAXNL=17
	CASE (103:112)
	  MAXNL=18
      END SELECT

    END SUBROUTINE

!*******************************************************************
!
! increase the main (n) and l quantum number (l)
! 
! functions returns false when all indexes have been exhaustet
! for the particular MAXNL, where MAXNL must be determined with
! the subroutine CALCULATE_MAX_N_L
! 
!*******************************************************************

      ! initialise the counters properly for the routine NEXT_N_L

    SUBROUTINE INIT_N_L( INDEX, N, L)
      IMPLICIT NONE
      INTEGER INDEX, N, L

      INDEX=0
      N=-1
      L=-1
      
    END SUBROUTINE


    FUNCTION NEXT_N_L( MAXNL, INDEX, N, L)
      USE prec
      IMPLICIT NONE
      LOGICAL  NEXT_N_L
      REAL(q) Z
      INTEGER MAXNL     ! maximum index
      INTEGER INDEX     ! index number
      INTEGER N,L       ! main n and l quantum number

      INDEX=INDEX+1
      IF (INDEX > MAXNL) THEN
        NEXT_N_L=.FALSE.
      ELSE
        NEXT_N_L=.TRUE.
        SELECT CASE (INDEX)
          CASE (1)
             N=1 ; L=0
          CASE (2)
             N=2 ; L=0
          CASE (3)
	     N=2 ; L=1
	  CASE (4)
	     N=3 ; L=0
	  CASE (5)
	     N=3 ; L=1
	  CASE (6)
	     N=4 ; L=0
	  CASE (7)
	     N=3 ; L=2
	  CASE (8)
	     N=4 ; L=1
	  CASE (9)
	     N=5 ; L=0
	  CASE (10)
	     N=4 ; L=2
	  CASE (11)
	     N=5 ; L=1
	  CASE (12)
	     N=6 ; L=0
	  CASE (13)
	     N=4 ; L=3
	  CASE (14)
	     N=5 ; L=2
	  CASE (15)
	     N=6 ; L=1
	  CASE (16)
 	     N=7 ; L=0
	  CASE (17)
	     N=5 ; L=3
	  CASE (18)
	     N=6 ; L=2
        END SELECT
      ENDIF
      
    END FUNCTION NEXT_N_L

!***********************************************************************
!
!	Berechnet Wellenfunktionen und Eigenenergien                  
!
!***********************************************************************

  SUBROUTINE CORE_WAVE_FKT(W, EIG, N, L, APOT, R ,Z,INDEX)

    USE radial
    USE constant
    TYPE (rgrid) :: R		! Grid
    REAL (q) :: W(:,:)		! Wellenfkt
    REAL (q) :: EIG(:)		! Eigenenergie
    INTEGER :: N		! Hauptquantenzahl
    INTEGER :: L		! Nebenquantenzahl
    INTEGER :: Z		! Ordnungszahl
    INTEGER :: FJ		! relat. quantenzahl
    INTEGER :: NTIMES		! zaehler
    INTEGER :: KI		! Endpunkt auswaertsintegration, aus OUTINT,
    INTEGER :: KJ               ! Startpunkt Einwaertsintegration, aus INWINT
    INTEGER :: NQNT		! Anzahl der nodes, die diewellefkt haben soll
    INTEGER :: NODES,INODES	! Anzahl der Knoten, aus OUTINT
    INTEGER :: NERR		! Fehler in MILNE-Procedure
    INTEGER :: INDEX  		! Index von Wellenfkt und Eigenenergie
    REAL (q):: APOT(R%NMAX) 	! Potential
    REAL (q):: EMIN, EMAX		! Energie unter- und obergrenze
    REAL (q):: E			! Energie in Integrationsvorgang
    REAL (q):: A(R%NMAX),B(R%NMAX)	! Wellenfkt fuer bestimmten Index
    REAL (q):: H3,RA,WN,RJ,RG      ! fuer Normierung
    INTEGER :: K,K1                     ! Zahler fuer Normierung
    LOGICAL :: DBL                    ! Hilfsvariable Normierung
    REAL(q), PARAMETER :: C = 137.03  ! Lichtgeschwindigkeit
    LOGICAL :: finish                 ! true, wenn E genau genug
    REAL(q) :: MA,MB         ! Quotient von A (bzw B) aus INWINT und OUTINT
    REAL(q) :: DG,QQQQ,QCOEF,WMIN,DELTE0,DE

1   FORMAT (//10X,"TOO MANY TRIES REQUIRED  ")
2   FORMAT (//10X,"TOO MANY NODES.  ",I5,3F4.1,3E15.7)
3   FORMAT (10X,"DE  TOO  LARGE.",3F4.1,4E15.7)
4   FORMAT (//100X,"TOO FEW NODES.   ",I5,3F4.1,3E15.7)
5   FORMAT (10X,5E14.4)
6   FORMAT (10X,"HARD",I3,"ERRORS OCCURED IN MILNE PROCEDURE")
7   FORMAT (10X,"mesh",8X,"R",13X,"V*r",7X,"A     ,6X,    ","B    ")
8   FORMAT (//100X,"CORRECT NODES.   ",I5,3F4.1,3E15.7)


    FJ = 2 * L
    EMIN = -Z*Z*RYTOEV*2/(N)**2
    EMAX = -1.0E-5
    E = 0.5 * (EMAX+EMIN)
    NQNT = N-L-1
    NERR=0
    H3=R%H/3.
    DELTE0=0
    finish=.false.

    DO NTIMES=1,200
      NERR=0
      CALL OUTINT (E,L,KI,NODES,R,APOT,A,B,NERR)	

      IF (NODES < NQNT) THEN		! too few nodes
        IF (E > EMIN) EMIN=E
        E = 0.5 * (E+EMAX)
        IF ((EMAX-EMIN) < (E * 1E-6)) THEN
          WRITE(*,4) NODES,N,L,FJ,EMIN,E,EMAX
	  STOP
	ENDIF
      
      ELSEIF (NODES > NQNT) THEN	! too many nodes
	IF (E < EMAX) EMAX=E
	E = 0.5 * (E+EMIN)
	IF ((EMAX-EMIN) < (E * 1E-6)) THEN	! warum -E ???
	  WRITE(*,2) NODES,N,L,FJ,EMIN,E,EMAX
	  STOP
	ENDIF
      
      ELSE				! correct Nr of nodes

        MA=A(KI)
        MB=B(KI)
        CALL INWINT (E,L,KI,KJ,R,APOT,A,B,INODES)
        MA=MA/A(KI)
        MB=MB/B(KI)

        IF (INODES /= 0) WRITE(*,*)"Nodes in INWINT!",INODES
        DO K=KI,KJ                      ! macht A stetig bei KI
          A(K)=A(K)*MA
          B(K)=B(K)*MA
        ENDDO

        RG=A(KJ)
        DG=B(KJ)
        QQQQ=DMIN1(ABS(RG),ABS(DG))
        IF(QQQQ < 1.0E-25_q) THEN       ! setzt A=0 ausserhalb KJ
          DO K=KJ+1,R%NMAX
            A(K)=0.0
            B(K)=0.0
          ENDDO
        ELSE
          WMIN=(1.E-35_q)/QQQQ          ! exponentieller Abfall
          QCOEF=SQRT(-(E+E))            ! ausserhalb KJ
          RJ=R%RSTART*R%D**(KJ-1)
          RA=RJ
          DO K=KJ+1,R%NMAX
            RA=RA*R%D
            WN=EXP(QCOEF*(RJ-RA))
            IF (WN < WMIN) THEN
              DO K1=K,R%NMAX
                A(K1)=0.0
                B(K1)=0.0
              ENDDO
              EXIT
            ENDIF
            A(K)=WN*RG
            B(K)=WN*DG
          ENDDO
        ENDIF

                                        ! Normierung von A und B
        RA=R%REND
        WN=RA*(A(R%NMAX)**2+(B(R%NMAX)/C)**2)/2
        K=R%NMAX
        RA=RA+RA
        RJ=1.0_q/R%D
        DBL=.FALSE.
        DO K=R%NMAX-1,2,-1
          RA=RA*RJ
          DBL=.NOT.DBL
          RG=RA*(A(K)*A(K)+(B(K)/C)**2)
          WN=WN+RG
          IF (DBL) WN=WN+RG
        ENDDO    ! IF (K.GT.2)  GO TO 32
        WN=H3*(WN+R%RSTART*(A(1)**2+(B(1)/C)**2)/2)
        
        DE=(MA-MB)/MA
        DE=((-(B(KI)*A(KI)))/WN)*DE
        IF (DE*DELTE0 < 0) DE=DE*0.5    ! bei Vorzeichenwechsel DE/2
        DELTE0=DE
        IF ((E > -0.001) .AND. (DE > 0.01)) finish=.true.
        IF ((EMAX-EMIN) < (-E*1E-7)) finish=.true.
        IF (DE > 0) EMIN=E
        IF (DE < 0) EMAX=E
        IF (.NOT. finish) E=E+DE
        IF (ABS(DE) < (-E*1E-6)) finish=.true.
        IF ((E >= EMAX) .OR. (E <= EMIN)) E=0.5*(EMAX+EMIN)
        
        IF (finish) THEN
          WN=1./SQRT(WN)
          A=A*WN
          B=B*WN/C    
          IF (NERR /= 0) WRITE(*,6) NERR 
          W(:,INDEX)=A*A+B*B
          EIG(INDEX)=E
          RETURN
        ENDIF
      ENDIF

    ENDDO

    WRITE(*,1)				! zu viele Durchlaeufe
    WRITE(*,5) N,L,FJ,EMIN,EMAX
    WRITE(*,1)
    WRITE(*,7)
    STOP    

  END SUBROUTINE


!*****************************************************************
!
!	Performs outward integration. The Runge-Kutte procedure
!	is used to provide the starting values for the MILNE
!	integration.
!
!******************************************************************

  SUBROUTINE OUTINT (E,L,KI,NODES,R,APOT,A,B,NERR)
	
    USE constant
    USE radial
    TYPE(rgrid) :: R
    INTEGER:: L
    REAL(q):: DX2
    REAL(q):: RA
    REAL (q) :: APOT(R%NMAX),VR(R%NMAX)         ! Potential, P. in atom units
    REAL(q):: A(R%NMAX),B(R%NMAX)	        ! Wellenfunktionen     
    INTEGER :: NODES				! Anzahl der Knoten
    REAL(q):: HOC,XK,P,XC,BGC,TR,TC,VC,UC,WC	! Hilfsvariable
    REAL(q):: X,UNP,WNP,UNP2,WNP2,XMFT,TEST     !  HILFSVARIABLE MILNE
    INTEGER :: KI,IK,K,KP1			!  HILFSVARIABLE
    INTEGER :: KM1,KM2,KM3,KM4,KM5,KIT		! HILFSVARIABLE MILNE
    INTEGER :: NERR				! Fehler in MILNE-Procedure
    REAL(q):: DA(5),DB(5)			! zur Berechnung der Wellenfkt
						! in Runge-Kutte Routine
    REAL(q):: AP(R%NMAX),BP(R%NMAX)	        ! Zwischenspeicher in
						! MILNE-Routine
    REAL(q), PARAMETER :: C = 137.037		! Lichtgeschwindigkeit
    REAL(q):: E,EAU				! Energie, E. in atom units
    REAL(q):: RN,RNOT				! R%REND, R%RSTART in atomunits
    LOGICAL :: INIT=.FALSE.

 ! convert to atomic units
    EAU = E / (2*RYTOEV)
    VR = APOT / (2*RYTOEV) *R%R /AUTOA

    RNOT=R%RSTART/AUTOA
    RN=R%REND/AUTOA

    DX2 = R%H / 2.
    XK=L*(L+1)
    XMFT=4.4444444444E-2_q * R%H
    TEST=1.0E6_q
 
    NODES=0

 !..... FIND THE CLASSICAL TURNING POINT........
    RA = RN
    DO KI=R%NMAX-1,11,-1
      RA=RA/R%D
      IF ((EAU*RA)>=VR(KI))  EXIT
    ENDDO
    KI=KI+1
    IF ((KI+10)>=R%NMAX) KI=R%NMAX-11

			!    	set up starting values
    IF (ABS(VR(1))>=0.1) THEN
      A(1)=1.0E-10_q
      HOC=-VR(1)/C
      P=(SQRT(XK-HOC*HOC+1)-1)/HOC
      IF ((ABS(HOC)<0.001) .AND. (L <= 0.01)) P=-HOC/2
      B(1)=C*P*A(1)
    ELSE
      A(1)=1.0E-10_q
      B(1)=L/2/RNOT*A(1)
    ENDIF

			!  Starting Runge-Kutte procedure
    X=LOG(RNOT)
    DO K=1,5
      KP1=K+1
      XC=X
      BGC=-VR(K)
      WC=B(K)
      UC=A(K)
      DO IK=1,4
        RA=EXP(XC)
        TR=2*RA
        TC=EAU*RA+BGC                        !
        VC=(1/C/C)*TC
        DA(IK)=DX2*((VC+TR)*WC+UC)
        DB(IK)=DX2*((-1.)*WC-(TC-XK/(VC+TR))*UC)
        SELECT CASE(IK)
           CASE(1)
	     XC=XC+DX2
             UC=UC+DA(1)
	     WC=WC+DB(1)
             BGC=0.5_q*(BGC-VR(KP1))
	   CASE(2)
	     UC=UC+DA(2)-DA(1)
	     WC=WC+DB(2)-DB(1)
	   CASE(3)
	     XC=XC+DX2
	     UC=UC+2*DA(3)-DA(2)
	     WC=WC+2*DB(3)-DB(2)
	     BGC=-VR(KP1)
	END SELECT
      ENDDO
      IK=4
      B(KP1)=B(K)+(DB(1)+DB(4)+2*(DB(2)+DB(3)))*0.33333333333333_q
      A(KP1)=A(K)+(DA(1)+DA(4)+2*(DA(2)+DA(3)))*0.33333333333333_q
      AP(KP1)=B(KP1)*(VC+TR)+A(KP1)
      BP(KP1)=-B(KP1)-(TC-XK/(VC+TR))*A(KP1)
      X=X+R%H
    ENDDO
			!  Begin Milne procedure
    RA=EXP(X)
    DO K=6,KI-1
      RA=RA*R%D
      KP1=K+1
      KM1=K-1
      KM2=K-2
      KM3=K-3
      KM4=K-4
      KM5=K-5
      TR=2*RA
      TC=EAU*RA-VR(KP1)
      VC=1/C/C*TC
      UNP=A(KM5)+0.3_q*R%H*(11*(AP(K)+AP(KM4))+26*AP(KM2)-14*(AP(KM1)+AP(KM3)))
      WNP=B(KM5)+0.3_q*R%H*(11*(BP(K)+BP(KM4))+26*BP(KM2)-14*(BP(KM1)+BP(KM3)))
      KIT=0
      
33    AP(KP1)=(VC+TR)*WNP+UNP
      BP(KP1)=(-1)*WNP-(TC-XK/(VC+TR))*UNP
      UNP2=A(KM3)+(7*(AP(KP1)+AP(KM3))+32*(AP(KM2)+AP(K))+12*AP(KM1))*XMFT
      WNP2=B(KM3)+(7*(BP(KP1)+BP(KM3))+32*(BP(KM2)+BP(K))+12*BP(KM1))*XMFT

      IF((ABS(TEST*(UNP2-UNP)) > ABS(UNP2)) &
      & .OR. (ABS(TEST*(WNP2-WNP)) > ABS(WNP2))) THEN
        IF (KIT < 10) THEN
	  KIT=KIT+1
          WNP=WNP2
          UNP=UNP2
	  GOTO 33
 	ELSE IF (NERR >= 0) THEN		!  nerr	
	  NERR=NERR+1
        ELSE
	  WRITE(*,83) E,K,UNP2,UNP,WNP2,WNP
83        FORMAT(1H,10X,"HARD TEST, E= ",E20.8,I5,4E20.8)
	ENDIF
      ENDIF

      B(KP1)=WNP2
      A(KP1)=UNP2

      NODES=NODES+0.501*ABS(SIGN(1._q,A(K))-SIGN(1._q,A(K-1)))
    ENDDO
    RETURN

  END SUBROUTINE



!************************************************************************
!*     PERFORM THE  INWARD INTEGRATION.   THIS ROUTINE IS A
!*     DERIVATIVE OF START2/DIFF2.   MUCH UNNECESSARY INDEXING
!*     HAS BEEN REMOVED AMONG OTHER THINGS.    DALE KOELLING
!*     MODIFICATION TO FOLLY13 EQUATIONS...DECEMBER 1, 1976...BNH.
!*************************************************************************


  SUBROUTINE INWINT (E,FL,KI,KJ,R,APOT,A,B,NODES)

    USE constant
    USE radial
 
    IMPLICIT NONE

    TYPE (rgrid) :: R                  ! Grid
    REAL(q) :: E,EAU                   ! Energie, E in atom units
    INTEGER :: FL                      ! Nebenquantenzahl
    INTEGER :: KJ,KI                   ! Start-, Endpunkt der Int.
    REAL(q) :: A(R%NMAX),B(R%NMAX)     ! Wellenfkt
    REAL(q) :: APOT(R%NMAX),VR(R%NMAX) ! Potential, P. in atom units

    REAL(q) :: DA(5),DB(5)     
    REAL(q),PARAMETER :: EPS=75.       ! VALUE USED TO DETERMINE THE PRACTICAL INFINITY
    REAL(q) :: RN,RNOT	               ! R%REND, R%RSTART in atomunits
    REAL(q) :: RA
    INTEGER :: M,K,I                   ! Hilfsvariablen
    REAL(q) :: DR                      ! Gridvariable
    REAL(q),PARAMETER :: C = 137.037   ! Lichtgeschwindigkeit
    REAL(q) :: H3                      ! Hilfsvariable
    REAL(q) :: RP12,RP21,ATK,BTK
    INTEGER :: NODES

    REAL(q),PARAMETER :: F1=1.045833333,F2=2.691666666,F3=1.1, &
           F4=0.4416666666, F5=0.07916666666 
    REAL(q),PARAMETER :: G1=0.9666666666,G2=4.133333333,G3=0.8, &
           G4=0.1333333333,G5=0.03333333333
    REAL(q),PARAMETER :: FG1=0.9333333333,FG2=4.266666666,FG3=1.6

    H3=-R%H/3.
    DR=1./R%D

  ! convert to atomic units
    EAU = E / (2*RYTOEV)
    VR = APOT / (2*RYTOEV) *R%R /AUTOA
    RNOT=R%RSTART/AUTOA
    RN=R%REND/AUTOA

    RA=RNOT*(R%D**(KI+10))       

    DO KJ=KI+11,R%NMAX
      IF(RA*(VR(KJ)-EAU*RA) > EPS) EXIT
      RA=RA*R%D
    ENDDO

    IF (KJ>R%NMAX) THEN
      KJ=R%NMAX 
      B(KJ)=1.0 ! MOD(L,2)
      A(KJ)=1.0_q-B(KJ)
    ELSE
      A(KJ)=1.0_q
      B(KJ)=SQRT(-EAU/(2.0_q*C*C+EAU))
    ENDIF

    DO M=1,4
      K=KJ-M
      A(K)=A(KJ)
      B(K)=B(KJ)
    ENDDO

    DO I=1,6
      K=KJ+1
      RA=RN*DR**(R%NMAX-K)
      DO M=1,5
        K=K-1
        RA=RA*DR
        RP12=RA+RA+(EAU*RA-VR(K))/C/C
        RP21=VR(K)-EAU*RA + FL*(FL+1.)/RP12
        DA(M)=H3*(1.*A(K)+RP12*B(K))
        DB(M)=H3*((-1.)*B(K)+RP21*A(K))
      ENDDO
      M=KJ-1
      A(M)=A(KJ)+F1*DA(1)+F2*DA(2)+F4*DA(4)-(F3*DA(3)+F5*DA(5))
      B(M)=B(KJ)+F1*DB(1)+F2*DB(2)+F4*DB(4)-(F3*DB(3)+F5*DB(5))
      M=M-1
      A(M)=A(KJ)+G1*DA(1)+G2*DA(2)+G3*DA(3)+G4*DA(4)-G5*DA(5)
      B(M)=B(KJ)+G1*DB(1)+G2*DB(2)+G3*DB(3)+G4*DB(4)-G5*DB(5)
      M=M-1
      A(M)=A(KJ)+1.0125_q*DA(1)+3.825_q*DA(2)+2.7_q*DA(3)+1.575_q*DA(4) &
               -0.1125_q*DA(5)
      B(M)=B(KJ)+1.0125_q*DB(1)+3.825_q*DB(2)+2.7_q*DB(3)+1.575_q*DB(4) &
               -0.1125_q*DB(5)
      M=M-1
      A(M)=A(KJ)+FG1*DA(1)+FG2*DA(2)+FG3*DA(3)+FG2*DA(4)+FG1*DA(5)
      B(M)=B(KJ)+FG1*DB(1)+FG2*DB(2)+FG3*DB(3)+FG2*DB(4)+FG1*DB(5)
    ENDDO 

    DA(1)=DA(2)
    DA(2)=DA(3)
    DA(3)=DA(4)
    DB(1)=DB(2)
    DB(2)=DB(3)
    DB(3)=DB(4)
    K=KJ-3
    RA=RN*DR**(R%NMAX-K)

    DO K=KJ-4,KI,-1
      RA=RA*DR
      RP12=RA+RA+(EAU*RA-VR(K))/C/C
      RP21=VR(K)-EAU*RA+FL*(FL+1.)/RP12
      ATK=A(K+4)+8.0_q*(DA(3)+DA(1)-0.5_q*DA(2))
      BTK=B(K+4)+8.0_q*(DB(3)+DB(1)-0.5_q*DB(2))
      DA(4)=H3*(1.*ATK+RP12*BTK)
      DB(4)=H3*((-1.)*BTK+RP21*ATK)
      A(K)=A(K+1)+1.125_q*DA(4)+2.375_q*DA(3)-0.625_q*DA(2)+0.125_q*DA(1)
      B(K)=B(K+1)+1.125_q*DB(4)+2.375_q*DB(3)-0.625_q*DB(2)+0.125_q*DB(1)
      DA(1)=DA(2)
      DA(2)=DA(3)
      DB(1)=DB(2)
      DB(2)=DB(3)
      DA(3) =H3*(1.*A(K)+RP12*B(K))
      DB(3) =H3*((-1.)*B(K)+RP21*A(K))      
    ENDDO

    NODES=0
    DO M=KI,KJ
      NODES=NODES+0.501*ABS(SIGN(1._q,A(K))-SIGN(1._q,A(K-1)))
    ENDDO

    RETURN
  END SUBROUTINE


!*******************************************************************
!
!  RAD_CL_SHIFT
!  calculate the core level shift on the radial PAW grid
!  this version uses a localised test charge stored in AUG
!  the full quantum mechanical test charge is applied in the
!  routine RAD_CL_SHIFT_AE
!
!*******************************************************************

    SUBROUTINE RAD_CL_SHIFT( POT, POTAE, R, CL_SHIFT, AUG )
      USE constant
      USE radial
      IMPLICIT NONE

      REAL(q) :: POT(:,:),POTAE(:,:)
      TYPE (rgrid) R
      REAL(q) :: CL_SHIFT
      REAL(q) :: AUG(:)        ! 1-normalized compensation charge

   ! local variables
      INTEGER LMMAIN,I
      REAL(q) SUM,SUM_AE
!-----------------------------------------------------------------------
! calculate \int V(r) aug(r)
!-----------------------------------------------------------------------
      LMMAIN=1
      SUM=0
      SUM_AE=0
      DO I=1,R%NMAX
         SUM   =SUM   +POT  (I,LMMAIN)*AUG(I)*R%SI(I)
         SUM_AE=SUM_AE+POTAE(I,LMMAIN)*AUG(I)*R%SI(I)
      ENDDO
      CL_SHIFT=CL_SHIFT+(SUM_AE-SUM)

    END SUBROUTINE RAD_CL_SHIFT

!*******************************************************************
!
!  RAD_CL_SHIFT_AE
!  calculate the core level shift on the radial PAW grid
!  AUG    is the test charge used on the plane wave grid
!          int POT(r) AUG(r)   is subtracted 
!  POTAE  is optional and specifies the change of the all electron 
!         potential w.r.t an atom
!          int POTAE(r) wavefunction(r) is added
!         if POTAE is supplied to the subroutine
!
!*******************************************************************

    SUBROUTINE RAD_CL_SHIFT_AE( POT, R, CL_SHIFT, AUG, MAXNL, W, EIG, POTAE )
      USE constant
      USE radial
      IMPLICIT NONE

      REAL(q) :: POT(:)      ! local reconstructed pseudopotential
      REAL(q), OPTIONAL :: POTAE(:)
      TYPE (rgrid) R         ! grid descriptor
      REAL(q) :: CL_SHIFT(:) ! on exit: core level eigenvalues 
       INTEGER :: MAXNL      ! number of (n,l) pairs
      REAL(q) :: AUG(:)      ! 1-normalized compensation charge
      REAL(q) :: W(:,:)      ! all electron wave functions
      REAL(q) :: EIG(:)      ! eigenvalues obtained by solving the radial Schroedinger equation
 

   ! local variables
      INTEGER MMAIN,I,LM
      REAL(q) SUM,SUM_AE
      REAL(q) TEST,NORM

      MMAIN=1
      SUM=0
      SUM_AE=0
      !
      ! int V_ps (r) rho_ps (r) r^2 dr
      !
      DO I=1,R%NMAX
         SUM   =SUM   +POT  (I)*AUG(I)*R%SI(I)
      ENDDO
      !
      ! int V_AE (r) rho_AE_nl (r) r^2 dr
      !
      IF (PRESENT(POTAE)) THEN
         DO LM=1,MAXNL
            SUM_AE=0
            NORM=0
            DO I=1,R%NMAX
               SUM_AE=SUM_AE+POTAE(I)*W(I,LM)*R%SI(I)
               NORM=NORM+W(I,LM)*R%SI(I)
            ENDDO
            CL_SHIFT(LM)=CL_SHIFT(LM)+(SUM_AE-SUM)+EIG(LM)
         ENDDO
      ELSE
         DO LM=1,MAXNL
            CL_SHIFT(LM)=CL_SHIFT(LM)-SUM+EIG(LM)
         ENDDO
      ENDIF
    END SUBROUTINE RAD_CL_SHIFT_AE



!*******************************************************************
!
! this test routine can be use the set up the array AUG
! which will afterwards hold the augmentation charges 
!
!*******************************************************************

      SUBROUTINE SETAUG_CL( AUG , NT, P )
      USE prec
      USE constant
      USE radial
      USE pseudo
      IMPLICIT NONE

      REAL(q)              :: AUG(:)
      TYPE (potcar),TARGET :: P
! local variable
      INTEGER, PARAMETER :: NQ=2
      REAL(q)  QQ(NQ)    ! parameters of Besselfunctions
      REAL(q)  A(NQ)     ! parameters of Besselfunctions
      TYPE (rgrid), TARGET  :: R
      REAL(q)  RCORE
      REAL(q) RR,SUM,QR,BJ,INT
      INTEGER I,N,NT

      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS() .OR. NT /= NT_CL) RETURN

#ifdef RCORE_MAIN
      RCORE= RCORE_MAIN
#else
!      WRITE(0,*)'internal error RCORE_MAIN undefined'
!     STOP
      RCORE=N_CL/4.0
#endif
      WRITE(0,*) 'RCORE is set to ',RCORE
      R%NMAX  =400
      R%RSTART=1E-2
      R%REND  =RCORE
      R%D     =(R%REND/R%RSTART)**(1._q/(R%NMAX-1)) 
      R%H     =LOG(R%D)

      ALLOCATE(R%R(R%NMAX))
         
      DO I=1,R%NMAX
         R%R(I)=R%RSTART*EXP((I-1)*R%H)
      ENDDO

      R%RMAX=  R%R(R%NMAX)
      R%REND  =RCORE
         
      CALL SET_SIMP(R)

      CALL AUG_SETQ(0,R,QQ,A,LCOMPAT=.FALSE.)
      AUG=0

      INT=0
      DO N=1,P%R%NMAX
         SUM=0
         IF (P%R%R(N) <= R%RMAX) THEN
            DO I=1,NQ
               QR=QQ(I)*P%R%R(N)
               CALL SBESSEL( QR, BJ, 0)
               SUM=SUM+BJ*A(I)*P%R%R(N)*P%R%R(N)
            ENDDO
         ENDIF
         AUG(N)=SUM
         INT=INT+SUM*P%R%SI(N)
      ENDDO
      AUG=AUG*Z_CL*(1._q/(2*SQRT(PI)))

      DEALLOCATE(R%R)

    END SUBROUTINE SETAUG_CL

!*******************************************************************
!
! this subroutine performs an FFT of the augmentation charges
! to reciprocal space
!    rho(q) =  4 pi /q \int sin(qr) rho(r) r dr
!
!
!*******************************************************************

      SUBROUTINE AUGTOQ( RHOQ, P, ENAUG )
      USE prec
      USE constant
      USE radial
      USE pseudo
      IMPLICIT NONE
      
      REAL(q) RHOQ(:)    ! on return: rho(q)
      TYPE (potcar),TARGET :: P
      REAL(q) ENAUG      ! cutoff for augmentation charge
! local variable
      INTEGER, PARAMETER :: NFFT=32768
      INTEGER NQL        ! number of grid points in local potential          
      REAL(q) DELQL      ! step width for the array RHOQ
      REAL(q) RMAXF      ! maximum r values after FFT
      REAL(q) QMAX       ! maximum q value for FFT
      REAL(q) RDEL       ! spacing in real space
      TYPE (rgrid), TARGET  :: R
      TYPE (rgrid), POINTER :: RP
      INTEGER, PARAMETER :: NQ=2
      REAL(q)  QQ(NQ)    ! parameters of Besselfunctions
      REAL(q)  A(NQ)     ! parameters of Besselfunctions
      REAL(q)  RCORE
      REAL(q) FFT(NFFT)  ! array for FFT
      REAL(q) RR,SUM,QR,BJ
      INTEGER I,K

!=======================================================================
! first set up the radial grid
! (either PAW or a new test charge grid)
!=======================================================================
   ! no PAW potential, chose a proper test charge
   ! and set a temporary grid
#ifndef RCORE_MAIN
      IF (.NOT. ASSOCIATED(P%QPAW)) THEN
         RCORE= 6/SQRT(ENAUG /RYTOEV)
#else
         RCORE= RCORE_MAIN
#endif
         R%NMAX  =400
         R%RSTART=1E-2
         R%REND  =RCORE
         R%D     =(R%REND/R%RSTART)**(1._q/(R%NMAX-1)) 
         R%H     =LOG(R%D)

         ALLOCATE(R%R(R%NMAX))
         
         DO I=1,R%NMAX
            R%R(I)=R%RSTART*EXP((I-1)*R%H)
         ENDDO

         R%RMAX=  R%R(R%NMAX)
         R%REND  =RCORE
         
         CALL SET_SIMP(R)
         RP=>R

#ifndef RCORE_MAIN
      ELSE

   ! PAW case, simply use PAW grid
         RCORE =P%PSDMAX
         RP=>   P%R
      ENDIF
#endif

      ! set up radial grid
      CALL AUG_SETQ(0,RP,QQ,A,LCOMPAT=.FALSE.)
      A=A/PI/4

      IF (.NOT. ASSOCIATED(P%QPAW)) THEN
         DEALLOCATE(R%R)
      ENDIF
!=======================================================================
! now set up test charge and perform FFT of the test charge 
! to reciprocal space
!=======================================================================
      NQL  = NPSPTS
      DELQL= P%PSGMAX/NQL
      RMAXF= 2*PI/DELQL
      RDEL = RMAXF/NFFT

      QMAX = DELQL*NFFT
! first set up the values in the real space array
      DO I=1,NFFT
         RR=RDEL*(I-1)
         IF (RR>RP%RMAX) THEN
            FFT(I)=0
         ELSE
            SUM=0
            DO K=1,NQ
               QR=QQ(K)*RR
               CALL SBESSEL( QR, BJ, 0)
               SUM=SUM+BJ*A(K)
            ENDDO
            FFT(I)=SUM*RR*4*PI
         ENDIF
      ENDDO

      CALL REALFT(FFT,NFFT/2,1)

      DO I=2,NQL
         RHOQ(I)=FFT(I*2)*RDEL/((I-1)*DELQL)
      ENDDO
      RHOQ(1)=1.0_q

    END SUBROUTINE AUGTOQ

    SUBROUTINE CL_MODIFY_PP( NTYP, P, ENAUG )
      USE prec
      USE pseudo
      USE constant
      IMPLICIT NONE
      
      INTEGER NTYP       ! number of types
      TYPE (potcar) :: P(NTYP)
      REAL(q) ENAUG      ! cutoff for augmentation charge
  !  local 
      INTEGER NQL        ! number of grid points in local potential          
      REAL(q) DELQL      ! step width for the array RHOQ
      REAL(q) QQ,FAKT
      REAL(q), ALLOCATABLE ::  RHOQ(:)  ! on return: rho(q)
      INTEGER I
      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS()) RETURN
      NQL   = NPSPTS
      ALLOCATE(RHOQ(NQL))
      
      CALL AUGTOQ( RHOQ, P(NT_CL), ENAUG )
#ifndef new
      RETURN
#endif

      DELQL= P(NT_CL)%PSGMAX/NQL
  ! potential is 4 pi e^2 / q^2 
      DO I=2,NQL 
         QQ=(I-1)*DELQL
         P(NT_CL)%PSP (I,2) =P(NT_CL)%PSP (I,2) - Z_CL*(RHOQ(I) - 1.0)* 4*PI*FELECT/QQ/QQ
      ENDDO
  ! this has to be worked out correctly 
  ! currently I use the potential at grid point 2 and add it to grid point 1
      I=2
      QQ=(I-1)*DELQL
      P(NT_CL)%PSP (1,2) =P(NT_CL)%PSP (1,2) -  Z_CL*(RHOQ(I) - 1.0)* 4*PI*FELECT/QQ/QQ

      P(NT_CL)%PSCORE=P(NT_CL)%PSP(1,2)

  ! refit spline
      CALL SPLCOF(P(NT_CL)%PSP(1,1) ,NPSPTS,NPSPTS,0._q)
      FAKT =  (P(NT_CL)%ZVALF + Z_CL) /  P(NT_CL)%ZVALF
      IF (FAKT==0) FAKT = 1
      P(NT_CL)%PSPRHO =  P(NT_CL)%PSPRHO* FAKT

      P(NT_CL)%ZVALF = P(NT_CL)%ZVALF + Z_CL

    END SUBROUTINE CL_MODIFY_PP

END MODULE

!*******************************************************************
!
! add core level potential to a local potential
!
!*******************************************************************

    SUBROUTINE CL_TEST_POT(CVTOT, CHTOT, CVION, GRIDC, E_CL_ION)
      USE mgrid
      USE cl
      IMPLICIT NONE
      TYPE (grid_3d) :: GRIDC    ! grid descriptor
      COMPLEX(q)     :: CVTOT(GRIDC%RC%NP)
      COMPLEX(q)     :: CHTOT(GRIDC%RC%NP)
      COMPLEX(q)     :: CVION(GRIDC%RC%NP)
      REAL(q)        :: E_CL_ION ! ion-test charge interaction
      REAL(q)        :: E_CL_ION_ ! ion-test charge interaction
      REAL(q)        :: E_CL_ION__ ! ion-test charge interaction

      E_CL_ION=0      
      IF (.NOT. FINALST_CORE_LEVEL_SHIFTS()) RETURN

#ifdef new
      RETURN
#endif

      IF (.NOT. ALLOCATED(CVTEST) .OR.SIZE(CVTEST) /= GRIDC%RC%NP ) THEN
         WRITE(0,*) "internal error: CL_TEST_POT not allocated"
         STOP
      ENDIF
      CALL CL_ION_TEST_CHARGE( CHTEST, CVION, GRIDC , E_CL_ION)
      CALL CL_ION_TEST_CHARGE( CHTEST, CVTOT, GRIDC , E_CL_ION_)
      CALL CL_ION_TEST_CHARGE( CHTOT, CVTEST, GRIDC , E_CL_ION__)

      CVTOT(1:GRIDC%RC%NP)=CVTOT(1:GRIDC%RC%NP)-CVTEST(1:GRIDC%RC%NP)

    END SUBROUTINE CL_TEST_POT

    SUBROUTINE CL_ION_TEST_CHARGE( CHTEST, CVION, GRIDC, E)
      USE prec
      USE mgrid
      IMPLICIT NONE
      TYPE (grid_3d) :: GRIDC    ! grid descriptor
      COMPLEX(q)     :: CHTEST(GRIDC%RC%NP)
      COMPLEX(q)     :: CVION(GRIDC%RC%NP)
      REAL(q)        :: E
    ! local
      INTEGER        :: NG,NC,N1,N2,N3
      REAL(q)        :: FACTM

      E=0
      NG=1
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        SETFACT1
        SETFACT

        E=E- MULFACT CVION(NG)*CONJG(CHTEST(NG))
        NG=NG+1
      ENDDO row
      ENDDO col

      CALLMPI( M_sum_d(GRIDC%COMM,E,1))

    END SUBROUTINE


!*******************************************************************
!
! this subroutine calculates the average electrostatic
! potential at a particular site
! as test charge either the augmentation charges (PAW) or
! a spherical test charge with an automatically chosen
! radius is used
!
!*******************************************************************

    SUBROUTINE CL_SHIFT_PW( GRIDC, LATT_CUR, IRDMAX, &
           T_INFO, P, NCDIJ, CVTOT, ENAUG, IU6)
      USE prec
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE asa
      USE radial
      USE constant
      USE pseudo
      USE cl

      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC    ! grid descriptor
      TYPE (type_info)   T_INFO   ! type descriptor
      TYPE (latt)        LATT_CUR ! lattice descriptor
      TYPE (potcar),TARGET :: P(T_INFO%NTYPD)
      INTEGER NCDIJ               ! number of dimensions of local potential
      RGRID :: CVTOT(DIMREAL(GRIDC%MPLWV),NCDIJ)  ! local potential
      INTEGER  IRDMAX             ! allocation required for augmentation
      REAL(q)  ENAUG              ! cutoff for augmentation grid
      INTEGER  IU6                ! IO unit
!  work arrays
      REAL (q) :: RADIUS(T_INFO%NTYP)
      TYPE (rgrid), TARGET  :: R
      TYPE (rgrid), POINTER :: RP
      REAL(q)   ,ALLOCATABLE ::   DIST(:),DEP(:),POT(:),YLM(:,:),XS(:),YS(:),ZS(:)
      REAL(q) :: VCORE(T_INFO%NIONS)
      INTEGER,ALLOCATABLE ::   NLI(:)
      REAL(q) :: NORM,SUM1,SUM2, RINPL, RCORE
      INTEGER :: LMYDIM, NTYP, LYMAX, NI, NIS, NPS, I, IRDMAX_LOCAL, N, INDMAX
      INTEGER :: INDEX, N1, L1, IND, MAXNL
      CHARACTER (LEN=1) :: LC(5) = (/ "s","p","d","f","g"  /)

      INTEGER, PARAMETER :: NQ=2
      REAL(q)  QQ(NQ)    ! parameters of Besselfunctions
      REAL(q)  A(NQ)     ! parameters of Besselfunctions
      
     
      LMYDIM=1        ! number of lm pairs
!=======================================================================
! first set default grid (used for non PAW potentials)
!=======================================================================
      RCORE = 6/SQRT(ENAUG /RYTOEV)
#ifdef RCORE_MAIN
         RCORE= RCORE_MAIN
#endif

      IRDMAX_LOCAL =4*PI*RCORE**3/3/(LATT_CUR%OMEGA/ &
     &     (GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ))+200

      IRDMAX_LOCAL=MAX(IRDMAX,IRDMAX_LOCAL)

      R%NMAX  =400
      R%RSTART=1E-2
      R%REND  =RCORE
      R%D     =(R%REND/R%RSTART)**(1._q/(R%NMAX-1)) 
      R%H     =LOG(R%D)

      ALLOCATE(R%R(R%NMAX))

      DO I=1,R%NMAX
         R%R(I)=R%RSTART*EXP((I-1)*R%H)
      ENDDO

      R%RMAX=  R%R(R%NMAX)
      R%REND  =RCORE

      CALL SET_SIMP(R)

      ALLOCATE( DIST(IRDMAX_LOCAL),DEP(IRDMAX_LOCAL),POT(IRDMAX_LOCAL), &
     &        YLM(IRDMAX_LOCAL,LMYDIM),NLI(IRDMAX_LOCAL), &
     &        XS(IRDMAX_LOCAL),YS(IRDMAX_LOCAL),ZS(IRDMAX_LOCAL))
!=======================================================================
! loop over all types
!=======================================================================
      VCORE = 0
      NORM  = 0
      NIS =1

      type: DO NTYP=1,T_INFO%NTYP
      !
      ! PAW use the conventional grids (as supplied)
      !
         IF (.NOT. ASSOCIATED(P(NTYP)%QPAW)) THEN
            RP=>R
         ELSE
            RP=> P(NTYP)%R
#ifdef RCORE_MAIN
            RP=>R
#endif            
         ENDIF

         CALL AUG_SETQ(0,RP,QQ,A,LCOMPAT=.FALSE.)

         A=A/PI/4
!=======================================================================
! loop over all ions
!=======================================================================
         RINPL=1._q/GRIDC%NPLWV

         ion: DO NI=NIS,T_INFO%NITYP(NTYP)+NIS-1
!-----------------------------------------------------------------------
! calculate the spherical harmonics YLM and the distance DIST  between
! grid-points and central atom (DEP and POT are work-arrays)
!-----------------------------------------------------------------------
         LYMAX=0

         NPS=1000000  ! (cutoff sphere is given by R%RMAX *(NPS-1)/NPS
                      ! use all points
         RADIUS(NTYP)=RP%RMAX
         CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),RP%RMAX,NPS, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX_LOCAL,INDMAX, &
     &        0.0_q,0.0_q,0.0_q,DIST(1),NLI(1), &
     &        XS(1),YS(1),ZS(1))

         !
         ! in the spin polarised case the up and down potential 
         ! is averaged
         ! 
         DO N=1,INDMAX
            POT(N)=CVTOT(NLI(N),1)
         ENDDO
         IF (NCDIJ==2) THEN
           DO N=1,INDMAX
             POT(N)=(POT(N)+CVTOT(NLI(N),2))/2
           ENDDO
         ENDIF
         IF (NCDIJ==4) THEN
           DO N=1,INDMAX
             POT(N)=(POT(N)+CVTOT(NLI(N),4))/2
           ENDDO
         ENDIF
!=======================================================================
! calculate the integral  int V Y(0,0) Q(0)
!=======================================================================
         CALL SETDEP_B(R,QQ,A,LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))

         SUM1=0
         SUM2=0
!DIR$ IVDEP
!OCL NOVREC
         DO IND=1,INDMAX
            SUM1=SUM1+POT(IND)*DEP(IND)
            SUM2=SUM2+DEP(IND)
         ENDDO
         VCORE(NI)=SUM1*RINPL
         NORM=NORM+SUM2*RINPL
!=======================================================================
      ENDDO ion
      NIS = NIS+T_INFO%NITYP(NTYP)
      ENDDO type
!-----------------------------------------------------------------------
      DEALLOCATE(DIST,DEP,POT,YLM,XS,YS,ZS,NLI)

! reduce VCORE and NORM
      CALLMPI( M_sum_d(GRIDC%COMM , VCORE ,T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM , NORM ,1))

      DO I=1,T_INFO%NIONS
         CL_SHIFT(:,I)=VCORE(I)+CL_SHIFT(:,I)
      ENDDO

      IF (IU6>0 .AND. .NOT.  ACCURATE_CORE_LEVEL_SHIFTS()) THEN
         WRITE(IU6,1) RADIUS
         WRITE(IU6,2) NORM/T_INFO%NIONS,(I,CL_SHIFT(1,I),I=1,T_INFO%NIONS)
         WRITE(IU6,*)
         WRITE(IU6,*)
         WRITE(IU6,*)
      ELSE  IF (IU6>0 ) THEN
         WRITE(IU6,11)
         NIS =1
         DO NTYP=1,T_INFO%NTYP
            DO NI=NIS,T_INFO%NITYP(NTYP)+NIS-1
               IF (ASSOCIATED(P(NTYP)%QPAW)) THEN
                  CALL CALCULATE_MAX_N_L(P(NTYP)%ZCORE, MAXNL)
                  CALL INIT_N_L(INDEX,N1,L1)
                  IND=0
                  WRITE(IU6,'(I4,"-")',ADVANCE='NO') NI
                  DO 
                    IF (.NOT. NEXT_N_L( MAXNL, INDEX, N1, L1)) EXIT
                    WRITE(IU6,'(1I3,A1,1F12.4)',ADVANCE='NO') N1,LC(L1+1),CL_SHIFT(INDEX,NI) 
                    IND=IND+1
                    IF (IND>=5) THEN
                      WRITE(IU6,'(/5X)',ADVANCE='NO')
                      IND=0
                    ENDIF
                  ENDDO
                  WRITE(IU6,*)
            !      WRITE(IU6,12) NI,(CL_SHIFT(J,NI),J=1,MAXNL)
                ENDIF
            ENDDO
            NIS = NIS+T_INFO%NITYP(NTYP)
         ENDDO
         WRITE(IU6,*)
      ENDIF

1     FORMAT(//' average (electrostatic) potential at core'/ & 
          '  the test charge radii are   ',10F8.4,')')
2     FORMAT(  '  (the norm of the test charge is            ',F8.4,')'/ &
          5(I8,F9.4))
11    FORMAT(//' the core state eigenenergies are')
12    FORMAT(I8,5F12.4,/8X,5F12.4,/8X,5F12.4,/8X,5F12.4,/8X)

       
      RETURN
    END SUBROUTINE


!************************ SUBROUTINE SETDEP ****************************
!
! this subroutine interpolates calculates a normalized
! test charge (linear combination of spherical Besselfunctions)
!
!***********************************************************************

    SUBROUTINE SETDEP_B(R,QQ,A,OMEGA,INDMAX,DIST,DEP)
      USE prec
      USE constant
      USE radial

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (rgrid) R
      INTEGER, PARAMETER :: NQ=2
      REAL(q)  QQ(NQ)
      REAL(q)  A(NQ)
       
      REAL(q) DIST(INDMAX)
      REAL(q) DEP(INDMAX)


      FACT= OMEGA
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
         SUM=0
         DO I=1,NQ
            QR=QQ(I)*DIST(IND)
            CALL SBESSEL( QR, BJ, 0)
            SUM=SUM+BJ*A(I)
         ENDDO
         DEP(IND) =SUM*FACT
      ENDDO

    END SUBROUTINE



