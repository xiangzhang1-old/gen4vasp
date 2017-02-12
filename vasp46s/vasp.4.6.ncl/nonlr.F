!#define nonlr_single
#include "symbol.inc"
      MODULE nonlr
      USE prec
      INCLUDE "nonlr.inc"

      LOGICAL LREAL_COMPAT
      INTEGER LREAL_COMPAT_ADD_ONE_GRID_POINT

      CONTAINS

!****************** subroutine NONLR_ALLOC *****************************
! RCS:  $Id: nonlr.F,v 1.7 2003/06/27 13:22:21 kresse Exp kresse $
!
! allocate required arrays
!***********************************************************************

      SUBROUTINE  NONLR_SETUP(NONLR_S,T_INFO,P, LREAL, LSPIRAL)
      USE prec
      USE pseudo
      USE poscar
      IMPLICIT NONE


      TYPE (nonlr_struct) NONLR_S
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      LOGICAL LREAL
      INTEGER, EXTERNAL :: MAXL1
!-MM- spin spiral stuff
      LOGICAL LSPIRAL
!-MM- end of addition

! local var
      INTEGER NIONS,LMDIM,NT

      NIONS=T_INFO%NIONS
      LMDIM=P(1)%LMDIM

      NONLR_S%LREAL  =LREAL
      NONLR_S%NK     =0
      NONLR_S%NTYP   =T_INFO%NTYP
      NONLR_S%NIONS  =T_INFO%NIONS
      NONLR_S%IRALLOC=0
      NONLR_S%NITYP  =>T_INFO%NITYP
      NONLR_S%POSION =>T_INFO%POSION
!-MM- spin spiral stuff
      NONLR_S%LSPIRAL=LSPIRAL
!-MM- end of addition

      ALLOCATE(NONLR_S%LMAX  (NONLR_S%NTYP), &
               NONLR_S%LMMAX (NONLR_S%NTYP), &
               NONLR_S%CHANNELS(NONLR_S%NTYP), &
               NONLR_S%PSRMAX(NONLR_S%NTYP), &
               NONLR_S%BETA  (NONLR_S%NTYP))

      DO NT=1,T_INFO%NTYP
        NONLR_S%LMAX(NT)    = MAXL1(P(NT))
        NONLR_S%LMMAX(NT)   =P(NT)%LMMAX
        NONLR_S%CHANNELS(NT)=P(NT)%LMAX
        NONLR_S%PSRMAX(NT)  =P(NT)%PSRMAX
        NONLR_S%BETA(NT)%PSPRNL=>P(NT)%PSPRNL
        NONLR_S%BETA(NT)%LPS   =>P(NT)%LPS
      ENDDO

      RETURN
      END SUBROUTINE

!****************** subroutine LREAL_COMPAT_MODE ***********************
!
! unfortunately vasp.4.4.X has a bug, insofar that the correct
! radial cutoff is not used in the real space projection
! scheme
! the real space projectors go through zero at a radius r_c corresponding
! to the grid point NPSRNL.
! vasp.4.4.X and vasp.4.5 interpolate the real space projectors one point
! beyond this point!
! The precise behavior can be determined by the flag LREAL_COMPAT
!   LREAL_COMPAT = .TRUE.  vasp.4.4 behavior
!   LREAL_COMPAT = .FALSE. correct new vasp.4.6 behavior
! in vasp.4.6 the default for LREAL_COMPAT is LCOMPAT.
!
!***********************************************************************

      SUBROUTINE LREAL_COMPAT_MODE(IU5, IU0, LCOMPAT)
        USE prec
        IMPLICIT NONE 
        
        INTEGER               :: IU5,IU0    ! input unit
        LOGICAL :: LCOMPAT

        LOGICAL :: LOPEN,LDUM
        INTEGER :: IDUM, N, IERR
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER*1 :: CHARAC
        
        LOPEN=.FALSE.
        OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')

        LREAL_COMPAT=LCOMPAT
        CALL RDATAB(LOPEN,'INCAR',IU5,'LREAL_COMPAT','=','#',';','L', &
             &  IDUM,RDUM,CDUM,LREAL_COMPAT,CHARAC,N,1,IERR)

        IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'Error reading item ''LREAL_COMPAT'' from file INCAR.'
           ENDIF
           LREAL_COMPAT=.FALSE.
        ENDIF

        IF (LREAL_COMPAT) THEN
           LREAL_COMPAT_ADD_ONE_GRID_POINT=1
        ELSE
           LREAL_COMPAT_ADD_ONE_GRID_POINT=0
        ENDIF
        CALL XML_INCAR('LREAL_COMPAT','L',IDUM,RDUM,CDUM,LREAL_COMPAT,CHARAC,N)

      END SUBROUTINE LREAL_COMPAT_MODE

      SUBROUTINE XML_WRITE_LREAL_COMPAT_MODE
        USE prec
        IMPLICIT NONE
        INTEGER :: IDUM
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER*1 :: CHARAC

        CALL XML_INCAR('LREAL_COMPAT','L',IDUM,RDUM,CDUM,LREAL_COMPAT,CHARAC,1)

      END SUBROUTINE XML_WRITE_LREAL_COMPAT_MODE


!****************** subroutine NONLR_ALLOC *****************************
! allocate required arrays
!***********************************************************************

      SUBROUTINE  NONLR_ALLOC(NONLR_S)
      USE pseudo
      IMPLICIT NONE

      TYPE (nonlr_struct) NONLR_S
      INTEGER NIONS,IRMAX

      NIONS = NONLR_S%NIONS
      IRMAX = NONLR_S%IRMAX

      IF (NONLR_S%LREAL) &
      ALLOCATE(NONLR_S%NLIMAX(NIONS), &
               NONLR_S%NLI   (IRMAX,NIONS), &
               NONLR_S%RPROJ (NONLR_S%IRALLOC))
#ifndef gammareal
!-MM- Original allocation statement
!     IF (NONLR_S%LREAL) &
!      ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS))
      IF (.NOT.NONLR_S%LSPIRAL .AND. NONLR_S%LREAL) &
      ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,1))
! Change in allocation of phaser array to accomodate spin spirals
      IF (NONLR_S%LSPIRAL .AND. NONLR_S%LREAL) &
      ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,2))
!-MM- end of alteration
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine NONLR_DEALLOC ***************************
! deallocate arrays
!***********************************************************************

      SUBROUTINE  NONLR_DEALLOC(NONLR_S)
      USE pseudo
      IMPLICIT NONE

      TYPE (nonlr_struct) NONLR_S

      IF (NONLR_S%LREAL) &
      DEALLOCATE(NONLR_S%NLIMAX, &
                 NONLR_S%NLI   , &
                 NONLR_S%RPROJ)
#ifndef gammareal
      IF (NONLR_S%LREAL) &
      DEALLOCATE(NONLR_S%CRREXP)
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine REAL_OPTLAY *****************************
!
! set IRMAX and IRALLOC in the non local PP structure
!
! for the parallel version the subroutine also optimizes
! the layout (i.e. data distribution) of the real space projectior 
! the routine tries to give the same number of grid points
! on all nodes
!
!***********************************************************************

      SUBROUTINE REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,LNOREDIS, &
                             LREALLOCATE,IU6,IU0)

      USE prec
      USE lattice
      USE mgrid
      USE constant
      USE pseudo
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)      GRID
      TYPE (latt)         LATT_CUR
      TYPE (nonlr_struct) NONLR_S
! local work arrays
      INTEGER, ALLOCATABLE :: USED_ROWS(:,:) ! counts how many elements
                          ! must be allocated for one row
      INTEGER, ALLOCATABLE :: REDISTRIBUTION_INDEX(:)
      LOGICAL  LNOREDIS   ! no redistribution allowed
      LOGICAL  LREALLOCATE

      LREALLOCATE=.FALSE.

      IF (.NOT. NONLR_S%LREAL) RETURN
!=======================================================================
! loop over all ions
!=======================================================================
      NLIIND=0
      IRMAX =0
      NIS=1
#ifdef MPI
      ALLOCATE(USED_ROWS(GRID%NGX,GRID%NGY), &
               REDISTRIBUTION_INDEX(GRID%NGX*GRID%NGY) )
      USED_ROWS=0
      IRALLOC=0  ! number of real space proj on local node
#endif

      type: DO NT=1,NONLR_S%NTYP
      IF (NONLR_S%LMMAX(NT)==0) GOTO 600
      LMMAXC=NONLR_S%LMMAX(NT)
      ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!-----------------------------------------------------------------------
! check some quantities
!-----------------------------------------------------------------------
      ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
      D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
      D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
      N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
      N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

      N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
      N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
      N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX
!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
#ifdef MPI
      IND=1

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)+1

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))
      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)+1

      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN
        IND=IND+1
        USED_ROWS(N1P,N2P)=USED_ROWS(N1P,N2P)+LMMAXC
        ! if on local processor add to IRALLOC
        IF (GRID%RL%INDEX(N1P,N2P) /=0) IRALLOC=IRALLOC+LMMAXC
      ENDIF
      ENDDO
      ENDDO; ENDDO
#else
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
      IND=1
      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))
      N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      NCOL=GRID%RL%INDEX(N2P,N3P)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN
        N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
        NCHECK=N1P+(NCOL-1)*GRID%NGX+1
        IF (NCHECK /= 1+N1P+GRID%NGX*(N2P+GRID%NGY* N3P)) THEN
          WRITE(*,*)'REAL_OPT: internal ERROR:',N1P,N2P,N3P, NCOL
          STOP
        ENDIF
        IND=IND+1
      ENDIF
      ENDDO; ENDDO; ENDDO
#endif

      INDMAX=IND-1
      IRMAX =MAX(IRMAX,INDMAX)
      NLIIND=NONLR_S%LMMAX(NT)*INDMAX+NLIIND
!=======================================================================
! end of loop over ions and types
      ENDDO ions
  600 NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO type
!=======================================================================


#ifndef MPI
! to avoid too often reallocation increase values by 10 %
      IF (IRMAX >  NONLR_S%IRMAX) THEN
         NONLR_S%IRMAX   =IRMAX*1.1
         LREALLOCATE=.TRUE.
      ENDIF

      IF( NLIIND > NONLR_S%IRALLOC) THEN
         NONLR_S%IRALLOC =NLIIND*1.1
         LREALLOCATE=.TRUE.
      ENDIF
#else
!=======================================================================
! now redistribute rows in MPI version
!=======================================================================
! first check whether everything is ok
      NLISUM=SUM(USED_ROWS)
      IF (NLISUM /= NLIIND) THEN
        WRITE(*,*)'REAL_OPTLAY: internal error (1)',NLISUM,NLIIND
        STOP
      ENDIF
! setup redistribution index
      NCOL_TOT=GRID%NGY*GRID%NGX
      DO I=1,NCOL_TOT
        REDISTRIBUTION_INDEX(I)=I
      ENDDO
      IF (.NOT. LNOREDIS) THEN
        CALL SORT_REDIS(NCOL_TOT,USED_ROWS(1,1),REDISTRIBUTION_INDEX(1))
        CALL REAL_OPTLAY_GRID(GRID,     REDISTRIBUTION_INDEX,USED_ROWS(1,1),IRALLOC)
      ENDIF
      IF (IRMAX >  NONLR_S%IRMAX) THEN
    ! IRMAX  is the maximum global number, could be improved !!!!
         NONLR_S%IRMAX   =IRMAX  *1.1
         LREALLOCATE=.TRUE.
      ENDIF
      IF( IRALLOC > NONLR_S%IRALLOC) THEN
    ! more savety on parallel machines increase by 40 %
         NONLR_S%IRALLOC =IRALLOC*1.4
         LREALLOCATE=.TRUE.
      ENDIF

      IALLOC_MAX=  IRALLOC
      IALLOC_MIN= -IRALLOC
      CALLMPI (M_max_i( GRID%COMM, IALLOC_MAX, 1))
      CALLMPI (M_max_i( GRID%COMM, IALLOC_MIN, 1))
      CALLMPI (M_sum_i( GRID%COMM, IRALLOC, 1))
      IALLOC_MIN=-IALLOC_MIN
      AKBYTES=1024/8  ! conversion from 8 byte words to kbytes

      IF (.NOT. LNOREDIS .AND. IU6>=0) &
         WRITE(IU6,*)'redistribution in real space done'
      IF (.NOT. LNOREDIS .AND. IU0>=0) &
         WRITE(IU6,*)'redistribution in real space done'
      IF (IU6>=0) &
        WRITE(IU6,1) IRALLOC/AKBYTES,IALLOC_MAX/AKBYTES,IALLOC_MIN/AKBYTES
 1    FORMAT(/' real space projection operators:'/ &
              '  total allocation   :',F10.2,' KBytes'/ &
             '  max/ min on nodes  :',2F10.2/)

      IF (NLISUM /= NLIIND) THEN
        WRITE(*,*)'REAL_OPTLAY: internal error (2)',IRALLOC,NLIIND
        STOP
      ENDIF
      DEALLOCATE(USED_ROWS,REDISTRIBUTION_INDEX)

#endif
      END SUBROUTINE


#ifdef MPI
!
! step through all columns and distribute them onto proc.
! in the manner 1 ... NCPU - NCPU ... 1 - 1 ... NCPU - etc.
!
      SUBROUTINE REAL_OPTLAY_GRID(GRID,REDISTRIBUTION_INDEX,USED_ROWS,IRALLOC)
      USE prec
      USE mgrid

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID
      INTEGER REDISTRIBUTION_INDEX(GRID%NGX*GRID%NGY), &
              USED_ROWS(GRID%NGX*GRID%NGY)
      LOGICAL LUP

      GRID%RL%INDEX= 0
      GRID%RL%NCOL = 0

      NODE_TARGET=0  ! NODE onto which column has to go
      LUP=.TRUE.     ! determines whether NODE_TARGET is increased or decreased
      IRALLOC=0

      NCOL_TOT=GRID%NGY*GRID%NGX
      DO NCOL=1,NCOL_TOT
        IND_REDIS=REDISTRIBUTION_INDEX(NCOL)
        N2=MOD(IND_REDIS-1,GRID%NGX)+1         ! x index (is fast)
        N3=   (IND_REDIS-1)/GRID%NGX+1         ! y index
        IF (LUP) THEN
          IF (NODE_TARGET == GRID%COMM%NCPU) THEN
            LUP=.FALSE.
          ELSE
            NODE_TARGET=NODE_TARGET+1
          ENDIF
        ELSE
          IF (NODE_TARGET == 1) THEN
            LUP=.TRUE.
          ELSE
            NODE_TARGET=NODE_TARGET-1
          ENDIF
        ENDIF

        IND_ON_CPU=(NCOL-1)/GRID%COMM%NCPU+1

        ! element on local node
        ! set up required elements
        IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
          GRID%RL%NCOL=GRID%RL%NCOL+1
          IF (IND_ON_CPU /= GRID%RL%NCOL) THEN
              WRITE(*,*)'REAL_OPTLAY: internal error(3) ',GRID%COMM%NODE_ME,IND_ON_CPU,GRID%RL%NCOL
              STOP
           ENDIF
          GRID%RL%INDEX(N2-1,N3-1)=IND_ON_CPU
          GRID%RL%I2(IND_ON_CPU)=N2 ! I2 contains x index
          GRID%RL%I3(IND_ON_CPU)=N3 ! I3      the y index
          IRALLOC=IRALLOC+USED_ROWS(NCOL)
        ENDIF
      ENDDO

      END SUBROUTINE
#endif


!****************** subroutine RSPHER  *********************************
!
!  subroutine RSPHER calculates the sperical harmonics multiplied
!  by the radial projection operators in real space
!  the result is the full real space projection operator  PROJ
!
!  all ions can be displaced by  a constant shift
!  this makes the calculation of finite differences possible
!
!  full nonlocal pseudopotential is given by
!    RPROJ = 1/Omega ^(1/2) Xi(r-R(N)) Y_lm(r-R(N) Exp(i k r-R(N))
!
!  IZERO :
!    -1 calculate minus the projection operator
!     0 set to projection operator
!     1 add to projection operator
!
!***********************************************************************

      SUBROUTINE RSPHER(GRID,NONLR_S, LATT_CUR )
      USE prec
      USE pseudo
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (potcar)      P(NONLR_S%NTYP)
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER NK

      IZERO=0
      CALL RSPHER_ALL(GRID,NONLR_S, LATT_CUR, LATT_CUR, LATT_CUR, &
                 0.0_q, 0.0_q, 0.0_q, 0)
      RETURN
      END SUBROUTINE

      SUBROUTINE RSPHER_ALL(GRID,NONLR_S,LATT_FIN1, LATT_FIN2, LATT_CUR, &
                 DISX,DISY,DISZ, IDISPL)
      USE prec
      USE pseudo
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      USE wave
      USE asa
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (potcar)      P1(NONLR_S%NTYP)
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_FIN1,LATT_FIN2,LATT_FIN
      INTEGER NK
      INTEGER IDISPL      ! 0 no finite differences, 1 finite differences
! work arrays
      REAL(q),ALLOCATABLE :: DIST(:),XS(:),YS(:),ZS(:),VPS(:),YLM(:,:),VYLM(:)

      LYDIM=MAXVAL(NONLR_S%LMAX)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      LMMAX =MAXVAL(NONLR_S%LMMAX) ! number of nlm indices in the non local potential
      IRMAX=NONLR_S%IRMAX

      ALLOCATE(DIST(IRMAX),XS(IRMAX),YS(IRMAX),ZS(IRMAX),VPS(IRMAX),YLM(IRMAX,LMYDIM), &
          VYLM(IRMAX*LMMAX))

!=======================================================================
! loop over all ions
!=======================================================================
      NLIIND=0
      NIS=1

      type: DO NT=1,NONLR_S%NTYP
      IF (NONLR_S%LMMAX(NT)==0) GOTO 600
      ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1

      ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)
!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
      D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
      D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
      N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
      N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

      N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
      N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
      N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX
      
      VYLM= 0

 dis: DO IDIS=-IDISPL,IDISPL,2

      IF (IDIS==-1) THEN
         LATT_FIN=LATT_FIN1
      ELSE IF (IDIS==1) THEN
         LATT_FIN=LATT_FIN2
      ELSE
         LATT_FIN=LATT_CUR
      ENDIF
!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
#ifdef MPI
      IND=1

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))
      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

      NCOL=GRID%RL%INDEX(N1P,N2P)
      IF (NCOL==0) CYCLE ! not on local node go on
      IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
        WRITE(*,*)'RSPHER: internal ERROR:', &
           GRID%RL%I2(NCOL),N1P+1, GRID%RL%I3(NCOL),N2P+1
        STOP
      ENDIF
!OCL SCALAR
      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))

      XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(XC*XC+YC*YC+ZC*ZC)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN
        X= X1*LATT_FIN%A(1,1)+X2*LATT_FIN%A(1,2)+X3*LATT_FIN%A(1,3)
        Y= X1*LATT_FIN%A(2,1)+X2*LATT_FIN%A(2,2)+X3*LATT_FIN%A(2,3)
        Z= X1*LATT_FIN%A(3,1)+X2*LATT_FIN%A(3,2)+X3*LATT_FIN%A(3,3)
        D=SQRT(X*X+Y*Y+Z*Z)

        N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)
        NONLR_S%NLI (IND,NI) =1+N3P+ GRID%NGZ*(NCOL-1)

        ZZ=Z-DISZ*IDIS
        YY=Y-DISY*IDIS
        XX=X-DISX*IDIS
        ! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
        ! is done using the well known formula  | R+d | = | R | + d . R/|R|
        ! this improves the stability of finite differences considerable
        IF (D<1E-4_q) THEN
          DIST(IND)=1E-4_q
        ELSE
          DIST(IND)=MAX(D-IDIS*(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
        ENDIF

        XS(IND)  =XX/DIST(IND)
        YS(IND)  =YY/DIST(IND)
        ZS(IND)  =ZZ/DIST(IND)
        IND=IND+1
      ENDIF
      ENDDO
      ENDDO
      ENDDO
#else
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
      IND=1
      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))
      N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      NCOL=GRID%RL%INDEX(N2P,N3P)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))

      XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(XC*XC+YC*YC+ZC*ZC)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN
        X= X1*LATT_FIN%A(1,1)+X2*LATT_FIN%A(1,2)+X3*LATT_FIN%A(1,3)
        Y= X1*LATT_FIN%A(2,1)+X2*LATT_FIN%A(2,2)+X3*LATT_FIN%A(2,3)
        Z= X1*LATT_FIN%A(3,1)+X2*LATT_FIN%A(3,2)+X3*LATT_FIN%A(3,3)
        D=SQRT(X*X+Y*Y+Z*Z)

        N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
        NONLR_S%NLI (IND,NI) =N1P+(NCOL-1)*GRID%NGX+1
        IF (NONLR_S%NLI (IND,NI) /= 1+N1P+GRID%NGX*(N2P+GRID%NGY* N3P)) THEN
          WRITE(*,*)'RSHPER internal ERROR:',N1P,N2P,N3P, NCOL
          STOP
        ENDIF
        ZZ=Z-DISZ*IDIS
        YY=Y-DISY*IDIS
        XX=X-DISX*IDIS

        ! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
        ! is done using the well known formula  | R+d | = | R | + d . R/|R|
        ! this improves the stability of finite differences considerable
        IF (D<1E-4_q) THEN
          DIST(IND)=1E-4_q
        ELSE
          DIST(IND)=MAX(D-IDIS*(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
        ENDIF

        XS(IND)  =XX/DIST(IND)
        YS(IND)  =YY/DIST(IND)
        ZS(IND)  =ZZ/DIST(IND)
        IND=IND+1
      ENDIF
      ENDDO
      ENDDO
      ENDDO
#endif
!-----------------------------------------------------------------------
!  compare maximum index with INDMAX
!-----------------------------------------------------------------------
      INDMAX=IND-1
      IF (INDMAX>NONLR_S%IRMAX) THEN
        WRITE(*,*)'internal ERROR: RSPHER:  NONLR_S%IRMAX must be increased to', &
     &            INT(INDMAX*1.1_q)
        STOP
      ENDIF
      NONLR_S%NLIMAX(NI)=INDMAX
!=======================================================================
! now calculate the tables containing the spherical harmonics
! multiplied by the pseudopotential
!=======================================================================
      LYDIM=NONLR_S%LMAX(NT)
      CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)

      LMIND=1
      l_loop: DO L=1,NONLR_S%CHANNELS(NT)
!-----------------------------------------------------------------------
! interpolate the non-local pseudopotentials
! and multiply by (LATT_CUR%OMEGA)^(1/2)
! interpolation is done here using spline-fits this inproves the
! numerical stability of the forces the MIN operation takes care
! that the index is between  1 and NPSRNL
!-----------------------------------------------------------------------
        FAKT= SQRT(LATT_FIN%OMEGA)

!DIR$ IVDEP
!OCL NOVREC
        DO IND=1,INDMAX
          I  =MIN(INT(DIST(IND)*ARGSC)+1,NPSRNL-1)

          REM=DIST(IND)-NONLR_S%BETA(NT)%PSPRNL(I,1,L)
          VPS(IND)=(NONLR_S%BETA(NT)%PSPRNL(I,2,L)+REM*(NONLR_S%BETA(NT)%PSPRNL(I,3,L)+ &
     &         REM*(NONLR_S%BETA(NT)%PSPRNL(I,4,L)+REM*NONLR_S%BETA(NT)%PSPRNL(I,5,L))))*FAKT
        ENDDO

        LL=NONLR_S%BETA(NT)%LPS(L)
        MMAX=2*LL

        ! invert sign for first displacement
        IF (IDIS==-1) THEN
           DO IND=1,INDMAX
              VPS(IND)=-VPS(IND)
           ENDDO
        ENDIF

        LMBASE=LL**2+1

        DO LM=0,MMAX
        DO IND=1,INDMAX
           IBAS = (LMIND-1+LM)*INDMAX
           VYLM(IBAS+IND)=VYLM(IBAS+IND)+VPS(IND)*YLM(IND,LM+LMBASE)
        ENDDO
        ENDDO

        LMIND=LMIND+MMAX+1
        ENDDO l_loop

      IF (LMIND-1/=NONLR_S%LMMAX(NT)) THEN
        WRITE(*,*)'internal ERROR: SPHER:  NONLR_S%LMMAX is wrong',LMIND-1,NONLR_S%LMMAX(NT)
        STOP
      ENDIF

      IF ( NONLR_S%LMMAX(NT)*INDMAX+NLIIND >= NONLR_S%IRALLOC) THEN
        WRITE(*,*)'internal ERROR RSPHER:', &
           'running out of buffer ',NLIIND,INDMAX,NONLR_S%LMMAX(NT),NT,NONLR_S%IRALLOC
        STOP
      ENDIF

      ENDDO dis
!-----------------------------------------------------------------------
! finally store the coefficients
!-----------------------------------------------------------------------
      DO LMIND=1,NONLR_S%LMMAX(NT)
         DO IND=1,INDMAX
            IBAS = (LMIND-1)*INDMAX
            NONLR_S%RPROJ(IND+IBAS+NLIIND)=VYLM(IND+IBAS)
         ENDDO
      ENDDO

      NLIIND= NONLR_S%LMMAX(NT)*INDMAX+NLIIND
!=======================================================================
! end of loop over ions
!=======================================================================
      ENDDO ions
  600 NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO type

      DEALLOCATE(DIST,XS,YS,ZS,VPS,YLM,VYLM)

      RETURN
      END SUBROUTINE

!****************** subroutine PHASER  *********************************
! subroutine PHASER
! recalculates the phase factor for the real-space projectors
! the recalculation is only done if the k-point changes
!***********************************************************************

      SUBROUTINE PHASER(GRID,LATT_CUR,NONLR_S, &
                    NK,WDES,DISX,DISY,DISZ)
      USE prec

      USE lattice
      USE mpimy
      USE mgrid
      USE constant
      USE pseudo
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes)     WDES

#ifdef gammareal
      RETURN
#else
!-----------------------------------------------------------------------
! k-point in kartesian coordiantes
!-----------------------------------------------------------------------
      VKX= WDES%VKPT(1,NK)*LATT_CUR%B(1,1)+WDES%VKPT(2,NK)*LATT_CUR%B(1,2)+WDES%VKPT(3,NK)*LATT_CUR%B(1,3)
      VKY= WDES%VKPT(1,NK)*LATT_CUR%B(2,1)+WDES%VKPT(2,NK)*LATT_CUR%B(2,2)+WDES%VKPT(3,NK)*LATT_CUR%B(2,3)
      VKZ= WDES%VKPT(1,NK)*LATT_CUR%B(3,1)+WDES%VKPT(2,NK)*LATT_CUR%B(3,2)+WDES%VKPT(3,NK)*LATT_CUR%B(3,3)

!-MM- spin spiral stuff
!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply zero when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
      QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
      QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
      QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!     QX=0;QY=0;QZ=0
!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
      IF (NONLR_S%LSPIRAL) THEN 
         NSPINORS=2
      ELSE
         NSPINORS=1
      ENDIF
      
      spinor: DO ISPINOR=1,NSPINORS
!-MM- end of addition
      
!=======================================================================
! loop over all ions
!=======================================================================
      NIS=1

!OCL SCALAR
      type: DO NT=1,NONLR_S%NTYP
      IF (NONLR_S%LMMAX(NT)==0) GOTO 600
      ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!-----------------------------------------------------------------------
! check some quantities
!-----------------------------------------------------------------------
      ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
      D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
      D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+GRID%NGZ+.99_q)-GRID%NGZ
      N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+GRID%NGY+.99_q)-GRID%NGY
      N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+GRID%NGX+.99_q)-GRID%NGX

      N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3)
      N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2)
      N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1)

!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
#ifdef MPI
      IND=1

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))
      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

      NCOL=GRID%RL%INDEX(N1P,N2P)
      IF (NCOL==0) CYCLE
      IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
        WRITE(*,*)'internal ERROR PHASER:', &
           GRID%RL%I2(NCOL),N1P+1, GRID%RL%I3(NCOL),N2P+1
        STOP
      ENDIF

      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN

        ZZ=Z-DISZ
        YY=Y-DISY
        XX=X-DISX

!-MM- begin alteration
! original statement
!       NONLR_S%CRREXP(IND,NI)=EXP(CITPI*(X*VKX+Y*VKY+Z*VKZ))
! change in phaser array to accomodate spin spirals
        NONLR_S%CRREXP(IND,NI,ISPINOR)=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))
!-MM- end of alteration

        IND=IND+1
      ENDIF
      ENDDO
      ENDDO; ENDDO
#else
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
      IND=1
      DO N3=N3LOW,N3HI
      X3=(N3*F3-NONLR_S%POSION(3,NI))
      N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

      DO N2=N2LOW,N2HI
      X2=(N2*F2-NONLR_S%POSION(2,NI))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-NONLR_S%POSION(1,NI))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL+LREAL_COMPAT_ADD_ONE_GRID_POINT) THEN

        ZZ=Z-DISZ
        YY=Y-DISY
        XX=X-DISX

!-MM- begin alteration
! original statement
!       NONLR_S%CRREXP(IND,NI)=EXP(CITPI*(X*VKX+Y*VKY+Z*VKZ))
! change in phaser array to accomodate spin spirals
        NONLR_S%CRREXP(IND,NI,ISPINOR)=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))
!-MM- end of alteration

        IND=IND+1
      ENDIF
      ENDDO; ENDDO; ENDDO
#endif
!=======================================================================
! end of loop over ions
!=======================================================================
      ENDDO ions
  600 NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO type
      
!-MM- conjugate phase alteration for spin down: -q/2 -> q/2
      QX=-QX
      QY=-QY
      QZ=-QZ
      ENDDO spinor
!-MM- end of addition

#endif
      RETURN
      END SUBROUTINE


!****************** subroutine RPRO1    ******************************
!
! this subroutine calculates the scalar product of one wavefunction with
! all projector functions in real space
! thesis gK Equ. (10.36)
!
!*********************************************************************


      SUBROUTINE RPRO1(NONLR_S,WDES1,W1)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1

      REAL(q) RP
! work array
      INTEGER :: IP
      REAL(qn),PARAMETER :: ONE=1,ZERO=0
#ifdef avoidalloc
      REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
      REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE   :: CPROJ(:)

      ALLOCATE(WORK(NONLR_S%IRMAX*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif

      CPROJ=0
!=======================================================================
! loop over ions
!=======================================================================
      LMBASE= 0
      
!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NLIIND= 0
      NIS=1

      typ: DO NT=1,NONLR_S%NTYP
      LMMAXC=NONLR_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

      ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
      INDMAX=NONLR_S%NLIMAX(NI)
      IF (INDMAX == 0) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
         IP=NONLR_S%NLI(IND,NI)+ISPINOR*WDES1%MPLWV
#ifdef gammareal
        WORK(IND)      = REAL( W1%CR(IP) ,KIND=qn)
#else
!-MM- changes to accomodate spin spirals
! original statement
!       CTMP=  W1%CR(IP)*NONLR_S%CRREXP(IND,NI)
        CTMP=  W1%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL)
!-MM- end of alteration               
        WORK(IND)      = REAL( CTMP ,KIND=qn)
        WORK(IND+NONLR_S%IRMAX)=AIMAG(CTMP)
#endif
      ENDDO
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
      DO LM=1,LMMAXC
        SUMR=0
        SUMI=0
!DIR$ IVDEP
!OCL NOVREC
        DO IND=1,INDMAX
          RP   =  NONLR_S%RPROJ(IND+ (LM-1)*INDMAX+NLIIND)
          SUMR = SUMR + WORK(IND)       * RP
#ifndef gammareal
          SUMI = SUMI + WORK(IND+NONLR_S%IRMAX) * RP
#endif
        ENDDO
        CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
      ENDDO
#else
      CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND), &
                   INDMAX, WORK(1) , 1 , ZERO ,  TMP(1,1), 1)
#ifndef gammareal
      CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND), &
                   INDMAX, WORK(1+NONLR_S%IRMAX) , 1 , ZERO ,  TMP(1,2), 1)
#endif

      l_loop: DO LM=1,LMMAXC
#ifdef gammareal
        CPROJ(LM+LMBASE)=TMP(LM,1)*WDES1%RINPL
#else
        SUMR=TMP(LM,1)
        SUMI=TMP(LM,2)
        CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
      ENDDO l_loop
#endif

  100 LMBASE= LMMAXC+LMBASE
      NLIIND= LMMAXC*INDMAX+NLIIND
      ENDDO ion

  600 NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONLR_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition     
      ENDDO spinor
! distribute the projected wavefunctions to nodes
      CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine RPROMU   ******************************
!
!  this subroutine  calculates the projection of a set of
!  bands onto the
!  real space projection operators
!
!*********************************************************************


      SUBROUTINE RPROMU(NONLR_S,WDES1,W1,NSIM,LDO)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1(NSIM)
      LOGICAL            LDO(NSIM)

! work array
      PARAMETER(NLM=101)
      REAL(qn),PARAMETER :: ONE=1,ZERO=0
#if defined (use_malloc)
      REAL(qn) :: WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM)
      GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
      POINTER       (work_p, work), (tmp_p, tmp), (cproj_p, CPROJ)
#else
#ifdef avoidalloc
      REAL(qn) :: WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM)
      GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
      REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE   :: CPROJ(:,:)
#endif
#endif
 
#if defined (use_cray_ptr)
      INTEGER            :: irmax, mplwv
      INTEGER            :: ns_nli(*)
      COMPLEX (q)        :: w1_cr(*), ns_crrexp(*)
      POINTER               (ns_nli_p, ns_nli), (w1_cr_p, w1_cr), &
                            (ns_crrexp_p, ns_crrexp)
#endif

#if defined (use_malloc)
      work_p  = malloc(q*ndata*nonlr_s%irmax*nsim)
      tmp_p   = malloc(q*nlm*ndata*nsim)
      cproj_p = malloc(2*q*wdes1%npro_tot*nsim)
#else
#ifndef avoidalloc
      ALLOCATE(WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM), &
               CPROJ(WDES1%NPRO_TOT,NSIM))
#endif
#endif

      CPROJ=0
!=======================================================================
! loop over ions
!=======================================================================
      LMBASE= 0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

#if defined (use_cray_ptr)
      mplwv = wdes1%mplwv
      irmax = nonlr_s%irmax
#endif

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NLIIND= 0
      NIS=1

      typ: DO NT=1,NONLR_S%NTYP
      LMMAXC=NONLR_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

      ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
#if defined (use_cray_ptr)
      ns_nli_p    = LOC (nonlr_s%nli(1,ni))
      ns_crrexp_p = LOC (nonlr_s%crrexp(1,ni,ispiral))
#endif

!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
      INDMAX=NONLR_S%NLIMAX(NI)
      IF (INDMAX == 0) GOTO 100
      IND0=0
      NPFILL=0

      DO NP=1,NSIM

      IF (LDO(NP)) THEN
#if defined (use_cray_ptr)
      w1_cr_p = LOC (w1(np)%cr)
#endif

!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
#if defined (use_cray_ptr)
         ip = ns_nli(ind) + ispinor * mplwv
#if defined (gammareal)
         work(ind+ind0) = REAL (w1_cr(ip), kind = qn)
#else
         ctmp = w1_cr(ip) * ns_crrexp(ind)
         work(ind+ind0) = REAL (ctmp, kind = q)
         work(ind+irmax+ind0) = AIMAG (ctmp)
#endif

#else
         IP=NONLR_S%NLI(IND,NI)+ISPINOR*WDES1%MPLWV
#ifdef gammareal
         WORK(IND+IND0) = REAL( W1(NP)%CR(IP) ,KIND=qn)
#else
!-MM- changes to accomodate spin spirals
! original statement
!       CTMP=  W1(NP)%CR(IP)*NONLR_S%CRREXP(IND,NI)
        CTMP=  W1(NP)%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL)
!-MM- end of alteration
        WORK(IND+IND0) = REAL( CTMP ,KIND=q)
        WORK(IND+(NONLR_S%IRMAX+IND0))=AIMAG(CTMP)
#endif
#endif
      ENDDO
      IND0=IND0+ndata * NONLR_S%IRMAX
      NPFILL=NPFILL+1
      ENDIF

      ENDDO

!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#ifdef RPROMU_DGEMV
      DO NP=0,ndata*NPFILL-1
         CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND), &
            INDMAX, WORK(1+NONLR_S%IRMAX*NP) , 1 , ZERO ,  TMP(1,NP+1), 1)
      ENDDO
#else
      CALL DGEMMn( 'T', 'N' , LMMAXC,  ndata*NPFILL, INDMAX, ONE, &
                   NONLR_S%RPROJ(1+NLIIND), INDMAX, WORK(1), NONLR_S%IRMAX, &
                   ZERO,  TMP(1,1), NLM )
#endif
      IND0=0
      DO NP=1,NSIM
      IF (LDO(NP)) THEN
      l_loop: DO LM=1,LMMAXC
#ifdef gammareal
        CPROJ(LM+LMBASE,NP)=TMP(LM,1+IND0)*WDES1%RINPL
#else
        SUMR=TMP(LM,1+IND0)
        SUMI=TMP(LM,2+IND0)
        CPROJ(LM+LMBASE,NP)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
      ENDDO l_loop
      IND0=IND0+ndata
      ENDIF
      ENDDO

  100 LMBASE= LMMAXC+LMBASE
      NLIIND= LMMAXC*INDMAX+NLIIND
      ENDDO ion


  600 NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONLR_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition 
      ENDDO spinor

! distribute the projected wavefunctions to nodes
      DO NP=1,NSIM
       IF (LDO(NP)) THEN
         CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
       ENDIF
      ENDDO

#if defined (use_malloc)
      CALL free (work_p)
      CALL free (tmp_p)
      CALL free (cproj_p)
#else
#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine RACCT    ******************************
!
!  this subroutine  calculates the non local part of the gradient for
!  all bands.
!  it is only for performance testing
!
!*********************************************************************

      SUBROUTINE RACCT(NONLR_S,WDES,W,GRID,CDIJ,CQIJ,LMDIM, NK)
      USE prec

      USE wave
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (wavespin)    W
      TYPE (grid_3d)     GRID
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! local variables
      TYPE (wavefun1)  W1(WDES%NSIM)
      LOGICAL ::       LDO(WDES%NSIM)
      COMPLEX(q),ALLOCATABLE :: CWORK(:,:),CWORK2(:)
      REAL(q) :: EVALUE(WDES%NSIM)

      LDO=.TRUE.
      NSIM = WDES%NSIM
      DO NT=1,NONLR_S%NTYP
       IF (NONLR_S%LMMAX(NT)/=0) GOTO 300
      ENDDO
      RETURN

 300  CONTINUE
      ALLOCATE(CWORK(GRID%MPLWV,NSIM),CWORK2(GRID%MPLWV))


    ! setup descriptor
      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)

      NPL=WDES%NPLWKP(NK)
      NGVECTOR=WDES%NGVECTOR(NK)

      DO ISP=1,WDES%ISPIN
      DO N=1,WDES%NBANDS,NSIM
        NUP=MIN(N+NSIM-1,WDES%NBANDS)
	CWORK=0
        DO NN=N,NUP
           NNP=NN-N+1
           CALL SETWAV_(W,W1(NNP),NN,NK,ISP)
           EVALUE(NNP)=W%CELEN(N,1,ISP)
        ENDDO
        CALL RACCMU(NONLR_S,WDES1,W1, LMDIM,CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP),EVALUE,CWORK(1,1), &
                  WDES1%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
        DO NN=N,NUP
           NNP=NN-N+1
           DO  ISPINOR=0,WDES%NRSPINORS-1
           CALL FFTEXT(NGVECTOR,WDES%NINDPW(1,NK),CWORK(1+ISPINOR*GRID%MPLWV,NNP),CWORK2(1+ISPINOR*NGVECTOR),GRID,.FALSE.)
           ENDDO
        ENDDO


      ENDDO
      ENDDO
      DEALLOCATE(CWORK,CWORK2)

      RETURN
      END SUBROUTINE

!****************** subroutine RLACC    ******************************
!
!  subroutine for calculating the non local contribution of
!  the Hamiltonian, using real space projection scheme
!  the result of the wavefunction projected on the projection operatores
!  must be given in CPROJ
!  the result is added to  CRACC
!                !!!!!
!*********************************************************************

      SUBROUTINE RACC(NONLR_S,WDES1,W1,LMDIM,CDIJ,CQIJ,EVALUE,  CRACC)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1
      TYPE (wavedes)     WDES

      COMPLEX(q)  CRACC(WDES1%NPLWVL)
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
               CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
! work arrays
      GDEF :: CRESUL(WDES1%NPROD)

      CALL OVERL1(WDES1, LMDIM,CDIJ,CQIJ, EVALUE, W1%CPROJ(1),CRESUL(1))
      CALL RACC0(NONLR_S,WDES1,CRESUL(1),CRACC(1))

      RETURN
      END SUBROUTINE

!****************** subroutine RACC0   ******************************
!
! this subroutine calculates a linear combination of
! projection operatores in real space
! the result is added to CRACC
!
!*********************************************************************

      SUBROUTINE RACC0(NONLR_S,WDES1,CPROJ_LOC,CRACC)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)     WDES1
      COMPLEX(q) CRACC(WDES1%NPLWVL)
      GDEF   CPROJ_LOC(WDES1%NPROD)

! work array
      REAL(q) RP
      INTEGER IP

      REAL(qn),PARAMETER :: ONE=1,ZERO=0
#ifdef avoidalloc
      REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
      REAL(qn),ALLOCATABLE :: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE    :: CPROJ(:)

      ALLOCATE(WORK(NONLR_S%IRMAX*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif

! merge projected wavefunctions from all nodes (if distributed over 
!   plane wave coefficients)
      CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))
!=======================================================================
! loop over ions
!=======================================================================
      LMBASE= 0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NLIIND= 0
      NIS=1

      typ: DO NT=1,NONLR_S%NTYP
      LMMAXC=NONLR_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

      ion: DO  NI=NIS,NONLR_S%NITYP(NT)+NIS-1
      INDMAX=NONLR_S%NLIMAX(NI)
      IF (INDMAX == 0) GOTO 100
!=======================================================================
! set TMP
!=======================================================================
!DIR$ IVDEP
!OCL NOVREC
      DO L=1,LMMAXC
       CTMP= CPROJ(LMBASE+L)*WDES1%RINPL
       TMP(L,1)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
       TMP(L,2)=AIMAG(CTMP)
#endif
      ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
#ifdef do_loops

! initialise accelerations to 0 (real &  imaginary part seperately)

!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        WORK(IND)=0
#ifndef gammareal
        WORK(IND+NONLR_S%IRMAX)=0
#endif
      ENDDO

      DO LM=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        RP  =  NONLR_S%RPROJ(IND+(LM-1)*INDMAX+NLIIND)
        WORK(IND)      =WORK(IND)      +TMP(LM,1)*RP
#ifndef gammareal
        WORK(IND+NONLR_S%IRMAX)=WORK(IND+NONLR_S%IRMAX)+TMP(LM,2)*RP
#endif
      ENDDO
      ENDDO
#else
      CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND), &
                   INDMAX, TMP(1,1) , 1 , ZERO , WORK(1), 1)
#ifndef gammareal
      CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND), &
                   INDMAX, TMP(1,2) , 1 , ZERO , WORK(1+NONLR_S%IRMAX), 1)
#endif

#endif
!=======================================================================
!  add the non local contribution to the accelerations in real space
!=======================================================================
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
         IP=NONLR_S%NLI(IND,NI)+ISPINOR*WDES1%MPLWV
#ifdef gammareal
        CRACC(IP)= CRACC(IP)+WORK(IND)
#else
!-MM- changes to accomodate spin spirals
! original statement
!       CRACC(IP)= CRACC(IP)+ &
!             CMPLX( WORK(IND) , WORK(IND+NONLR_S%IRMAX) ,KIND=q) *CONJG(NONLR_S%CRREXP(IND,NI,))
        CRACC(IP)= CRACC(IP)+ &
              CMPLX( WORK(IND) , WORK(IND+NONLR_S%IRMAX) ,KIND=q) *CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL))
!-MM- end of alteration
#endif
      ENDDO
      
 100  LMBASE= LMMAXC+LMBASE
      NLIIND= LMMAXC*INDMAX+NLIIND
      ENDDO ion
 600  NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONLR_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition      
      ENDDO spinor

#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine RLACCMU  ******************************
!
!  subroutine for calculating the non local contribution of
!  the Hamiltonian, using real space projection scheme
!  for a set of bands simultaneously
!  the result is added to  CRACC
!                !!!!!
!*********************************************************************

      SUBROUTINE RACCMU(NONLR_S,WDES1,W1, &
     &     LMDIM,CDIJ,CQIJ,EVALUE, CRACC,LD, NSIM, LDO)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)    WDES1
      TYPE (wavedes)    WDES
      TYPE (wavefun1)    W1(NSIM)

      COMPLEX(q) CRACC(LD, NSIM)
      OVERLAP    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      REAL(q)    EVALUE(NSIM)
      LOGICAL    LDO(NSIM)
! work arrays
      GDEF :: CRESUL(WDES1%NPROD,NSIM)

      DO NP=1,NSIM
      IF (LDO(NP)) THEN
         CALL OVERL1(WDES1, LMDIM,CDIJ,CQIJ, EVALUE(NP), W1(NP)%CPROJ(1),CRESUL(1,NP))
!         CALL RACC0(NONLR_S,WDES1,CRESUL(1,NP),CRACC(1,NP))
      ENDIF
      ENDDO
      IF (NSIM/=1) THEN
         CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC,LD, NSIM,LDO)
      ELSE
         CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
      ENDIF

      RETURN
      END SUBROUTINE

!****************** subroutine RACCMU   ******************************
!
! this subroutine calculates a set of linear combination of
! projection operatores in real space
! the result is added to CRACC
!
!*********************************************************************

      SUBROUTINE RACC0MU(NONLR_S,WDES1,CPROJ_LOC,CRACC,LD, NSIM, LDO)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes1)     WDES1
      COMPLEX(q) CRACC(LD,NSIM)
      GDEF   CPROJ_LOC(WDES1%NPROD,NSIM)
      LOGICAL LDO(NSIM)

! work array
      INTEGER, PARAMETER  :: NLM=101

      REAL(qn),PARAMETER :: ONE=1,ZERO=0
#if defined (use_malloc)
      GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
      REAL(qn) :: WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM)
      POINTER       (work_p, work), (tmp_p, tmp), (cproj_p, CPROJ)
#else
#ifdef avoidalloc
      GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
      REAL(qn) :: WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM)
#else
      REAL(qn),ALLOCATABLE :: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#endif

#if defined (use_cray_ptr)
      INTEGER            :: irmax, mplwv
      INTEGER            :: ns_nli(*)
      COMPLEX (q)        :: w1_cr(*), ns_crrexp(*)
      POINTER               (ns_nli_p, ns_nli), (w1_cr_p, w1_cr), &
                            (ns_crrexp_p, ns_crrexp)
#endif

#if defined (use_malloc)
      work_p  = malloc(q*ndata*nonlr_s%irmax*nsim)
      tmp_p   = malloc(q*nlm*2*ndata*nsim)
      cproj_p = malloc(2*q*wdes1%npro_tot*nsim)
#else
#ifndef avoidalloc
      ALLOCATE(WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM), &
               CPROJ(WDES1%NPRO_TOT,NSIM))
#endif
#endif

! merge projected wavefunctions from all nodes

      DO NP=1,NSIM
       IF (LDO(NP)) THEN
         CALL MRG_PROJ(WDES1,CPROJ(1,NP),CPROJ_LOC(1,NP))
       ENDIF
      ENDDO
!=======================================================================
! loop over ions
!=======================================================================
      LMBASE= 0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

#if defined (use_cray_ptr)
      mplwv = wdes1%mplwv
      irmax = nonlr_s%irmax
#endif

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NLIIND= 0
      NIS=1

      typ: DO NT=1,NONLR_S%NTYP
      LMMAXC=NONLR_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

      ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
#if defined (use_cray_ptr)
      ns_nli_p    = LOC (nonlr_s%nli(1,ni))
      ns_crrexp_p = LOC (nonlr_s%crrexp(1,ni,ispiral))
#endif
      INDMAX=NONLR_S%NLIMAX(NI)
      IF (INDMAX == 0) GOTO 100
!=======================================================================
! set TMP
!=======================================================================
      IND0=0
      NPFILL=0
      DO NP=1,NSIM
      IF (LDO(NP)) THEN
!DIR$ IVDEP
!OCL NOVREC

      DO L=1,LMMAXC
       CTMP= CPROJ(LMBASE+L,NP)*WDES1%RINPL
       TMP(L,1+IND0)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
       TMP(L,2+IND0)=AIMAG(CTMP)
#endif
      ENDDO
      IND0=IND0+ndata
      NPFILL=NPFILL+1
      ENDIF

      ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
#ifdef RACCMU_DGEMV
      DO NP=0,ndata*NPFILL-1
         CALL DGEMVn( 'N' , INDMAX, LMMAXC,  ONE , NONLR_S%RPROJ(1+NLIIND), &
            INDMAX,  TMP(1,1+NP) , 1 , ZERO ,  WORK(1+NONLR_S%IRMAX*NP), 1)
      ENDDO
#else
      CALL DGEMMn( 'N' , 'N', INDMAX, ndata*NPFILL, LMMAXC, ONE, &
                   NONLR_S%RPROJ(1+NLIIND), INDMAX, TMP(1,1) , NLM , &
                   ZERO , WORK(1), NONLR_S%IRMAX)
#endif
!=======================================================================
!  add the non local contribution to the accelerations
!=======================================================================
      IND0=0
      DO NP=1,NSIM
      IF (LDO(NP)) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
#if defined (use_cray_ptr)
         ip = ns_nli(ind) + ispinor * mplwv
#if defined (gammareal)
         cracc(ip,np) = cracc(ip,np) + work(ind+ind0)
#else
         cracc(ip,np) = cracc(ip,np) + CMPLX (work(ind+ind0), &
                        work(ind+(irmax+ind0)), q) * &
                        CONJG (ns_crrexp(ind))
#endif
#else
         IP=NONLR_S%NLI(IND,NI)+ISPINOR*WDES1%MPLWV
#ifdef gammareal
        CRACC(IP,NP)= CRACC(IP,NP)+WORK(IND+IND0)
#else
!-MM- changes to accomodate spin spirals
! original statement
!       CRACC(IP,NP)= CRACC(IP,NP)+ &
!           CMPLX( WORK(IND+IND0), WORK(IND+(NONLR_S%IRMAX+IND0)) ,KIND=q)* &
!           CONJG(NONLR_S%CRREXP(IND,NI)
        CRACC(IP,NP)= CRACC(IP,NP)+ &
            CMPLX( WORK(IND+IND0), WORK(IND+(NONLR_S%IRMAX+IND0)) ,KIND=q)* &
            CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL))
!-MM- end of alteration
#endif
#endif
      ENDDO
      IND0=IND0+ndata * NONLR_S%IRMAX
      ENDIF
      ENDDO

 100  LMBASE= LMMAXC+LMBASE
      NLIIND= LMMAXC*INDMAX+NLIIND
      ENDDO ion
 600  NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONLR_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition
      ENDDO spinor

#if defined (use_malloc)
      CALL free (work_p)
      CALL free (tmp_p)
      CALL free (cproj_p)
#else
#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
#endif
      RETURN
      END SUBROUTINE

!****************** subroutine RNLPR     *******************************
! subroutine for calculating the non-local energy per ion
!
! E(ION,k) = SUM(BAND,L,M) Z(ION,BAND,L,M,k) CONJG( Z(ION,BAND,L,M,k))
!
!***********************************************************************

      SUBROUTINE RNLPR(GRID,NONLR_S,P,LATT_FIN1,LATT_FIN2,LATT_CUR,W,WDES, &
     &    LMDIM,NIOND,CDIJ,CQIJ, DISX,DISY,DISZ,ENL)
      USE prec

      USE pseudo
      USE wave
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (potcar)      P(NONLR_S%NTYP)
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W,WTMP
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_FIN1,LATT_FIN2

      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      DIMENSION ENL(NONLR_S%NIONS)
! allocate required work space
      GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)

      ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!=======================================================================
!  calculate the projection operator
!=======================================================================
      NK=1

      CALL RSPHER_ALL(GRID,NONLR_S,LATT_FIN1,LATT_FIN2,LATT_CUR, DISX,DISY,DISZ, 1)

      ENL=0
      WTMP=W
      WTMP%CPROJ => CPROW  ! relink the CPROJ array to temporary workspace

      kpoint: DO NK=1,WDES%NKPTS

      CALL PHASER(GRID,LATT_CUR,NONLR_S, NK,WDES, DISX,DISY,DISZ)
      CALL RPRO(NONLR_S,WDES,WTMP,GRID,NK)

      spin: DO ISP=1,WDES%ISPIN
!=======================================================================
!  sum up to give the non-local energy per ion
!=======================================================================

      band: DO N=1,WDES%NBANDS
      EVALUE=W%CELEN(N,NK,ISP)
      WEIGHT=WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LBASE =ISPINOR *WDES%NPRO/2
      LBASE_=ISPINOR_*WDES%NPRO/2

      NIS=1
      typ: DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 510

         ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
!DIR$ IVDEP
!OCL NOVREC
            DO L=1 ,LMMAXC
               DO LP=1,LMMAXC
                    ENL(NI)=ENL(NI)+WEIGHT*W%CPROJ(LBASE_+LP,N,NK,ISP)*GCONJG(CPROW(LBASE+L,N,NK,ISP))* &
                         WDES%RSPIN*(CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
               ENDDO
            ENDDO
            LBASE = LMMAXC+LBASE
            LBASE_= LMMAXC+LBASE_
         ENDDO ion
510      NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO spin

      ENDDO kpoint

      DEALLOCATE(CPROW)
      RETURN
      END SUBROUTINE

!****************** subroutine STRNLR    *******************************
!
!  subroutine for calculating the non-local contributions to stress,
!  use central differences
!  all components to the stress tensor are calculated
!  except if ISIF = 1
!
!  uncomment CTEST-lines if you want to test finit-differences
!
!***********************************************************************

      SUBROUTINE STRNLR(GRID,NONLR_S,P,LATT_CUR,W,WDES, &
     &    LMDIM,NIOND,CDIJ,CQIJ, ISIF,FNLSIF)
      USE prec

      USE pseudo
      USE wave
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (potcar)      P(NONLR_S%NTYP)
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_FIN1,LATT_FIN2

      DIMENSION FNLSIF(3,3)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

!-----non local part
      DIMENSION ENL(NONLR_S%NIONS)

      DIS=1E-5_q
!TEST
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
!=======================================================================
! initialise non-local forces to zero
!=======================================================================
      FNLSIF=0

!=======================================================================
! calculate the contribution to the energy from the nonlocal
! pseudopotential for undistorted lattice
!=======================================================================
      DISX=0
      DISY=0
      DISZ=0
!=======================================================================
! calculate the contribution to the energy from the nonlocal
! pseudopotential for elongation of each basis-vector
!=======================================================================
      DO IDIR=1,3
      DO JDIR=1,3

      LATT_FIN1=LATT_CUR
      LATT_FIN2=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3; DO J=1,3
          LATT_FIN1%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
          LATT_FIN2%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
        ENDDO; ENDDO
      ELSE
!  all directions
        DO I=1,3
          LATT_FIN1%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
          LATT_FIN2%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
        ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN1)
      CALL LATTIC(LATT_FIN2)

      CALL RNLPR(GRID,NONLR_S,P,LATT_FIN1,LATT_FIN2,LATT_CUR,W,WDES, &
     &    LMDIM,NIOND,CDIJ,CQIJ, DISX,DISY,DISZ,ENL)

      DO NI=1,NONLR_S%NIONS
        FNLSIF(IDIR,JDIR)=FNLSIF(IDIR,JDIR)+ENL(NI)
      ENDDO
!
!  only isotrop pressure terminate loop
!
      IF (ISIF==1) THEN
        FNLSIF(2,2)= FNLSIF(1,1)
        FNLSIF(3,3)= FNLSIF(1,1)
        GOTO 400 ! terminate (not very clean but who cares)
      ENDIF

      ENDDO
      ENDDO
!=======================================================================
! calculation finished  scale pressure
!=======================================================================
  400 CONTINUE
      CALLMPI( M_sum_d(WDES%COMM, FNLSIF, 9))

      FNLSIF=FNLSIF/DIS
!TEST
!      WRITE(*,'(E10.3,3E14.7)')DIS,((FNLSIF(I,J),I=1,3),J=1,3)
!      IF (DIS>1E-10) GOTO 1000
!TEST

      RETURN
      END SUBROUTINE


!****************** subroutine FORNLR    *******************************
!
!  subroutine for calculating the non local contribution
!  to the forces acting onto the ions (using simple finite
!  differences)
!  uncomment CTEST-lines if you want to test finit-differences
!
!***********************************************************************

      SUBROUTINE FORNLR(GRID,NONLR_S,P,LATT_CUR,W,WDES, &
     &    LMDIM,NIOND,CDIJ,CQIJ, FORNL)
      USE prec

      USE pseudo
      USE wave
      USE mpimy
      USE mgrid
      USE lattice
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (potcar)      P(NONLR_S%NTYP)
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W,WTMP
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR

      DIMENSION FORNL(3,NONLR_S%NIONS)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!-----some temporary arrays
      DIMENSION ENL(NONLR_S%NIONS)
      DIMENSION DISPL(3)
! allocate required work space
      GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)

      GDEF,POINTER :: CPROT(:,:,:,:)

      ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))

      DIS=1E-5_q
!TEST
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
!=======================================================================
! initialise non-local forces to zero
!=======================================================================
      FORNL=0
!=======================================================================
! calculate the contribution to the force from the nonlocal
! projection functions for displacement X using central (semianlaytical)
! finite differences (about 9 digits precision)
!=======================================================================
      dir: DO IDIR=1,3
      ENL=0

      DISPL=0
      NK=1
      ! operator= beta(r-R-dis)
      DISPL(IDIR)= DIS
      CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR, DISPL(1),DISPL(2),DISPL(3),1)

      WTMP=W
      WTMP%CPROJ => CPROW       ! relink the CPROJ array to temporary workspace

      kpoint: DO NK=1,WDES%NKPTS
        DIS0=0

        CALL PHASER(GRID,LATT_CUR,NONLR_S, NK,WDES, DIS0,DIS0,DIS0)
        CALL RPRO(NONLR_S,WDES,WTMP,GRID,NK)

        spin: DO ISP=1,WDES%ISPIN

        band: DO N=1,WDES%NBANDS
        EVALUE=W%CELEN(N,NK,ISP)
        WEIGHT=WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)

        spinor: DO ISPINOR=0,WDES%NRSPINORS-1
        DO ISPINOR_=0,WDES%NRSPINORS-1

        LBASE =ISPINOR *WDES%NPRO/2
        LBASE_=ISPINOR_*WDES%NPRO/2
          
        NIS=1
        typ: DO NT=1,WDES%NTYP
           LMMAXC=WDES%LMMAX(NT)
           IF (LMMAXC==0) GOTO 510

           ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
!DIR$ IVDEP
!OCL NOVREC
              DO L=1 ,LMMAXC
                 DO LP=1,LMMAXC
                    ENL(NI)=ENL(NI)+WEIGHT*W%CPROJ(LBASE_+LP,N,NK,ISP)*GCONJG(CPROW(LBASE+L,N,NK,ISP))* &
                         WDES%RSPIN*(CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                 ENDDO
              ENDDO

              LBASE = LMMAXC+LBASE
              LBASE_= LMMAXC+LBASE_
           ENDDO ion
510        NIS = NIS+WDES%NITYP(NT)
        ENDDO typ

        ENDDO
        ENDDO spinor
        ENDDO band
        ENDDO spin
        ENDDO kpoint

        DO NI=1,WDES%NIONS
           NIP=NI_GLOBAL(NI, WDES%COMM_INB)
           FORNL(IDIR,NIP)=FORNL(IDIR,NIP)-ENL(NI)/DIS
        ENDDO

      ENDDO dir

      CALLMPI( M_sum_d(WDES%COMM, FORNL(1,1),NONLR_S%NIONS*3))

!TEST
!      WRITE(*,'(4E20.12)') DIS,FORNL(1,1),FORNL(2,1),FORNL(3,1)
!      STOP
!      GOTO 1000
!TEST

      DEALLOCATE(CPROW)

      RETURN
      END SUBROUTINE

      END MODULE

!****************** subroutine RPRO     ******************************
!
!  this subroutine  calculates the projection of all bands onto the
!  real space projection operators doing a set of
!  bands at the same time
!  SGI crashes when compiling steep.F if this subroutine is
!  inside the MODULE
!  
!*********************************************************************

      SUBROUTINE RPRO(NONLR_S,WDES,W,GRID,NK)
      USE prec

      USE wave
      USE mpimy
      USE mgrid
      USE nonlr
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonlr_struct) NONLR_S
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (wavespin)    W
      TYPE (grid_3d)     GRID
! local variables
      TYPE (wavefun1)  W1(WDES%NSIM)
      LOGICAL ::       LDO(WDES%NSIM)

      LDO=.TRUE.
      NSIM = WDES%NSIM
      DO NT=1,NONLR_S%NTYP
       IF (NONLR_S%LMMAX(NT)/=0) GOTO 300
      ENDDO
      RETURN

 300  CONTINUE
      DO N=1,NSIM
        ALLOCATE(W1(N)%CR(GRID%MPLWV*WDES%NRSPINORS))
      ENDDO

    ! setup descriptor
      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)

      NPL=WDES%NPLWKP(NK)
      NGVECTOR=WDES%NGVECTOR(NK)

      DO ISP=1,WDES%ISPIN
      DO N=1,WDES%NBANDS,NSIM
        NUP=MIN(N+NSIM-1,WDES%NBANDS)
        DO NN=N,NUP
           NNP=NN-N+1
           CALL SETWAV_(W,W1(NNP),NN,NK,ISP)
           DO ISPINOR=0,WDES%NRSPINORS-1
              CALL FFTWAV(NGVECTOR,WDES%NINDPW(1,NK),W1(NNP)%CR(1+ISPINOR*WDES1%MPLWV),W1(NNP)%CW(1+ISPINOR*NGVECTOR),GRID)
           ENDDO
        ENDDO
        IF (NSIM/=1) THEN
           CALL RPROMU(NONLR_S,WDES1,W1,NUP-N+1,LDO)
        ELSE
           CALL RPRO1(NONLR_S,WDES1,W1(1))
        ENDIF
      ENDDO
      ENDDO

      DO N=1,NSIM
        DEALLOCATE(W1(N)%CR)
      ENDDO

      RETURN
      END SUBROUTINE

