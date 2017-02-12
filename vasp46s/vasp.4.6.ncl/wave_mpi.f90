!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part         --------------------------------
!
!   charge density: full grid mode
!
!
!   charge density complex
!
!
!   wavefunctions: full grid mode
!
!
!   wavefunctions complex
!
!
!   common definitions
!





!************************************************************************
! RCS:  $Id: wave_mpi.F,v 1.3 2001/04/03 10:43:18 kresse Exp $
!
!  this module contains the routines required to communicate
!  wavefunctions and projected wavefunctions
!
!  there are also two quite tricky routines SET_WPOINTER
!  which allow to generate pointer to F77-sequenced arrays
!  There is no guarantee that this will work on all computers, but
!  currently it seems to be ok
!
!***********************************************************************

  MODULE wave_mpi
    USE prec
    USE mpimy
!
! I need an interface block here because I pass the
! first element of a pointer to CPTWFP
! most F90 comiler does not allow such a construct
!
      INTERFACE
      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP
!      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      END SUBROUTINE
      END INTERFACE

      INTERFACE
      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP
!      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      END SUBROUTINE
      END INTERFACE


      ! each communcation package needs a unique identifier (ICOMM)
      ! this is handled by the global variable ICOMM_BASE_HANDLE
      ! which is incremented / decremented by ICOMM_INCREMENT
      ! whenever a redistribution handle is allocated

      INTEGER, PARAMETER :: ICOMM_BASE=10000, ICOMM_INCREMENT=1000
      INTEGER, SAVE :: ICOMM_BASE_HANDLE=ICOMM_BASE

      ! is assyncronous communication allowed or not

      LOGICAL :: LASYNC=.FALSE.

      TYPE REDIS_PW_CTR
         INTEGER :: NB                   ! number of bands to be redistributed
         INTEGER :: NBANDS               ! maximum band index
         INTEGER :: NEXT                 ! next vacant slot
         COMPLEX(q), POINTER :: CPTWFP(:,:)  ! storage for redistribution
         INTEGER,POINTER :: BAND(:)      ! bands that are currently redistributed
         INTEGER,POINTER :: SREQUEST(:,:)! handles for outstanding send requests
         INTEGER,POINTER :: RREQUEST(:,:)! handles for outstanding receive requests
         INTEGER :: ICOMM                ! identifier for communication
         TYPE(communic),POINTER :: COMM  ! communication handle
      END TYPE REDIS_PW_CTR

    CONTAINS
!************************ SUBROUTINE REDIS_PW_ALLOC ********************
! 
! This subroutine allocates the storage required for asyncronous
! redistribution of wave function coeffients
! the redistribution can be initiated for up to NSIM bands 
! at the same time
!
!***********************************************************************

    SUBROUTINE REDIS_PW_ALLOC(WDES, NSIM, H)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)   :: WDES
      INTEGER :: NSIM                  ! number of bands done simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H ! handle
    END SUBROUTINE REDIS_PW_ALLOC


    SUBROUTINE REDIS_PW_DEALLOC(H)
      USE prec
      USE wave
      IMPLICIT NONE

      INTEGER :: NSIM                  ! number of bands done simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H ! handle
      INTEGER :: IND
    END SUBROUTINE REDIS_PW_DEALLOC

  END MODULE wave_mpi

!=======================================================================
!
! this routine returns a pointer to an SEQUENCED F77 like array
! with a given storage convention
! 1. dimension is N1 second (1._q,0._q) N2
!
!=======================================================================

      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      CW_P => CPTWFP(:,:)

      END SUBROUTINE


      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CPTWFP)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CPTWFP(N1,N2)
      CW_P => CPTWFP(:,:)

      END SUBROUTINE

!=======================================================================
!  distribute projector part of wavefunctions over all 
!=======================================================================

      SUBROUTINE DIS_PROJ(WDES1,CPROJ,CPROJ_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPROJ(WDES1%NPRO_TOT)
      COMPLEX(q) CPROJ_LOCAL(WDES1%NPRO)
      CPROJ_LOCAL(1:WDES1%NPRO)=CPROJ(1:WDES1%NPRO)
      END SUBROUTINE

!=======================================================================
!  merge projector part wavefunctions from all 
!  definitely not the most efficient implementation
!  but it works
!=======================================================================

      SUBROUTINE MRG_PROJ(WDES1,CPROJ,CPROJ_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPROJ(WDES1%NPRO_TOT)
      COMPLEX(q) CPROJ_LOCAL(WDES1%NPRO)
      CPROJ(1:WDES1%NPRO)=CPROJ_LOCAL(1:WDES1%NPRO)

      END SUBROUTINE

!=======================================================================
!  distribute plane wave part of (1._q,0._q) wavefunction over all COMM_INB 
!  the supplied wavefunction must have the standard serial layout and
!  must be defined on the node COMM_INB%IONODE
!=======================================================================

      SUBROUTINE DIS_PW(WDES1,CPTWFP,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      CW_LOCAL(1:WDES1%NPL)=CPTWFP(1:WDES1%NPL)
      END SUBROUTINE

!=======================================================================
! distribute (1._q,0._q) specific band to 
! only required after i.e. reading wavefunctions
! it is sufficient if wavefunctions are defined on master node
! but mind, that all  must call this communication routine
!=======================================================================

      SUBROUTINE DIS_PW_BAND(WDES1, NB, CPTWFP, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      CALL DIS_PW(WDES1,CPTWFP,CW_LOCAL(1,NB))
      END SUBROUTINE


!=======================================================================
!  merge plane wave part of wavefunctions from all COMM_INB 
!=======================================================================

      SUBROUTINE MRG_PW(WDES1,CPTWFP,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      CPTWFP(1:WDES1%NPL)=CW_LOCAL(1:WDES1%NPL)
      END SUBROUTINE

!=======================================================================
! merge  band NB from all  into a local copy CPTWFP
! only required for i.e. writing wavefunctions
!=======================================================================

      SUBROUTINE MRG_PW_BAND(WDES1, NB, CPTWFP, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CPTWFP(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      CALL MRG_PW(WDES1,CPTWFP,CW_LOCAL(1,NB))
      END SUBROUTINE

!=======================================================================
! merge eigenvalues from all 
!=======================================================================

      SUBROUTINE MRG_CEL(WDES,W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
      END SUBROUTINE

!=======================================================================
! merge partial occupancies from all 
!=======================================================================

      SUBROUTINE MRG_FERWE(WDES,W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
      END SUBROUTINE



!************************ SUBROUTINE REDIS_PW **************************
!
! redistribute plane wave coefficients from over band to over
! plane wave coefficient or vice versa
! this operation is done in place to reduce storage demands
!
! mind that if the routine is called twice the original distribution
! is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PW(WDES1, NBANDS, CPTWFP)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      COMPLEX(q) :: CPTWFP(WDES1%NRPLWV*NBANDS)
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
! local variables
      INTEGER :: NRPLWV,N,NB,INFO

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_ALL **********************
!
! redistribute all coefficients 
!
!***********************************************************************

      SUBROUTINE REDIS_PW_ALL(WDES, W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_START ********************
!
! redistribute plane wave coefficients from over band to over
! plane wave coefficient or vice versa asyncronously 
! (1._q,0._q) band is initiated
!
!***********************************************************************

      SUBROUTINE REDIS_PW_START(WDES, CPTWFP, BANDINDEX, H)
      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND

      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_STOP  ********************
!
! finish redistribution of  plane wave coefficients from over band to over
! plane wave coefficient
!
!***********************************************************************

      SUBROUTINE REDIS_PW_STOP(WDES, CPTWFP, BANDINDEX, H)
      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CPTWFP(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND,N
      END SUBROUTINE


!************************ SUBROUTINE REDIS_PROJ ************************
!
! redistribute projector part of wave function from over band to over
! plane wave coefficient or vice versa
! this operation is done in place to reduce storage demands
!
! mind that if the routine is called twice the original distribution
! is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PROJ(WDES1, NBANDS, CPROJ)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      COMPLEX(q) :: CPROJ(WDES1%NPROD*NBANDS)
      COMPLEX(q), ALLOCATABLE :: CWORK(:)

      INTEGER, PARAMETER :: MCOMP=2
! local variables
      END SUBROUTINE


!************************ SUBROUTINE SET_NPL_NPRO **********************
!
! set the local number of plane waves on a node after redistribution
! of plane wave coefficients and projected wavefunctions
! on entry the global number of plane waves must be given
!
!***********************************************************************

      SUBROUTINE SET_NPL_NPRO(WDES1, NPL, NPRO)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NPL,NPRO,NPL_REM,NPRO_REM,I
! local variables
      NPL_REM =NPL
      NPRO_REM=NPRO
      NPL=WDES1%NRPLWV
      NPRO=WDES1%NPROD

      END SUBROUTINE
