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





!**********************************************************************
!
! this module can be used if time for allocate and
! deallocation turns out to be a problem
! 
!**********************************************************************


MODULE smart_allocate
  USE prec
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE SMART_ALLOCATE_REAL(A,N)
    REAL(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_REAL

  SUBROUTINE SMART_ALLOCATE_WAVE(A,N)
    COMPLEX(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_WAVE

  SUBROUTINE SMART_ALLOCATE_COMPLEX(A,N)
    COMPLEX(q),POINTER :: A(:)
    INTEGER N

    IF (ASSOCIATED(A)) THEN
       IF (SIZE(A) < N) THEN
          DEALLOCATE(A)
       ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(A)) THEN
       ALLOCATE(A(N))
    ENDIF

  END SUBROUTINE SMART_ALLOCATE_COMPLEX

END MODULE smart_allocate
      
