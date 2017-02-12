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





!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! Module containing wrapper for scaLAPACK
! written by Gilles de Wijs (gD) and Georg Kresse (gK)
! modified to run on any number of  by Dario Alfe
!
!=======================================================================


 MODULE scala
   USE prec
   USE mpimy
   LOGICAL, PUBLIC :: LscaLAPACK = .FALSE.
   LOGICAL, PUBLIC :: LscaLU     = .FALSE.  ! use parallel LU  decomposition


 CONTAINS

      SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN,N)
      TYPE (communic) COMM
      INTEGER N              ! NxN matrix to be distributed
      COMPLEX(q)    Amatin(N,N)    ! input/output matrix
      END SUBROUTINE

      SUBROUTINE pDSSYEX_ZHEEVX(COMM,AMATIN,W,N)
      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      COMPLEX(q)    AMATIN(N,N)     ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
      END SUBROUTINE

      SUBROUTINE INIT_scala(N, MALLOC)

      IMPLICIT NONE
      INTEGER N,MALLOC
      MALLOC=0

      END SUBROUTINE

      SUBROUTINE INIT_scala_t3d
      END SUBROUTINE

 END MODULE


