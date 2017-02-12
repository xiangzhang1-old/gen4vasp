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






!************************ SUBROUTINE CHGLOC ****************************
! RCS:  $Id: chgloc.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  calculate local charge density n(r)
!
!***********************************************************************

      SUBROUTINE CHGLOC(NBANDS,NKDIM,LDIMP,NIONS,ISPIN, &
     &     PAR,FERWE)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION PAR(NBANDS,NKDIM,LDIMP,NIONS,ISPIN)
      DIMENSION FERWE(NBANDS,NKDIM,ISPIN)

      RETURN
      END
