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






!======================================================================
! RCS:  $Id: mpi.F,v 1.6 2003/06/27 13:22:20 kresse Exp kresse $
!
! dummy module if MPI is not used
! a few files will not compile with this dummy module
! i.e. fftmpi.F fftmpi_map.F
!
!======================================================================
      MODULE mpimy
      TYPE communic
        INTEGER nup
      END TYPE
      CONTAINS
      SUBROUTINE mpi_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
      END SUBROUTINE
      END MODULE


