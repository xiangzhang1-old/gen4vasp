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
! RCS:  $Id: chain.F,v 1.2 2002/08/14 13:59:37 kresse Exp $
!
! Module which implements the elastic band and the nudged
! elastic band method (for references see below)
! module becomes active if IMAGES tag is read from the INCAR files
!
!**********************************************************************

  MODULE chain
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    IMPLICIT NONE

    REAL(q),ALLOCATABLE,SAVE :: posion_all(:,:,:)
    REAL(q) :: spring=10

!**********************************************************************
!
!  routine for forces between the images on the elastic band
!
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions,posion,toten,force,a,b,iu6)
      INTEGER :: nions
      INTEGER :: iu6
      REAL(q) :: posion(3,nions)
      REAL(q) :: force(3,nions),toten
      REAL(q) :: tangent(3,nions)
      REAL(q) :: a(3,3),b(3,3)
! local variables
      REAL(q) :: x(3),x1(3),x2(3)
      REAL(q) :: e_chain,e_image,norm,norm1,norm2,proj,proj2,d1,d2
      REAL(q) :: force_chain(3,nions),d,const,ftangent
      INTEGER node,i,ni

      IF (images==0) RETURN
      IF (spring==-1000) RETURN
    END SUBROUTINE chain_force

!**********************************************************************
!
! there are several points where a global sum between chains
! is required in order to get correct results
! terms like the total energy / or the kinetic energy
! make only sense in a global and not per image sense
!
!**********************************************************************

    SUBROUTINE sum_chain( value )
      REAL(q) :: value
      IF (images==0) RETURN
      IF (spring==-1000) RETURN

      
    END SUBROUTINE sum_chain

!**********************************************************************
!
! also the logical break conditionions must be
!
!**********************************************************************

    SUBROUTINE and_chain( value )
      LOGICAL :: value
      REAL(q) :: sum

      IF (images==0) RETURN
      IF (spring==-1000) RETURN

! if (1._q,0._q) node is .FALSE., .FALSE. is returned on all 
      IF (value) THEN
         sum=0
      ELSE
         sum=1
      ENDIF
      

      IF (sum>=1) THEN
         value=.FALSE.
      ELSE
         value=.TRUE.
      ENDIF
    END SUBROUTINE and_chain

!**********************************************************************
!
! initialize the chain (repeated image mode)
! read the spring constant
! and  the two outer images, these images are kept fixed
! during the entire simulation
!
!**********************************************************************

    SUBROUTINE chain_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

! needed only temporarily
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics)  :: DYN
      INTEGER     IDUM,IERR,N,idir,node
      CHARACTER*1   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!R.S
      integer tiu6, tiu0
        tiu6 = IO%IU6
        tiu0 = IO%IU0

! quick return, if we are not running in image mode
      IF (images==0) RETURN
    END SUBROUTINE chain_init

!**********************************************************************
!
! returns true if hyper nudged elastic band method is used
!
!**********************************************************************

    FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
      IF (images==0 .OR. spring /= 0 ) THEN
         LHYPER_NUDGE=.FALSE.
      ELSE
         LHYPER_NUDGE=.TRUE.
      ENDIF

    END FUNCTION LHYPER_NUDGE
END MODULE chain
