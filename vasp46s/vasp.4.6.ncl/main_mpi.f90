      MODULE main_mpi
      USE prec
      USE base
      USE mpimy
      IMPLICIT NONE

!***********************************************************************
! RCS:  $Id: main_mpi.F,v 1.3 2002/04/16 07:28:45 kresse Exp $
!
! This module initializes the communication univeres used in VASP
!
! we have five communicator (see below)
!
!***********************************************************************

      TYPE (communic),TARGET :: COMM_WORLD ! our world wide communicator
      TYPE (communic),TARGET :: COMM_CHAIN ! communication between images
      TYPE (communic),TARGET :: COMM       ! one image communicator
      TYPE (communic),TARGET :: COMM_INTER ! between-band communicator
      TYPE (communic),TARGET :: COMM_INB   ! in-band communicator
      INTEGER :: IMAGES=0
      CHARACTER(LEN=10) ::  DIR_APP
      INTEGER           ::  DIR_LEN=0

!***********************************************************************
!
! initialize MPI and read entry IMAGES and NPAR from INCAR
! and sub divide communicators
! because redistribution of nodes
! might be done one must do this as early as possible
!
!***********************************************************************

      CONTAINS
!***********************************************************************
!
! make the  directory entry which is used for
! fileio
!
!***********************************************************************

      SUBROUTINE MAKE_DIR_APP(node)
      INTEGER node

! in principle one can chose here any string one wants to use
! only DIR_LEN must be adjusted
      WRITE (DIR_APP  , "(I1,I1,'/')") MOD(node/10,10),MOD(node,10)
      IF (IMAGES==0) THEN
         DIR_LEN=0
      ELSE
         DIR_LEN=3
      ENDIF

      END SUBROUTINE MAKE_DIR_APP

!***********************************************************************
!
! once unit 6 is open write number of nodes and
! all other parameters
!
!***********************************************************************
      SUBROUTINE WRT_DISTR(IU6)
      INTEGER IU6
      WRITE(IU6,*)'serial version'
      END SUBROUTINE WRT_DISTR

      END MODULE
