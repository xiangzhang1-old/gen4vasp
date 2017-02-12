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
#if defined(MPI) || defined(MPI_CHAIN)
      SUBROUTINE INIT_MPI(NPAR,IO)

      TYPE (in_struct)   IO
      INTEGER NPAR
! local
      REAL(q)     RDUM  ; CHARACTER*1   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
      INTEGER     N,IERR
      INTEGER     node
!R.S
      integer tiu6, tiu0
        tiu6 = IO%IU6
        tiu0 = IO%IU0

      CALL M_init(COMM_WORLD)
      IMAGES=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'IMAGES','=','#',';','I', &
     &            IMAGES,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
         WRITE(TIU0,*)'Error reading item ''IMAGES'' from file INCAR.'
         STOP
      ENDIF
! NPAR number of bands distributed over processors, default to a large value
! so that VASP gets best performance in most cases
      NPAR=256
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'NPAR','=','#',';','I', &
     &            NPAR,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
         WRITE(TIU0,*)'Error reading item ''NPAR'' from file INCAR.'
         STOP
      ENDIF

!----------------------------------------------------------------------
! M_divide: creates a 2 dimensional cartesian topology
!  for seperate images or work groups
!  each work group (image) will run VASP independently in one
!  sub directory (01-99) of the current directory
! this mode is required to support either independent calculations
! on parallel machines or the nudged elastic band method
!----------------------------------------------------------------------
      IF (IMAGES>0) THEN
         CALL M_divide( COMM_WORLD, IMAGES, COMM_CHAIN, COMM, .TRUE. )
         CALL M_initc( COMM_WORLD)
         CALL M_initc( COMM)
         CALL M_initc( COMM_CHAIN)
      ELSE
         CALL M_initc( COMM_WORLD)
         COMM=COMM_WORLD
         COMM_CHAIN=COMM
      ENDIF
      IF ( COMM_WORLD%NODE_ME == 1 ) &
      WRITE(TIU0,'(" running on ",I4," nodes")') COMM_WORLD%NCPU

      IF ( COMM_WORLD%NODE_ME == 1 .AND. IMAGES>0 ) &
      WRITE(TIU0,'(" each image running on ",I4," nodes")') COMM%NCPU

      IF (NPAR>=1) THEN
!----------------------------------------------------------------------
! M_divide: creates a 2 dimensional cartesian topology within one
!  work group (image)
!  this is required for simulatneous distribution over bands and plane
!  wave coefficients
!  communicators are created for inter-band in intra-band communication
!  NPAR is the number of bands over which is parallelized
! the resulting layout will be the following (NPAR=4, NCPU=8):
! (MPI uses allways Row-major layout)
! wave1           0 (0,0)        1 (0,1)
! wave2           2 (1,0)        3 (1,1)
! wave3           4 (2,0)        5 (2,1)
! wave4           6 (3,0)        7 (3,1)
! wave5           0 (0,0)        1 (0,1)
! etc.
!
! the sub-communicators COMM_INB are one dimensional communicators
! which allow communication within one row i.e. nodes are grouped to
!  0-1     2-3        4-5      6-7
! for shortness these groups will we called in-band-groups
! there communicators are called in-band-communicator
!
! the sub-communicators COMM_INTER are one dimensional communicators
! which allow communication within one column i.e nodes are grouped
!  0-2-4-6     and    1-3-5-7
! these groups will we called inter-band-groups
!
! the most complicated thing is the FFT of soft chargedensities
! the following algorithm is used:
! the soft chargedensity is calculated in real space,
! charge from all bands is merged to processor 0 and 1 using
! the COMM_INTER
! than a FFT involving all processors is done (on processors 2-7
! no components exist in real space)
! the final result in reciprocal space is defined on all processors
! (see SET_RL_GRID for more information)
!----------------------------------------------------------------------
         CALL M_divide( COMM, NPAR, COMM_INTER, COMM_INB, .FALSE.)
         CALL M_initc( COMM)   ! propably not required but who knows
         CALL M_initc( COMM_INTER)
         CALL M_initc( COMM_INB)
         IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(TIU0,'(" distr:  one band on ",I4," nodes, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU
      ELSE
        COMM_INTER = COMM
        COMM_INB   = COMM
      ENDIF

      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU6 = -1
      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU0 = -1

      IF ( IO%IU0/=-1 .AND. IO%IU6 == -1) THEN
         WRITE(*,*)'internal ERROR: io-unit problem',IO%IU0,IO%IU6
         STOP
      ENDIF

      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)

! if all nodes should write (giving complete mess) do not use the
! following line
      IF (COMM_WORLD%NODE_ME /= COMM_WORLD%IONODE .AND. IO%IU0>0) THEN
         OPEN(UNIT=IO%IU0,FILE=DIR_APP(1:DIR_LEN)//'stdout',STATUS='UNKNOWN')
      ENDIF

      END SUBROUTINE
#endif
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
#ifdef MPI
      IF (IU6>=0) THEN
        WRITE(IU6,'(" running on ",I4," nodes")') COMM_WORLD%NCPU
        IF (IMAGES>0 ) &
        WRITE(IU6,'(" each image running on ",I4," nodes")') COMM%NCPU
        WRITE(IU6,'(" distr:  one band on ",I4," nodes, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU

      ENDIF
#else
      WRITE(IU6,*)'serial version'
#endif
      END SUBROUTINE WRT_DISTR

      END MODULE
