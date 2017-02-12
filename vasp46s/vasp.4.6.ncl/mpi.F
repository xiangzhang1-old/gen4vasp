#include "symbol.inc"
#undef STOP

#if ! defined(MPI) && ! defined(MPI_CHAIN)
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


#else

!=======================================================================
! 
! in most MPI implementations the collective communcation
! is slow, hence by default we avoid them 
! if you want to use them define use_collective (here or in the makefile)
!
! use_collective use MPI_alltoall 
! avoid_async    tries to post syncronised send and read operations
!                such that collisions are avoided, this is usually slower
! PROC_GROUP     does communication is block wise in a group of
!                roughly PROC_GROUP processors
!                this reduces collisions
!
! in addition our own implementation of the collective communication
! routines allow 
!=======================================================================

!#define use_collective
!#define avoid_async

#ifndef MPI_BLOCK
#define MPI_BLOCK   1000
!
! blocking over group of processors 
! on a Gigabit ethernet this can improve performance substantially 
! without hardware flowcontrol
! when hardware flowcontrol is enabled PROC_GROUP is usually not
! required
!
!#define PROC_GROUP 4
#endif

!=======================================================================
!
! MPI communication routines for VASP
! all communication should be done using this interface to allow
! adaption of other communication routines
! routines were entirely rewritten by Kresse Georg,
! but functionallity is similar to a module written by Peter
! Lockey at Daresbury
!
!======================================================================
      MODULE mpimy
      USE prec

      INCLUDE "mpimy.inc"
! Standard MPI include file.
! I would like to have everything in the header but "freaking" SGI
! compiler can not handle this, thus I have to use an include file
      INCLUDE "mpif.h"
! There are no global local sum routines in MPI, thus some workspace
! is required to store the results of the global sum
      INTEGER,PARAMETER ::  NZTMP=MPI_BLOCK/2, NDTMP=MPI_BLOCK, NITMP=MPI_BLOCK
! workspace for integer, complex, and real
      COMPLEX(q),SAVE :: ZTMP_m(NZTMP)
      REAL(q),SAVE    :: DTMP_m(NDTMP)
      INTEGER,SAVE    :: ITMP_m(NITMP)

#ifndef IFC
      EQUIVALENCE (DTMP_m,ZTMP_m)
      EQUIVALENCE (ITMP_m,ZTMP_m)
#endif

      CONTAINS
!----------------------------------------------------------------------
!
! M_init: initialise the basic communications
! (number of nodes, determine ionode)
!
!----------------------------------------------------------------------
!
      SUBROUTINE M_init( COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

     TYPE(communic) COMM
      INTEGER i, ierror

      call MPI_init( ierror )
      IF ( ierror /= MPI_success ) THEN
         WRITE(*,*) 'Initpm: Error in MPI_init'
        STOP
      ENDIF
!
! initial communicator is world wide
! set only NCPU, NODE_ME and IONODE
! no internal setup done at this point
!
      COMM%MPI_COMM= MPI_comm_world

      call MPI_comm_rank( COMM%MPI_COMM, COMM%NODE_ME, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_rank',ierror)
      COMM%NODE_ME= COMM%NODE_ME+1

      call MPI_comm_size( COMM%MPI_COMM, COMM%NCPU , ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_size',ierror)

      COMM%IONODE = 1

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_divide: creates a 2 dimensional cartesian topology
!  and a communicator along rows and columns of the process matrix
!
!----------------------------------------------------------------------

      SUBROUTINE M_divide( COMM, NPAR, COMM_INTER, COMM_INB, reorder)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM, COMM_INTER, COMM_INB, COMM_CART
      INTEGER NPAR,NPAR_2
      INTEGER, PARAMETER :: ndims=2
      INTEGER :: dims(ndims)
      LOGICAL :: periods(ndims), reorder, remain_dims(ndims)
      INTEGER :: ierror

      IF (NPAR >= COMM%NCPU) NPAR=COMM%NCPU
      dims(1)       = NPAR
      dims(2)       = COMM%NCPU/ NPAR
      IF (dims(1)*dims(2) /= COMM%NCPU ) THEN
         WRITE(0,*) 'M_divide: can not subdivide ',COMM%NCPU,'nodes by',NPAR
      ENDIF

      periods(ndims)=.FALSE.

      CALL MPI_Cart_create( COMM%MPI_COMM , ndims, dims, periods, reorder, &
                COMM_CART%MPI_COMM , ierror)
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_create', ierror)
! create the in-band communicator
      remain_dims(1)= .FALSE.
      remain_dims(2)= .TRUE.

      CALL MPI_Cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INB%MPI_COMM, ierror )
      IF ( ierror /= MPI_success )&
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_sub (1) ', ierror)

! create the inter-band communicator
      remain_dims(1)= .TRUE.
      remain_dims(2)= .FALSE.

      CALL MPI_Cart_sub( COMM_CART%MPI_COMM, remain_dims, COMM_INTER%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Dividepm: Error in MPI_Cart_sub (2) ', ierror)
! overwrite initial communicator by new one
      COMM=COMM_CART

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_initc: initialise a communicator
!
!----------------------------------------------------------------------

      SUBROUTINE M_initc( COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror
      INTEGER id_in_group

      call MPI_comm_rank( COMM%MPI_COMM, id_in_group, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_rank',ierror)

      call MPI_comm_size( COMM%MPI_COMM, COMM%NCPU, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_comm_size',ierror)

      COMM%NODE_ME = id_in_group + 1
      CALL init_hard_ids(COMM)

      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)
      COMM%IONODE =1
      END SUBROUTINE

!----------------------------------------------------------------------
!
! init_hard_ids: map the virtual (MPI node_id) to the real node_id
!  (T3E, T3D specific)
!
!----------------------------------------------------------------------

      SUBROUTINE init_hard_ids(COMM)
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, tmp, ierror, id_in_group
      INTEGER, ALLOCATABLE :: hid_tmp(:)
#ifdef T3D_SCA
      INTEGER,INTRINSIC :: MY_PE
!
!  map the virtual (MPI node_id) to the real node_id
!  this is required to support communicators defined on a sub-grid
!  using shmem
!
      ALLOCATE( COMM%hid(0: COMM%NCPU-1), hid_tmp(0: COMM%NCPU-1) )
      hid_tmp =0
      hid_tmp( COMM%NODE_ME-1) = MY_PE()

      call MPI_allreduce( hid_tmp(0), COMM%hid(0), COMM%NCPU, &
              MPI_integer, MPI_sum, COMM%MPI_COMM, ierror )
       IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('inittree: Error in MPI_allreduce',ierror)
      DEALLOCATE( hid_tmp )
#endif
      END SUBROUTINE


      END MODULE
!======================================================================
!
! all other routines are often called with either
! real or complex arrays, vectors or scalars, so I can not put
! them into the F90 module
!
!======================================================================
!----------------------------------------------------------------------
!
! M_exit: exit MPI, very simple just exit MPI
!
!----------------------------------------------------------------------

      SUBROUTINE M_exit()
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER i, ierror, id_in_group

!      call MPI_barrier(MPI_comm_world, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Exitpm: Error in MPI_barrier',ierror)

      call MPI_finalize( ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('Exitpm: Error in MPI_finalize',ierror)
      STOP

      RETURN
      END SUBROUTINE


!----------------------------------------------------------------------
!
! M_stop: exits MPI and program because of error, a message is
! printed
!
!----------------------------------------------------------------------

      SUBROUTINE M_stop(message)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER*(*) message
      INTEGER ierror

      WRITE (*,*) message

      call MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      RETURN
      END SUBROUTINE


      SUBROUTINE M_stop_ierr(message, ierror)
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      CHARACTER*(*) message
      INTEGER ierror

      WRITE (*,*) message, ierror

      call MPI_abort(MPI_comm_world , 1, ierror )
      STOP

      RETURN
      END SUBROUTINE

!======================================================================
!
! Send and Receive routines, map directly onto MPI
!
!======================================================================

!----------------------------------------------------------------------
!
! M_send_i: send n integers store in ivec to node node
!
!----------------------------------------------------------------------

      SUBROUTINE M_send_i (COMM, node, ivec, n)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)

      INTEGER status(MPI_status_size), ierror

      call MPI_send( ivec(1), n, MPI_integer, node-1, 200, &
     &               COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_send returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! M_recv_i: receive n integers into array ivec from node node
!
!----------------------------------------------------------------------

      SUBROUTINE M_recv_i(COMM, node, ivec, n )
      USE prec
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER node, n
      INTEGER ivec(n)
      INTEGER status(MPI_status_size), ierror

      call MPI_recv( ivec(1), n, MPI_integer , node-1, 200, &
     &               COMM%MPI_COMM, status, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_recv returns',ierror)

      RETURN
      END


!======================================================================
!
! global sum and maximum routines
!
!======================================================================

!----------------------------------------------------------------------
!
! M_sum_i: performs a global sum on n integers in vector ivec
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_i(COMM, ivec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! check whether n is sensible
      IF (n==0 .OR. COMM%NCPU==1 ) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sum_i ', n
         STOP
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         call MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
     &                       MPI_sum, COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_integer) returns',ierror)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_max_i: performs a global max on n integers in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_i(COMM, ivec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER ivec(n)
      INTEGER j,k

      INTEGER ierror, status(MPI_status_size), ichunk

! check whether n is sensible
      IF (n==0 .OR. COMM%NCPU==1 ) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_max_i ', n
         STOP
      END IF

!  there is no inplace global max in MPI, thus we have to use
!  a work array

      DO j = 1, n, NITMP
         ichunk = MIN( n-j+1 , NITMP)

         call MPI_allreduce( ivec(j), ITMP_m(1), ichunk, MPI_integer, &
                             MPI_max, COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_max, MPI_integer) returns',ierror)

         DO k = 0, ichunk-1
            ivec(j+k) = ITMP_m(k+1)
         ENDDO

      ENDDO

      RETURN
      END



!----------------------------------------------------------------------
!
! M_max_d: performs a global max search on n doubles in vector vec
!
!----------------------------------------------------------------------

      SUBROUTINE M_max_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk
! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_max_d ', n
         STOP
      END IF

!  there is no inplace global max in MPI, thus we have to use
!  a work array
      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         call MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_max, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_max, MPI_double_prec) returns',ierror)
#if defined (use_fastbcopy)
         CALL FASTBCOPY8(ichunk , DTMP_m(1), vec(j))
#else
         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
#endif

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! M_sumb_d: performs a global sum on n doubles in vector vec
!  uses MPI_allreduce which is usually very inefficient
!  faster alternative routines can be found below 
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      INTEGER j

      INTEGER  ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sumb_d ', n
         STOP
      END IF

!  there is no inplace global sum in MPI, thus we have to use
!  a work array

      DO j = 1, n, NDTMP
         ichunk = MIN( n-j+1 , NDTMP)

         call MPI_allreduce( vec(j), DTMP_m(1), ichunk, &
                             MPI_double_precision, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_double_prec) returns',ierror)

#if defined (use_fastbcopy)
         CALL FASTBCOPY8(ichunk, DTMP_m(1), vec(j))
#else
         CALL DCOPY(ichunk , DTMP_m(1), 1 ,  vec(j) , 1)
#endif

      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
!
! to make live easier, a global sum for scalars
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_s(COMM, n, v1, v2, v3, v4)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n),v1,v2,v3,v4

      vec=0

      IF (n>0) vec(1)=v1
      IF (n>1) vec(2)=v2
      IF (n>2) vec(3)=v3
      IF (n>3) vec(4)=v4
      IF (n>4) THEN
          WRITE(*,*) 'internal ERROR: invalid n in M_sum_s ', n
          STOP
      END IF

      CALL M_sumb_d(COMM, vec, n)

      IF (n>0) v1=vec(1)
      IF (n>1) v2=vec(2)
      IF (n>2) v3=vec(3)
      IF (n>3) v4=vec(4)
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_z: performs a global sum on n complex items in vector vec
!  uses MPI_allreduce which is usually very inefficient
!  faster alternative routines can be found below 
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumb_z(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)
      INTEGER j

      INTEGER ierror, status(MPI_status_size), ichunk

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

! check whether n is sensible
      IF (n==0) THEN
         RETURN
      ELSE IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid n in M_sumb_z ', n
         STOP
      END IF


!  there is no inplace global sum in MPI, thus we have to use
!  a work array
      DO j = 1, n, NZTMP
         ichunk = MIN( n-j+1 , NZTMP)

         call MPI_allreduce( vec(j), ZTMP_m(1), ichunk, &
                             MPI_double_complex, MPI_sum, &
                             COMM%MPI_COMM, ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: MPI_allreduce (MPI_sum, MPI_double_complex) returns',ierror)

#if defined (use_fastbcopy)
         CALL FASTBCOPY8(ichunk , ZTMP_m(1),  vec(j))
#else
         CALL ZCOPY(ichunk , ZTMP_m(1), 1 ,  vec(j) , 1)
#endif

      ENDDO

      RETURN
      END

!======================================================================
!
! Global Copy Routines
!
!======================================================================

!----------------------------------------------------------------------
!
! M_bcast_i: copy n integers from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n, MPI_integer, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_i) returns',ierror)


      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_i_from: copy n integers from node inode to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_i_from(COMM, vec, n , inode)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER inode
      INTEGER vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n, MPI_integer, inode-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_i) returns',ierror)

      RETURN
      END


!----------------------------------------------------------------------
!
! M_bcast_d: copy n double precession from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_d(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_precision, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_d) returns',ierror)


      RETURN
      END


!----------------------------------------------------------------------
!
! M_bcast_z: copy n double precession complex from root to all nodes
!
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z(COMM, vec, n )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_complex, COMM%IONODE-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_z) returns',ierror)


      RETURN
      END

!----------------------------------------------------------------------
!
! M_bcast_z_from: copy n double precession complex from inode to all nodes
! 
! needed for exchange operator
! this routine was inserted by Robin Hirschl on 29.8.2002
!----------------------------------------------------------------------

      SUBROUTINE M_bcast_z_from(COMM, vec, n, inode )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) vec(n)

      INTEGER inode,ierror, status(MPI_status_size)

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      IF (n<0) THEN
         WRITE(*,*) 'internal ERROR: invalid vector n in M_bcast_i'
         STOP
      END IF

      call MPI_bcast( vec(1), n,  MPI_double_complex, inode-1, COMM%MPI_COMM, &
     &                ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: MPI_bcast (M_bcast_z_from) returns',ierror)


      RETURN
      END


!======================================================================
!
! Global Exchange Routine
!
!======================================================================

!----------------------------------------------------------------------
!
! M_alltoallv_z: complex global exchange routine which maps directly onto
!  MPI_alltoallv
!  since MPI_alltoallv is usually very slow, an alternative implementation
!  optimised for clusters exists
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_z(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv, rprcv)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      COMPLEX(q) xsnd(*)         ! send buffer
      COMPLEX(q) xrcv(*)         ! receive buffer

      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
      INTEGER nsnd(COMM%NCPU+1) ! number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) ! number of data recv from each node

!----------------------------------------------------------------------
#ifdef T3D_SMA
!----------------------------------------------------------------------
      INTEGER i, ierror,jnode, j, ndata, shmem_st, iold
      INTEGER,EXTERNAL ::  shmem_put, shmem_udcflush, GET_D_STREAM


      ndata = nsnd(COMM%NODE_ME)
      IF (ndata >= 0) THEN
         ! local copy (if receiver == sender)
         DO j = 1,ndata
            xrcv(prcv(COMM%NODE_ME)+j) = xsnd(psnd(COMM%NODE_ME)+j)
         ENDDO
      ENDIF

      ! syncronize all nodes involved in communication
      iold = GET_D_STREAM()           ! read current streaming state
      CALL QUIET_D_STREAM()           ! disable  streaming
      CALL MPI_barrier( COMM%MPI_COMM, ierror )

      DO i = 1, COMM%NCPU

        ! send data to jnode (deadlock free)
        jnode = 1 + MOD(IEOR((i-1),(COMM%NODE_ME-1)), COMM%NCPU)

        ndata = nsnd(jnode)
        IF (ndata >= 0) THEN
           IF( jnode /= COMM%NODE_ME) THEN
          ! remote put
              shmem_st = shmem_put(xrcv(rprcv(jnode)+1), &
                             xsnd(psnd(jnode)+1), &
                             (ndata*2), COMM%hid(jnode-1))
            ELSE

            ENDIF
        ENDIF

      ENDDO

      CALL MPI_barrier( COMM%MPI_COMM, ierror ) ! wait for other nodes
      CALL SET_D_STREAM(iold)                   ! restore streaming state
      shmem_st = shmem_udcflush()               ! flush the data cache

!----------------------------------------------------------------------
#else
#if defined (use_collective) || defined (use_symmetric_heap)
!----------------------------------------------------------------------
      INTEGER ierror

      call MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_double_complex, &
                          xrcv(1), nrcv(1), prcv(1), MPI_double_complex, &
                          COMM%MPI_COMM, ierror )

      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z MPI_alltoallv returns',ierror)

!----------------------------------------------------------------------
#else
!----------------------------------------------------------------------
      INTEGER ierror,sndcount,rcvcount,prcv_,psnd_,i,in
      INTEGER :: tag=201
      INTEGER :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=MPI_BLOCK/2
      INTEGER       :: block, p, nstat
      INTEGER       :: maxsnd_rcvcount

      maxsnd_rcvcount=MAX(MAXVAL(nsnd(1:COMM%NCPU)),MAXVAL(nrcv(1:COMM%NCPU)))

    DO block = 0, maxsnd_rcvcount-1, max_
      p        = 1 + block   ! pointer to the current block base address
      nstat    = 0
      
      ! initiate the receive and send on all nodes
      ! local copy is done below
      ! handle remaining  NCPU-1 packages
      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network
         
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         rcvcount=MIN(max_, nrcv(i+1)-block)
         prcv_   =prcv(i+1)

         IF (rcvcount>0) THEN
            nstat=nstat+1
            call MPI_irecv( xrcv(prcv_+p), rcvcount, MPI_double_complex, &
                       i, tag,  COMM%MPI_COMM, request(nstat), ierror )
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoallv_d irecv returns',ierror)
         ENDIF
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         ! recommended for T3E 
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         sndcount=MIN(max_, nsnd(i+1)-block)
         psnd_   =psnd(i+1)

         IF (sndcount>0) THEN
            nstat=nstat+1
            call MPI_isend( xsnd(psnd_+p), sndcount, MPI_double_complex, &
     &                  i, tag,  COMM%MPI_COMM, request(nstat), ierror )
            IF ( ierror /= MPI_success ) &
               CALL M_stop_ierr('ERROR: M_alltoallv_d irecv returns',ierror)
         ENDIF
      ENDDO

     ! local memory copy for data kept on the local node

      sndcount=MIN(max_, nsnd(COMM%NODE_ME)-block)
      prcv_   =prcv(COMM%NODE_ME)
      psnd_   =psnd(COMM%NODE_ME)
      
      IF (sndcount>0) &
#if defined (use_fastbocpy)
         CALL FASTBCOPY8( sndcount, xsnd( psnd_+p), xrcv( prcv_+p))
#else
         CALL ZCOPY( sndcount, xsnd( psnd_+p), 1, xrcv( prcv_+p) , 1 )
#endif

      call MPI_waitall(nstat , request, A_OF_STATUSES, ierror)
      IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall in M_alltoall_d returns',ierror)

    ENDDO
!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)
!----------------------------------------------------------------------
#endif
#endif
!----------------------------------------------------------------------

      RETURN
      END

!----------------------------------------------------------------------
!
! M_alltoall_i: integer routine which maps directly onto
!  MPI_alltoallv
!  this is used only once by VASP 
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_i(COMM, xsnd, psnd, nsnd, xrcv, prcv, nrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER xsnd(*)           ! send buffer
      INTEGER xrcv(*)           ! receive buffer

      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
      INTEGER nsnd(COMM%NCPU+1) ! number of data send to each node
      INTEGER nrcv(COMM%NCPU+1) ! number of data recv from each node

! local data
      INTEGER ierror

      call MPI_alltoallv( xsnd(1), nsnd(1), psnd(1), MPI_integer, &
                          xrcv(1), nrcv(1), prcv(1), MPI_integer, &
                          COMM%MPI_COMM, ierror )

      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: imexch MPI_alltoallv returns',ierror)

      RETURN
      END

!----------------------------------------------------------------------
!
! on the T3E M_alltoallv_z can use shmemput instead of MPI
! M_alltoallv_z_prepare prepares the additional array rprcv
! which contains the remote locations
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_raddr(COMM, prcv, rprcv)
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM

      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
      INTEGER rprcv(COMM%NCPU)  ! remote location of data
! local variable
      INTEGER ierror

      ! only one simple MPI_alltoall is required
      CALL MPI_alltoall( prcv(1),  1, MPI_integer, &
                         rprcv(1), 1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z_prepare',ierror)


      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoallv_simple requires as only input the number of data nsnd send
! from each node to each other node
! it assumes a continous data arrangement on sender and receiver
! and sets up the arrays which are required for MPI_alltoall 
! nrcv, psnd, prcv
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoallv_simple(COMM, nsnd, nrcv, psnd, prcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
! input
      INTEGER nsnd(COMM%NCPU)   ! number of data send to each node
! output
      INTEGER psnd(COMM%NCPU+1) ! location of data in send buffer (0 based)
      INTEGER nrcv(COMM%NCPU)   ! number of data recv from each node
      INTEGER prcv(COMM%NCPU+1) ! location of data in recv buffer (0 based)
! local variable
      INTEGER ierror,i

      ! only one simple MPI_alltoall is required
      ! to find number of received data on each node
      CALL MPI_alltoall( nsnd(1),  1, MPI_integer, &
                         nrcv(1),  1, MPI_integer, &
                         COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoallv_z_prepare',ierror)
      ! now set the locations assuming linear arrangement of data

      psnd(1)=0
      prcv(1)=0

      DO i=1,COMM%NCPU
        psnd(i+1)=psnd(i)+nsnd(i)
        prcv(i+1)=prcv(i)+nrcv(i)
      ENDDO

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoall_d: complex and real global exchange routine
!     redistributes an array from distribution over bands to
!     distribution over coefficient (or vice versa)
!     original distribution            final distribution
!     |  1  |  2  |  3  |  4  |      |  1  |  1  |  1  |  1  |
!     |  1  |  2  |  3  |  4  |      |  2  |  2  |  2  |  2  |
!     |  1  |  2  |  3  |  4  | <->  |  3  |  3  |  3  |  3  |
!     |  1  |  2  |  3  |  4  |      |  4  |  4  |  4  |  4  |
!
!     xsnd is the array to be redistributed (having n elements)
!     xrcv is the result array with n/NCPU elements received
!     from each processore
!
!     mind that only (n/NCPU) *NCPU data are exchanged
!     it is the responsability of the user to guarantee that n is
!     correct
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st,i 
!----------------------------------------------------------------------
#if defined(T3D_SMA)
!----------------------------------------------------------------------

      INTEGER i, j, inode, jnode, iold
      INTEGER,EXTERNAL ::  shmem_put, shmem_udcflush, GET_D_STREAM

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU
      inode = COMM%NODE_ME-1

      ! do a local memory-memory copy (if inode == jnode)
#if defined (use_fastbcopy)
      CALL FASTBCOPY8( sndcount, xsnd(inode*sndcount + 1), xrcv(inode*sndcount + 1))
#else
      CALL DCOPY( sndcount, xsnd(inode*sndcount + 1), 1, xrcv(inode*sndcount + 1) , 1 )
#endif
      ! syncronize all nodes involved in communication

      iold = GET_D_STREAM()           ! read current streaming state
      CALL QUIET_D_STREAM()           ! disable  streaming
      CALL MPI_barrier( COMM%MPI_COMM, ierror )

      DO i = 0, COMM%NCPU-1

        ! send data to jnode
        jnode = MOD(IEOR(i,(COMM%NODE_ME-1)), COMM%NCPU)

        IF (jnode /= inode) THEN
        ! put the data
            shmem_st = shmem_put(xrcv(inode*sndcount + 1), &
                 xsnd(jnode*sndcount + 1), &
                 sndcount, COMM%hid(jnode))
        END IF


      ENDDO

      CALL MPI_barrier( COMM%MPI_COMM, ierror ) ! wait for other nodes
      CALL SET_D_STREAM(iold)                   ! restore streaming state
      shmem_st = shmem_udcflush()               ! flush the data cache

!----------------------------------------------------------------------
#else
!----------------------------------------------------------------------

      INTEGER, SAVE :: tag=201
      INTEGER       :: in
      INTEGER       :: request((COMM%NCPU-1)*2)
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,(COMM%NCPU-1)*2)
      INTEGER, PARAMETER :: max_=MPI_BLOCK
      INTEGER       :: block, p, sndcount_
      INTEGER       :: actual_proc_group, com_proc_group, &
           proc_group, group_base, i_in_group, irequests

! quick return if possible
      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

!----------------------------------------------------------------------
#if defined (use_collective) || defined (use_symmetric_heap)
!----------------------------------------------------------------------

      call MPI_alltoall( xsnd(1), sndcount, MPI_double_precision, &
     &                   xrcv(1), rcvcount, MPI_double_precision, &
     &                   COMM%MPI_COMM, ierror )
      IF ( ierror /= MPI_success ) &
         CALL M_stop_ierr('ERROR: M_alltoall_d returns',ierror)

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)

!----------------------------------------------------------------------
#else
#ifndef avoid_async
#ifndef PROC_GROUP
!
!  this version  throws all data in one go onto all nodes
!  then waits for all sends and receives to finish
!----------------------------------------------------------------------
      ! initiate the receive on all nodes
      ! local copy has already been done on each node handle remaining  NCPU-1 packages
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

      DO in = 1, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network
         
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         call MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(in), ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)
      ENDDO

      ! initiate the send on all nodes
      DO in = 1, COMM%NCPU-1
         ! recommended for T3E 
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based

         call MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(in+ COMM%NCPU-1), ierror )
         IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)
      ENDDO

      ! local memory copy for data kept on the local node
      ! overlaps with communication
#if defined (use_fastbcopy)
      CALL FASTBCOPY8( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), xrcv((COMM%NODE_ME-1)*sndcount + p))
#else
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )
#endif
        
      call MPI_waitall((COMM%NCPU-1)*2, request, A_OF_STATUSES, ierror)
      IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall in M_alltoall_d returns',ierror)
      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)
!----------------------------------------------------------------------
#else
!
!  this version  groups the processors into groups
!  with PROC_GROUP members each
!  communication is first done within one group and
!  then between any two groups
!----------------------------------------------------------------------
      actual_proc_group=MIN(PROC_GROUP,COMM%NCPU)
      ! comm_proc_group is allways larger than COMM%NCPU
      com_proc_group =((COMM%NCPU+actual_proc_group-1)/ actual_proc_group)*actual_proc_group
      
!      WRITE(*,*) COMM%NODE_ME, actual_proc_group,com_proc_group

      !
      ! blocking on the nodes
      !
      DO proc_group=0, com_proc_group -1, actual_proc_group
      !
      ! blocking on the data 
      ! 
      DO block = 0, sndcount-1, max_
      sndcount_= MIN(max_, sndcount-block)
      p        = 1 + block   ! pointer to the current block base address

         irequests=0 ! counts the number of issued requests

         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! within the group send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node from which to receive (i is zero based)

           i = MOD(i_in_group+group_base+proc_group , com_proc_group) 
!           WRITE(*,*) 'receive',i,COMM%NODE_ME-1

           ! local copies are handled below, and take care of sends from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              call MPI_irecv( xrcv(i*sndcount + p), sndcount_, MPI_double_precision, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              IF ( ierror /= MPI_success ) &
                 CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)
           ENDIF
         ENDDO

      ! initiate the send on all nodes
         DO in = 0, actual_proc_group-1
           ! base index of the present group
           group_base = ((COMM%NODE_ME-1)/actual_proc_group)*actual_proc_group

           ! send to node i_in_group
           i_in_group = MOD(in+COMM%NODE_ME-1 , actual_proc_group)

           ! actual node to which to send
           i = MOD(i_in_group+group_base-proc_group+com_proc_group , com_proc_group) 

!           WRITE(*,*) 'send',COMM%NODE_ME-1,i

           ! local receives are handled below, also take care of receives from smaller groups
           IF ( i/=COMM%NODE_ME-1 .AND. i<COMM%NCPU ) THEN
              irequests=irequests+1
              call MPI_isend( xsnd(i*sndcount + p), sndcount_, MPI_double_precision, &
                   i, tag,  COMM%MPI_COMM, request(irequests), ierror )
              IF ( ierror /= MPI_success ) &
                  CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)
           ENDIF
         ENDDO

      call MPI_waitall(irequests, request, A_OF_STATUSES, ierror)
      IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall in M_alltoall_d returns',ierror)

      ! local memory copy for data kept on the local node
      ! overlaps with communication
#if defined (use_fastbcopy)
      CALL FASTBCOPY8( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), xrcv((COMM%NODE_ME-1)*sndcount + p))
#else
      CALL DCOPY( sndcount_, xsnd((COMM%NODE_ME-1)*sndcount + p), 1, xrcv((COMM%NODE_ME-1)*sndcount + p) , 1 )
#endif
        
      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)

      ENDDO

!      CALL MPI_barrier( COMM%MPI_COMM, ierror )
!      IF ( ierror /= MPI_success ) &
!         CALL M_stop_ierr('Initpm: Error in MPI_barrier',ierror)
!----------------------------------------------------------------------
#endif
#else
! more moderate version copies data between two nodes and then blocks
! not so bad if we have full duplex point to point connections
!----------------------------------------------------------------------
      ! local memory copy for data kept on the local node
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

      ! initiate the receive and send on all nodes
      ! local copy has already been done on each node 
      ! handle remaining  NCPU-1 packages
      DO in = 1, COMM%NCPU-1
         ! receive from: (own node id) - in
         i = MOD(-in+COMM%NODE_ME-1 +COMM%NCPU, COMM%NCPU)  ! i zero based
         call MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(1), ierror )
         IF ( ierror /= MPI_success ) &
            CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)

         !
         ! sent to: (own node id) + in
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)  ! i zero based
         call MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, request(2), ierror )
         IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d irecv returns',ierror)

         call MPI_waitall(2, request, A_OF_STATUSES, ierror)
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall in M_alltoall_d returns',ierror)

!         call MPI_barrier( COMM%MPI_COMM, ierror )
!            IF ( ierror /= MPI_success ) &
!              CALL M_stop_ierr('ERROR: M_barrier in M_alltoall_d returns',ierror)

      ENDDO
!----------------------------------------------------------------------
#endif
#endif
!----------------------------------------------------------------------

#endif
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_alltoall_z: uses M_alltoall_d with twice as many elements
!
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_z(COMM, n, xsnd, xrcv )
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      COMPLEX(q) xsnd(n), xrcv(n)

      CALL M_alltoall_d(COMM, n*2, xsnd, xrcv)
      END SUBROUTINE

!----------------------------------------------------------------------
!
! z/M_alltoall_d: complex and real global exchange routine
!     redistributes an array from distribution over bands to
!     distribution over coefficient (or vice versa)
!     original distribution            final distribution
!     |  1  |  2  |  3  |  4  |      |  1  |  1  |  1  |  1  |
!     |  1  |  2  |  3  |  4  |      |  2  |  2  |  2  |  2  |
!     |  1  |  2  |  3  |  4  | <->  |  3  |  3  |  3  |  3  |
!     |  1  |  2  |  3  |  4  |      |  4  |  4  |  4  |  4  |
!
!     xsnd is the array to be redistributed (having n elements)
!     xrcv is the result array with n/NCPU elements received
!     from each processore
!
!     mind that only (n/NCPU) *NCPU data are exchanged
!     it is the responsability of the user to guarantee that n is
!     correct
!----------------------------------------------------------------------

      SUBROUTINE M_alltoall_d_async(COMM, n, xsnd, xrcv, tag, srequest, rrequest )
      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER n
      INTEGER tag
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      REAL(q) xsnd(n), xrcv(n)
! how many data are send / and received
      INTEGER sndcount, rcvcount, ierror, shmem_st
      INTEGER i,j,in

      IF (COMM%NCPU == 1) RETURN

      sndcount = n/ COMM%NCPU
      rcvcount = n/ COMM%NCPU

      ! local memory copy for data kept on the local node
      DO i = 1,sndcount
         xrcv((COMM%NODE_ME-1)*sndcount + i) = xsnd((COMM%NODE_ME-1)*sndcount + i)
      ENDDO

      ! initiate send and receive on all nodes
      ! local copy has already been done each node send NCPU-1 packages
      j=1
      DO in = 0, COMM%NCPU-1
         ! send to node in + own node id
         ! such a construct should allow for an efficient use of the network
         
         i = MOD(in+COMM%NODE_ME-1 , COMM%NCPU)
         IF ( COMM%NODE_ME-1 /= i) THEN
            call MPI_isend( xsnd(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, srequest(j), ierror )
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d_async MPI_alltoall returns',ierror)
            call MPI_irecv( xrcv(i*sndcount + 1), sndcount, MPI_double_precision, &
     &                  i, tag,  COMM%MPI_COMM, rrequest(j), ierror )
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_alltoall_d_async MPI_alltoall returns',ierror)
            j=j+1
         ENDIF
      ENDDO
        
      END SUBROUTINE

      SUBROUTINE M_alltoall_wait(COMM, srequest, rrequest )

      USE mpimy
      IMPLICIT NONE
      INCLUDE "pm.inc"

      TYPE(communic) COMM
      INTEGER srequest(COMM%NCPU), rrequest(COMM%NCPU)
      INTEGER ierror
      INTEGER A_OF_STATUSES(MPI_STATUS_SIZE,COMM%NCPU)

      ! wait for the NCPU-1 outstanding packages

      call MPI_waitall(COMM%NCPU-1, srequest, A_OF_STATUSES, ierror)
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall returns',ierror)
      call MPI_waitall(COMM%NCPU-1, rrequest, A_OF_STATUSES, ierror)
            IF ( ierror /= MPI_success ) &
              CALL M_stop_ierr('ERROR: M_waitall returns',ierror)

      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sumf_d: performs a fast global sum on n doubles in
! vector vec (algorithm by Kresse Georg)
!
! uses complete interchange algorithm (my own invention, but I guess
!  some people must know it)
! exchange data between nodes, sum locally and
! interchange back, this algorithm is faster than typical MPI based
! algorithms (on 8 nodes under MPICH a factor 4)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n,ncount,nsummed,ndo,i,j, info, n_,mmax
      REAL(q) vec(n)
!----------------------------------------------------------------------
#if defined(T3D_SMA)
!----------------------------------------------------------------------
      INTEGER MALLOC_DONE
      INTEGER, EXTERNAL :: ISHM_CHECK
      COMMON /SHM/ MALLOC_DONE, PBUF
      POINTER ( PBUF, vec_inter )
      REAL(q) :: vec_inter(n/COMM%NCPU*COMM%NCPU)
      INTEGER :: max_=n/COMM%NCPU

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN
      ! do we have sufficient shm workspace to use fast interchange algorithm
      ! no use conventional M_sumb_d
      IF (ISHM_CHECK(n) == 0) THEN
         CALL M_sumb_d(COMM, vec, n)
         RETURN
      ENDIF
!----------------------------------------------------------------------
#else
!----------------------------------------------------------------------
      REAL(q), ALLOCATABLE :: vec_inter(:)
    ! maximum work space for quick sum
!
! maximum communication blocks
! too large blocks are slower on the Pentium architecture
! probably due to caching
!
      INTEGER, PARAMETER :: max_=MPI_BLOCK

    ! quick return if possible
      IF (COMM%NCPU == 1) RETURN
      
      mmax=MIN(n/COMM%NCPU,max_)
      ALLOCATE(vec_inter(mmax*COMM%NCPU))
!----------------------------------------------------------------------
#endif
!----------------------------------------------------------------------

      nsummed=0
      n_=n/COMM%NCPU

      DO ndo=0,n_-1,mmax
     ! forward exchange
         ncount =MIN(mmax,n_-ndo)
         nsummed=nsummed+ncount*COMM%NCPU

         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
        ! sum localy
         DO i=2, COMM%NCPU
           CALL DAXPY(ncount, 1.0_q, vec_inter(1+(i-1)*ncount), 1, vec_inter(1), 1)
         ENDDO
     ! replicate data (will be send to each proc)
         DO i=1, COMM%NCPU
            DO j=1,ncount
               vec(ndo*COMM%NCPU+j+(i-1)*ncount) = vec_inter(j)
            ENDDO
         ENDDO
     ! backward exchange
         CALL M_alltoall_d(COMM, ncount*COMM%NCPU, vec(ndo*COMM%NCPU+1), vec_inter(1))
#if defined (use_fastbcopy)
         CALL FASTBCOPY8( ncount*COMM%NCPU, vec_inter(1), vec(ndo*COMM%NCPU+1))
#else
         CALL DCOPY( ncount*COMM%NCPU, vec_inter(1), 1, vec(ndo*COMM%NCPU+1), 1 )
#endif
      ENDDO
      
     ! that should be it
      IF (n_*COMM%NCPU /= nsummed) THEN
         WRITE(0,*) 'internal error in M_sumf_g',n_,nsummed
         STOP
      ENDIF

      IF (n-nsummed /= 0 ) &
        CALL M_sumb_d(COMM, vec(nsummed+1), n-nsummed)

#if defined(T3D_SMA)
     ! nup nothing to do here
#else
      DEALLOCATE(vec_inter)
#endif
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sumf_g: performs a fast global sum on n complex in
! vector 'vec' (see above)
!
!----------------------------------------------------------------------

      SUBROUTINE M_sumf_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)
      CALL M_sumf_d(COMM, vec, 2*n)
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_z: performs a sum on n double complex numbers
!  it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_z(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      IF ( 2*n>MPI_BLOCK) THEN
         CALL M_sumf_d(COMM, vec, 2*n)
      ELSE
         CALL M_sumb_d(COMM, vec, 2*n)
      ENDIF
            
      END SUBROUTINE

!----------------------------------------------------------------------
!
! M_sum_d:  performs a sum on n double prec numbers
!  it uses either sumb_d or sumf_d
!
!----------------------------------------------------------------------

      SUBROUTINE M_sum_d(COMM, vec, n)
      USE mpimy
      IMPLICIT NONE

      TYPE(communic) COMM
      INTEGER n
      REAL(q) vec(n)

      IF ( n>MPI_BLOCK) THEN
         CALL M_sumf_d(COMM, vec, n)
      ELSE
         CALL M_sumb_d(COMM, vec, n)
      ENDIF
            
      END SUBROUTINE

#endif
