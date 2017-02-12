#include "symbol.inc"
!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! Module containing wrapper for scaLAPACK
! written by Gilles de Wijs (gD) and Georg Kresse (gK)
! modified to run on any number of nodes by Dario Alfe
!
!=======================================================================

#ifndef scaLAPACK

 MODULE scala
   USE prec
   USE mpimy
   LOGICAL, PUBLIC :: LscaLAPACK = .FALSE.
   LOGICAL, PUBLIC :: LscaLU     = .FALSE.  ! use parallel LU  decomposition


 CONTAINS

      SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN,N)
      TYPE (communic) COMM
      INTEGER N              ! NxN matrix to be distributed
      GDEF    Amatin(N,N)    ! input/output matrix
      END SUBROUTINE

      SUBROUTINE pDSSYEX_ZHEEVX(COMM,AMATIN,W,N)
      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      GDEF    AMATIN(N,N)     ! input/output matrix
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

#else

 MODULE scala
    USE prec
    USE mpimy
    LOGICAL, PUBLIC :: LscaLAPACK = .TRUE.
    LOGICAL, PUBLIC :: LscaLU     = .FALSE.  ! use parallel LU  decomposition

!
! customize if required
!
      INTEGER,PRIVATE :: NB=80        ! blocking factor for distribution of matrices
      ! P4 optimal, larger values were even slightly better (160) but still
      ! slower on a Gigabit cluster than ZHEEVX

      ! customization of matrix diagonalization

      INTEGER,PRIVATE :: NCLUST=16    ! maximum cluster of eigenvector
      REAL(q),PRIVATE :: ABSTOL=1e-10 ! specifies eigenvector orthogonality tolerance
      REAL(q),PRIVATE :: ORFAC=-1.e0  ! controls reorthogonalisation of eigenvectors

      !  see manpage for pSSYEVX pZHEEVX:
      !    NCLUST is the number of eigenvalues in the largest cluster,
      !    where a cluster is defined as a set of close eigenvalues:
      !    {W(K),...,W(K+NCLUST-1)|W(J+1)<= W(J)+orfac*norm(A)}
!
! end customization
!
      LOGICAL,SAVE,PRIVATE :: LINIT=.FALSE.  ! tells whether MODULE is initialized
!
! Cray uses different layout for descriptors in the old
! scaLAPACK version, if you have access only to this version
! please change following lines accordingly
!#ifdef T3D_SCA
!      INTEGER,PARAMETER,PRIVATE :: CSRC_=6, CTXT_=7, DLEN_=8, LLD_=8, &
!                           MB_=3, M_=1, NB_=4, N_=2, RSRC_=5
!#else
      INTEGER,PARAMETER,PRIVATE :: BLOCK_CYCLIC_2D=1, DT_=1, &
                           CSRC_ =8, CTXT_=2, DLEN_=9, LLD_=9, &
                           MB_=5, M_=3, NB_=6, N_=4, RSRC_=7
!#endif
!
! stuff filled in by INIT_scala
!
      INTEGER,SAVE,PRIVATE :: IAM            ! processor number
      INTEGER,SAVE,PRIVATE :: NPROCS         ! number of PEs
      INTEGER,SAVE,PRIVATE :: ICTXT          ! context handle of grid
      INTEGER,SAVE,PRIVATE :: NPROW,NPCOL    ! processor grid dimensions
      INTEGER,SAVE,PRIVATE :: MYROW,MYCOL    ! processor coordinates in ps. grid
      INTEGER,SAVE,PRIVATE :: NP             ! number of rows on the processor
      INTEGER,SAVE,PRIVATE :: NQ             ! number of cols on the processor
      INTEGER,SAVE,PRIVATE :: DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER,SAVE,PRIVATE :: DESCZ( DLEN_ ) ! distributed matrix descriptor array
      INTEGER,SAVE,PRIVATE :: LWWORK,LIWORK,LRWORK,MALLOCPQ ! sizes of scalapack workarrays
      INTEGER,SAVE,PRIVATE :: NBANDS_        ! number of nbands, remembered for check

#ifdef gammareal
      INTEGER, PARAMETER,PRIVATE :: MCOMP=1
#else
      INTEGER, PARAMETER,PRIVATE :: MCOMP=2
#endif

 CONTAINS

!=======================================================================
!
! this subroutine determines the processor grid
!
!=======================================================================

!=======================================================================
!
! this subroutine does:
!   1. calls INIT_scala
!   2. calls the setup of distributed the distributed matrix
!   3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!   4. calls reconstruction of distributed data into patched matrix
! note: the sum over processors of the patched matrix is not realised
!       in this module
!
!=======================================================================


      SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN,N) ! ad hoc GILLES
      IMPLICIT NONE

      TYPE (communic) COMM                     ! ad hoc GILLES
      INTEGER N              ! NxN matrix to be distributed
      GDEF    Amatin(N,N)    ! input/output matrix
! local variables
      INTEGER MINB,MAXD,MINP,MAXB,LENGTH,LENGTHW,INFO
      INTEGER ierror         ! ad hoc GILLES
      REAL(q) ANOTHING
      INTEGER INOTHING
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP

#ifndef T3D_SCA
      GDEF, ALLOCATABLE ::  A(:)
#else
      POINTER ( A_P, A );         GDEF :: A(N*N)
      POINTER ( Z_P, Z );         GDEF :: Z(N*N)
      POINTER ( WWORK_P, WWORK ); GDEF :: WWORK(N*N)
      POINTER ( RWORK_P, RWORK ); REAL(q) :: RWORK(N*N)
      POINTER ( IWORK_P, IWORK ); INTEGER IWORK(N*N)
      POINTER ( WORK_P, WORK ) ;  INTEGER WORK(N)
      COMMON /SCASHM/ A_P,Z_P,WWORK_P,RWORK_P,IWORK_P,WORK_P
#endif

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      CALL INIT_scala(N, INFO)

!-----------------------------------------------------------------------
! allocate workarray (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------

    #ifndef T3D_SCA
      ALLOCATE(A(MALLOCPQ))
    #endif

!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( MYROW < NPROW .AND. MYCOL < NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,N,A,DESCA)
!       call scalapack routine
        CALL BLACS_GRIDINFO( ICTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

        INFO=0
      #ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, DESCA, INFO )
      #else
        CALL PZPOTRF( 'U', N, A, 1, 1, DESCA, INFO )
      #endif

        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_TRTRI, POTRF, INFO:',INFO
          WRITE(*,*) 'STOP'
          STOP
        ENDIF
      #ifdef gammareal
         CALL PDTRTRI &
      #else
         CALL PZTRTRI &
      #endif
     &    ('U', 'N', N, A, 1, 1, DESCA, INFO)
        IF (INFO.NE.0) THEN
          WRITE(*,*) 'pPOTRF_TRTRI, TRTRI, INFO:',INFO
          STOP
        ENDIF

      CALL RECON(AMATIN,N,A,DESCA)
   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      #ifndef T3D_SCA
        DEALLOCATE(A)
      #endif

      RETURN
      END SUBROUTINE


!=======================================================================
!
! call to pZHEEVX respectively pDSSYEX
! i.e. diagonalization of an symmetric or hermitian matrix
!
!=======================================================================

      SUBROUTINE pDSSYEX_ZHEEVX(COMM,AMATIN,W,N)
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      GDEF    AMATIN(N,N)     ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.
      INTEGER MINB,MAXD,MINP,MAXB,LENGTH,LENGTHW,INFO
      INTEGER M,NZ,I1,NNP,NN,NEIG,NP0,NQ0,CLUSTERSIZE
      INTEGER I2

#ifndef T3D_SCA
      GDEF, ALLOCATABLE    ::  A(:)     ! A matrix to be diagonalized
      GDEF, ALLOCATABLE    ::  Z(:)     ! Z eigenvector matrix
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
#else
      POINTER ( A_P, A );         GDEF :: A(N*N)
      POINTER ( Z_P, Z );         GDEF :: Z(N*N)
      POINTER ( WWORK_P, WWORK ); GDEF :: WWORK(N*N)
      POINTER ( RWORK_P, RWORK ); REAL(q) :: RWORK(N*N)
      POINTER ( IWORK_P, IWORK ); INTEGER IWORK(N*N)
      POINTER ( WORK_P, WORK ) ;  INTEGER WORK(N)
      COMMON /SCASHM/ A_P,Z_P,WWORK_P,RWORK_P,IWORK_P,WORK_P
#endif
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(N, INFO)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((NPROW*NPCOL) > COMM%NCPU) THEN
        WRITE(*,*) "pDSSYEX_ZHEEVX: too many processors ", &
                     NPROW*NPCOL,COMM%NCPU
        STOP
      ENDIF
#ifndef T3D_SCA
!-----------------------------------------------------------------------
! allocation (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
        ALLOCATE(A(MALLOCPQ))
        ALLOCATE(Z(MALLOCPQ))
        ALLOCATE(IWORK(LIWORK))
        ALLOCATE(WWORK(LWWORK))
        ALLOCATE(RWORK(LRWORK))
#endif

!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
   calc: IF( MYROW < NPROW .AND. MYCOL < NPCOL ) THEN

!     get parts of the matrix into the local arrays
      CALL DISTRI(AMATIN,N,A,DESCA)
      INFO=0

!     call scalapack routine
      #ifdef gammareal
        CALL PDSYEVX( 'V', 'A', 'U', N, A, 1, 1, DESCA, ZERO, ZERO, 13,  &
                    -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCZ, WWORK, &
                    LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
      #else
        CALL PZHEEVX( 'V', 'A', 'U', N, A, 1, 1, DESCA, ZERO, ZERO, 13,  &
                    -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCZ, &
                    WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )
      #endif
        IF (M.NE.N) THEN
          WRITE (*,*) 'LWWORK',LWWORK,LRWORK,LIWORK,DESCA(M_)
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvalues found",M,N
          STOP
        ENDIF
        IF (NZ.NE.N) THEN
          WRITE(*,*) "ERROR in subspace rotation PSSYEVX: not enough eigenvectors computed",M,N
          STOP
        ENDIF
        DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
            WRITE(*,*) "ERROR in subspace rotation PSSYEVX: I2,IFAIL= ",I2,IFAIL(I2)
            STOP
          ENDIF
        ENDDO

! cluster check, gap check ?
!        IF (INFO.NE.0) THEN
!          do_io WRITE(0,*) "WARNING in subspace rotation, PSSYEVX: INFO=",INFO
!        ENDIF
!     reconstruct from the Z matrix
      CALL RECON(AMATIN,N,Z,DESCZ)

   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      #ifndef T3D_SCA
        DEALLOCATE(A)
        DEALLOCATE(Z)
        DEALLOCATE(IWORK)
        DEALLOCATE(WWORK)
        DEALLOCATE(RWORK)
      #endif

      RETURN
      END SUBROUTINE


!=======================================================================
!
! setup a distributed matrix
!
!=======================================================================

      SUBROUTINE DISTRI(AMATIN,N,A,DESCA)
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array

      GDEF    AMATIN(N,N),A(*)
      INTEGER N,NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW,MYCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,desca(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,desca(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
      END SUBROUTINE


!=======================================================================
!
! merge distributed matrix into one large matrix
!
!=======================================================================

      SUBROUTINE RECON(amatin,n,a,desca)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array

      GDEF    AMATIN(N,N),A(*)
      INTEGER N,NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW,MYCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(desca(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,desca(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,desca(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!=======================================================================
!
! messy scalapack initialization, calculate all required workspace
! and set up matrix descriptors
!
!=======================================================================
      SUBROUTINE INIT_scala(N, MALLOC)

      USE main_mpi            ! to get world communicator
      IMPLICIT NONE

!     TYPE(communic) COMM
!     TYPE(communic) COMM_WORLD

      INTEGER IPDIM,NGUESS,NPROCS,N,INFO
      INTEGER MINB,MAXB,MINP,MAXD,MALLOC,LENGTH,LENGTHW,LENGTHW2
      INTEGER NP0,NQ0,NN,NNP,NEIG
      INTEGER ierror

#ifdef T3D_SCA
      POINTER ( A_P, A );         GDEF :: A(N*N)
      POINTER ( Z_P, Z );         GDEF :: Z(N*N)
      POINTER ( WWORK_P, WWORK ); GDEF :: WWORK(N*N)
      POINTER ( RWORK_P, RWORK ); REAL(q) :: RWORK(N*N)
      POINTER ( IWORK_P, IWORK ); INTEGER IWORK(N*N)
      POINTER ( WORK_P, WORK ) ;  INTEGER WORK(N)

      COMMON /SCASHM/ A_P,Z_P,WWORK_P,RWORK_P,IWORK_P,WORK_P
#endif

      INTEGER,INTRINSIC :: MAX,MIN
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      integer :: factors(4), pwr(4)        
      data factors /2, 3, 5, 7/                      ! a few factors... 
      integer :: nnpcol, nnprow, mr, ndum, nmin, i

      EXTERNAL BLACS_PINFO
      EXTERNAL DESCINIT

      IF ((LINIT).AND.(NBANDS_.NE.N)) THEN
        WRITE (*,*) 'INIT_scala: NUMBER OF BANDS IS NOT CONSTANT,'
        WRITE (*,*) 'NOT ALLOWED ON T3D, STOP'
        STOP
      ENDIF

      IF (.NOT. LINIT) THEN

      LINIT=.TRUE.
      NBANDS_=N

! number of nodes (for one image)
      NPROCS = COMM%NCPU
! determine processor grid, according to naive guess
      pwr=0
      mr = nprocs
      nprow=1; npcol=1
      do i=1,4
         do
           if( mod(mr,factors(i)) == 0 ) then
              pwr(i) = pwr(i) + 1
              mr = mr/factors(i)
           else
              exit
           endif
         enddo
         nnprow = pwr(i)/2
         nnpcol = pwr(i) - nnprow
         if(mod(i,2)==0) then         ! to make the grid more uniformly distributed (hopefully)
           ndum = nnpcol
           nnpcol = nnprow
           nnprow = ndum
         endif
         nprow = nprow*factors(i)**nnprow
         npcol = npcol*factors(i)**nnpcol
      enddo
      if( nprow*npcol < nprocs )then  ! final adjustment to the total number of procs
         nmin = MIN(nprow,npcol)
         ndum = nprocs/nprow/npcol
         if(nprow==nmin)then          ! again, hoping that in this way the grid is more uniform
           nprow = nprow * ndum       ! you might think to a better algorithm ( I dont think it
         else                         ! really matters )
           npcol = npcol * ndum
         endif
      endif

      DWRITE "NPROW",NPROW," NPCOL",NPCOL

! make processor grids

      CALL MPI_barrier( COMM_WORLD%MPI_COMM, ierror )
      CALL PROCMAP( ICTXT, 1, NPROW, NPCOL)

! calculate local size of matrices

      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      DWRITE "NPROW",NPROW," NPCOL",NPCOL,"MYROW ",MYROW," MYCOL",MYCOL

      NP = NUMROC(N,NB,MYROW,0,NPROW)   ! get number of rows on proc
      NQ = NUMROC(N,NB,MYCOL,0,NPCOL)   ! get number of cols on proc
      DWRITE "SETUP,NP,NQ",NP,NQ," MB,NB",desca(MB_),desca(NB_)," lld",desca(lld_)

      CALL DESCINIT(DESCA,N,N,NB,NB,0,0,ICTXT,MAX( 1, NP ),INFO ) ! setup descriptor
      IF (INFO.NE.0) THEN
        WRITE(*,*) 'pPOTRF_TRTRI, DESCINIT, INFO: ', INFO
        STOP
      ENDIF
      CALL DESCINIT(DESCZ,N,N,NB,NB,0,0,ICTXT,MAX( 1, NP ),INFO ) ! setup descriptor
      IF (INFO.NE.0) THEN
        WRITE(*,*) 'pSSYEZS_ETC, DESCINIT, A, INFO: ', INFO
        STOP
      ENDIF

! calculate the size of the BLACS workspace for CRAY and allocate

   ! first for pPOTRF_TRTRI

      MINB=MIN(DESCA(MB_),DESCA(NB_))
      MAXB=MAX(DESCA(MB_),DESCA(NB_))
      MINP=MIN(NPROW,NPCOL)
      MAXD=N
      ! pspotrf length=4*maxb*(((maxd/minb)/minp)*maxb+maxb)*8  bytes
      ! pcpotrf length=4*maxb*(((maxd/minb)/minp)*maxb+maxb)*16 bytes
      ! pstrtri length=2*maxb*(((maxd/minb)/minp)*maxb+maxb)*8  bytes
      !     pctrtri length=2*maxb*(((maxd/minb)/minp)*maxb+maxb)*16
      LENGTH=4*MAXB*(((MAXD/MINB)/MINP)*MAXB+8*MAXB)*8*MCOMP
      LENGTHW=LENGTH/8    ! length in words

   ! next for pDSSYEX_ZHEEVX

      MINB=MIN(DESCA(MB_),DESCA(NB_))
      MAXB=MAX(DESCA(MB_),DESCA(NB_))
      MINP=MIN(NPROW,NPCOL)
      MAXD=N
      LENGTH=4*MAXB*(((MAXD/MINB)/MINP)*MAXB+8*MAXB)*8*MCOMP
      LENGTHW2=LENGTH/8             ! length in words

   ! now determine the largest, allocate and initialise workspace

      LENGTHW=MAX(LENGTHW,LENGTHW2)
      DWRITE "lengthw", LENGTHW

#ifdef T3D_SCA
      CALLMPI ( M_max_i(COMM_WORLD, LENGTHW, 1))           ! just in case...
      CALL SHPALLOC( WORK_P, LENGTHW , INFO, 1)
      IF (INFO.NE.0) THEN
        WRITE (*,*) 'INIT_SCALA: ERROR: SHPALLOC A, INFO:',INFO
        STOP
      ENDIF

      CALL INITBUFF(WORK,LENGTH)  ! hand over workspace to BLACS
#endif

! calculate scalapack workspace and allocate, pDSSYEX_ZHEEVX only

   ! iwork
      NNP=MAX(N,NPROW*NPCOL+1,4)
      LIWORK=6*NNP
   ! wwork
      NN = MAX( N, NB, 2 )
      NEIG = N                ! number of eigenvectors requested
      IF ((DESCA(MB_).NE.NB).OR.(DESCA(NB_).NE.NB) .OR. &
          (DESCZ(MB_).NE.NB).OR.(DESCZ(NB_).NE.NB) .OR. &
          (DESCA(RSRC_).NE.0) .OR.(DESCA( CSRC_).NE.0) .OR.  &
          (DESCZ(RSRC_).NE.0) .OR.(DESCZ( CSRC_).NE.0)) THEN
        WRITE(*,*) 'INIT_SCALA: pSSYEZS_ETC, ERROR'
        STOP
      ENDIF
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
      NQ0 = MAX( NUMROC( NEIG, NB, 0, 0, NPCOL ), NB )
     #ifdef gammareal
      LRWORK=1
      LWWORK=5*N+MAX(5*NN,NP0*NQ0)+ICEIL(NEIG,NPROW*NPCOL)*NN+2*NB*NB+ &
         NCLUST*N
     #else
      LRWORK=4*N+MAX( 5*NN, NP0 * NQ0 )+ICEIL( NEIG, NPROW*NPCOL)*NN+ &
         NCLUST*N
      LWWORK=N + MAX((NP0 +NQ0 +NB ) * NB, 3)
     #endif

      MALLOCPQ=NP*NQ

      MALLOC  =0
#ifdef T3D_SCA
! just in case make

      CALLMPI ( M_max_i(COMM_WORLD, MALLOCPQ, 1))           ! just in case...
      CALLMPI ( M_max_i(COMM_WORLD, LWWORK, 1))           ! just in case...
      CALLMPI ( M_max_i(COMM_WORLD, LIWORK, 1))           ! just in case...
      CALLMPI ( M_max_i(COMM_WORLD, LRWORK, 1))           ! just in case...

      MALLOC=MALLOCPQ*MCOMP*2+LWWORK*MCOMP+LIWORK+LRWORK
#endif

! now hope everything is properly set up, try to return

      ENDIF

      END SUBROUTINE

      SUBROUTINE INIT_scala_T3D
      IMPLICIT NONE
#ifdef T3D_SCA
      POINTER ( A_P, A );         GDEF :: A(*)
      POINTER ( Z_P, Z );         GDEF :: Z(*)
      POINTER ( WWORK_P, WWORK ); GDEF :: WWORK(*)
      POINTER ( RWORK_P, RWORK ); REAL(q) :: RWORK(*)
      POINTER ( IWORK_P, IWORK ); INTEGER IWORK(*)
      POINTER ( WORK_P, WORK ) ;  INTEGER WORK(*)

      INTEGER MALLOC_DONE,MTOTAL
      INTEGER, EXTERNAL :: ISHM_CHECK
      COMMON /SCASHM/ A_P,Z_P,WWORK_P,RWORK_P,IWORK_P,WORK_P
      COMMON /SHM/  MALLOC_DONE, PBUF
      POINTER ( PBUF, BUF ); REAL(q) :: BUF(*)

      MTOTAL=0
      A_P=PBUF                     ; MTOTAL=MTOTAL+MALLOCPQ*MCOMP
      Z_P=A_P +MALLOCPQ*MCOMP      ; MTOTAL=MTOTAL+MALLOCPQ*MCOMP
      WWORK_P=Z_P+MALLOCPQ*MCOMP   ; MTOTAL=MTOTAL+LWWORK*MCOMP
      IWORK_P=WWORK_P+LWWORK*MCOMP ; MTOTAL=MTOTAL+LIWORK
      RWORK_P=IWORK_P+LIWORK       ; MTOTAL=MTOTAL+LRWORK

      IF ( ISHM_CHECK(MTOTAL) ==0 ) THEN
         WRITE(0,*) 'internal error INIT_scala_T3D: not enough shmem'
         STOP
      ENDIF

#endif
      END SUBROUTINE INIT_scala_T3D

!=======================================================================
!
!     this subroutine map processors onto a grid
!     and allows to performe scaLAPACK operations on a sub set of nodes
!     taken from BLACS example code (written by Clint Whaley 7/26/94)
!     modified by gD and gK
!
!=======================================================================

      SUBROUTINE PROCMAP(CONTEXT, MAPPING, NPROW, NPCOL)
      USE main_mpi     ! to get MPI communication handle
      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER :: CONTEXT, MAPPING, NPROW, NPCOL
      INTEGER :: IMAP(NPROW, NPCOL)
      INTEGER ierror
!     ..
!
!  Purpose
!  =======
!  PROCMAP maps NPROW*NPCOL onto an VASP-Communicator
!  It is a nice system independent implementation
!
!  Arguments
!  =========
!
!  CONTEXT      (output) INTEGER
!               This integer is used by the BLACS to indicate a context.
!               A context is a universe where messages exist and do not
!               interact with other contexts messages.  The context includes
!               the definition of a grid, and each processs coordinates in it.
!
!  MAPPING      (input) INTEGER
!               Way to map processes to grid.  Choices are:
!               1 : row-major natural ordering
!               2 : column-major natural ordering
!
!  NPROW        (input) INTEGER
!               The number of process rows the created grid
!               should have.
!
!  NPCOL        (input) INTEGER
!               The number of process columns the created grid
!               should have.
!
!  =====================================================================
!
      INTEGER,EXTERNAL :: BLACS_PNUM
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      INTEGER TMPCONTXT, NPROCS, I, J, K
      INTEGER BLACS_IDS(COMM%NCPU), MY_BLACS_ID
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP
!
!     See how many processes there are in the system
!     NPROCS here refers to all processors being used
      CALL BLACS_PINFO( I, NPROCS )
!
!     Temporarily map all processes into 1 x NPROCS grid
!
      CALL BLACS_GET( 0, 0, TMPCONTXT )
      CALL BLACS_GRIDINIT( TMPCONTXT, 'Row', 1, NPROCS )
!
!     If we want a row-major natural ordering
!
      CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

      MY_BLACS_ID= BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)

      BLACS_IDS = 0
      BLACS_IDS( COMM%NODE_ME ) = MY_BLACS_ID

      CALL M_sum_i ( COMM, BLACS_IDS, COMM%NCPU)

!     WRITE (*,*) 'MAPPING',NPCOL,NPROW,MY_BLACS_ID
      CALL MPI_barrier( COMM_WORLD%MPI_COMM, ierror )
      K=1
      IF (MAPPING .EQ. 1) THEN
         DO I = 1, NPROW
            DO J = 1, NPCOL
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
!
!     If we want a column-major natural ordering
!
      ELSE IF (MAPPING .EQ. 2) THEN
         DO J = 1, NPCOL
            DO I = 1, NPROW
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
      ELSE
         WRITE(*,*) 'INIT_SCALA, PROCMAP: unknown mapping, STOP'
         STOP
      END IF
!
!     Free temporary context
!
      CALL MPI_barrier( COMM_WORLD%MPI_COMM, ierror )
      CALL BLACS_GRIDEXIT(TMPCONTXT)
!      WRITE (*,*) 'MAPPING',NPCOL,NPROW
!
!     Apply the new mapping to form desired context
!
      CALL BLACS_GET( 0, 0, CONTEXT )
      CALL BLACS_GRIDMAP( CONTEXT, IMAP, NPROW, NPROW, NPCOL )

!      WRITE(*,'("local map",16I3)') ((IMAP(I,J),I=1,NPROW),J=1,NPCOL)

      RETURN
      END SUBROUTINE

END MODULE
#endif

