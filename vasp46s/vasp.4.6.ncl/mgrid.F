#include "symbol.inc"
! RCS:  $Id: mgrid.F,v 1.5 2003/06/27 13:22:19 kresse Exp kresse $

      MODULE mgrid
      USE prec
      USE mpimy
      INCLUDE "mgrid.inc"

      LOGICAL LCOMPAT   !   compatibility modus to vasp.4.4

!***********************************************************************
! this module sets up the GRID structure
! the GRID structure contains (among other things) the communication
! patterns for the 3d FFT.
! If LPLANE_WISE is set the data are distributed in real and reciprocal
! space plane by plane i.e. one processor holds all elements of
! a plane with a specific x index
! this  reduces the communication in the FFT considerably
! the default for LPLANE_WISE can be set in this file (see below),
! or using the flag LPLANE in the INCAR reader
!#define plane_wise
!***********************************************************************
    #ifdef plane_wise
      LOGICAL :: LPLANE_WISE=.TRUE.
     #else
      LOGICAL :: LPLANE_WISE=.FALSE.
     #endif
      CONTAINS

!****************  SUBROUTINE INIGRD ***********************************
!
! initialize the basic grid structure
! minimal setup
!***********************************************************************
      SUBROUTINE INILGRD(NGX,NGY,NGZ,GRID)
      USE prec
      USE base
      USE mpimy

      IMPLICIT NONE

      INTEGER NGX,NGY,NGZ
      TYPE (grid_3d) GRID

      GRID%NGX=NGX
      GRID%NGY=NGY
      GRID%NGZ=NGZ
      GRID%NPLWV =NGX*NGY*NGZ
      GRID%MPLWV =0
      ALLOCATE(GRID%LPCTX(NGX),GRID%LPCTY(NGY),GRID%LPCTZ(NGZ))
      ALLOCATE(GRID%LPCTX_(NGX),GRID%LPCTY_(NGY),GRID%LPCTZ_(NGZ))
      GRID%NGPTAR(1)=NGX
      GRID%NGPTAR(2)=NGY
      GRID%NGPTAR(3)=NGZ
      CALL INILPC(NGX,NGY,NGZ,GRID%LPCTX(1:NGX),GRID%LPCTY(1:NGY),GRID%LPCTZ(1:NGZ))
      CALL INILPC(NGX,NGY,NGZ,GRID%LPCTX_(1:NGX),GRID%LPCTY_(1:NGY),GRID%LPCTZ_(1:NGZ))
      GRID%LPCTX_(NGX/2+1) = 0
      GRID%LPCTY_(NGY/2+1) = 0
      GRID%LPCTZ_(NGZ/2+1) = 0
      RETURN
      END  SUBROUTINE


!****************  SUBROUTINE DEALLOC_GRD ******************************
!
!  deallocate a GRID strucuture
!
!***********************************************************************
      SUBROUTINE  DEALLOC_GRD(GRID)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      DEALLOCATE(GRID%RC%I2)
      DEALLOCATE(GRID%RC%I3)
      DEALLOCATE(GRID%RL%I2)
      DEALLOCATE(GRID%RL%I3)
      DEALLOCATE(GRID%RL%INDEX)

      END SUBROUTINE

!****************  SUBROUTINE INILPC ***********************************
!
! initialize the loop counters LPCTX,LPCTY,LPCTZ etc that
! label the number of the reciprocal lattice vectors in the x,y,z
! directions, respectively. for the x direction the reciprocal lattice
! vectors corresponding to the first,second,...,ngxth elements in all
! of the reciprocal lattice arrays are 0,1,..,(NGX/2),-((NGX/2-1),..,-1
! times the x reciprocal lattice vector
!
!***********************************************************************

      SUBROUTINE INILPC (NGX,NGY,NGZ,LPCTX,LPCTY,LPCTZ)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION LPCTX(NGX),LPCTY(NGY),LPCTZ(NGZ)

      DO 100 NX=1,(NGX/2)+1
        LPCTX(NX)=NX-1
 100  CONTINUE
      DO 110 NX=(NGX/2)+2,NGX
        LPCTX(NX)=NX-1-NGX
 110  CONTINUE
      DO 120 NY=1,(NGY/2)+1
        LPCTY(NY)=NY-1
 120  CONTINUE
      DO 130 NY=(NGY/2)+2,NGY
        LPCTY(NY)=NY-1-NGY
 130  CONTINUE
      DO 140 NZ=1,(NGZ/2)+1
        LPCTZ(NZ)=NZ-1
 140  CONTINUE
      DO 150 NZ=(NGZ/2)+2,NGZ
        LPCTZ(NZ)=NZ-1-NGZ
 150  CONTINUE
      RETURN
      END SUBROUTINE




!***********************************************************************
!  Generate layout for a 3d grid (optionally a subgrid of GRID)
!
! setup data distribution for a subgrid contained in a supergrid
! it is guaranteed that the new (small) grid has a data distribution that
! ensures that if one column (in reciprocal space) on the large grid
! (supergrid) is on proc x the corresponding column is also on proc x 
! on the small grid (subgrid)
! if a subgrid is generated a transition table
!  for going from one to the other grid is also created
!
! if the 3d grid has to be used also in real space (FFTs) 
! LREAL must be set to .TRUE. 
! grids in real space have the same date distribution 
! if and only if the dimension of the grids are the same
!
!***********************************************************************

!***********************************************************************
!
!  a RC grid is a grid possibly real in the direct space
!  and complex in the reciprocal space (depends on the selection in 
!   the makefile and symbol.inc)
!
!***********************************************************************

      SUBROUTINE GEN_RC_GRID(GRID)
      USE prec
      USE mpimy
      IMPLICIT NONE

      TYPE (transit)     TRANS
      TYPE (grid_3d)     GRID

      GRID%NGX_rd=xm(GRID%NGX)
      GRID%NGY_rd=GRID%NGY
      GRID%NGZ_rd=zm(GRID%NGZ)
#ifdef realmode
      GRID%LREAL =.TRUE.
#else
      GRID%LREAL =.FALSE.
#endif      
      CALL GEN_SUB_GRID_(GRID, GRID, TRANS,.TRUE.,.FALSE.)

      END SUBROUTINE


!***********************************************************************
!
!  full grid, complex in direct space
!  and complex in the reciprocal space
!
!***********************************************************************

      SUBROUTINE GEN_GRID(GRID)
      USE prec
      USE mpimy

      IMPLICIT NONE

      TYPE (transit)     TRANS
      TYPE (grid_3d)     GRID

      GRID%NGX_rd=GRID%NGX
      GRID%NGY_rd=GRID%NGY
      GRID%NGZ_rd=GRID%NGZ

      GRID%LREAL=.FALSE.

      CALL GEN_SUB_GRID_(GRID,GRID, TRANS,.TRUE.,.FALSE.)

      END SUBROUTINE

!***********************************************************************
!
!  generate RC sub grids
!
!***********************************************************************

      SUBROUTINE GEN_RC_SUB_GRID(GRID_SUB, GRID, TRANS, LREAL, LSUB)
      USE prec
      USE mpimy

      IMPLICIT NONE

      LOGICAL LREAL,LSUB
      TYPE (grid_3d)     GRID
      TYPE (grid_3d)     GRID_SUB
      TYPE (transit)     TRANS

      GRID_SUB%NGX_rd=xm(GRID_SUB%NGX)
      GRID_SUB%NGY_rd=GRID_SUB%NGY
      GRID_SUB%NGZ_rd=zm(GRID_SUB%NGZ)

#ifdef realmode
      GRID_SUB%LREAL =.TRUE.
#else
      GRID_SUB%LREAL =.FALSE.
#endif      

      CALL GEN_SUB_GRID_(GRID_SUB, GRID, TRANS, LREAL, LSUB)
      END SUBROUTINE

!***********************************************************************
!
!  generate sub grid
!
!***********************************************************************

      SUBROUTINE GEN_SUB_GRID(GRID_SUB, GRID, TRANS, LREAL, LSUB)
      USE prec
      USE mpimy

      IMPLICIT NONE

      LOGICAL LREAL,LSUB
      TYPE (grid_3d)     GRID
      TYPE (grid_3d)     GRID_SUB
      TYPE (transit)     TRANS

      GRID_SUB%NGX_rd=GRID_SUB%NGX
      GRID_SUB%NGY_rd=GRID_SUB%NGY
      GRID_SUB%NGZ_rd=GRID_SUB%NGZ
      
      GRID_SUB%LREAL=.FALSE.

      CALL GEN_SUB_GRID_(GRID_SUB, GRID, TRANS, LREAL, LSUB)
      END SUBROUTINE

!***********************************************************************
!
!  generate the sub grid
!
!***********************************************************************

      SUBROUTINE GEN_SUB_GRID_(GRID_SUB,GRID, TRANS,LREAL,LSUB)
      USE prec
      USE mpimy

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      LOGICAL LREAL,LSUB
      TYPE (grid_3d)     GRID
      TYPE (grid_3d)     GRID_SUB
      TYPE (transit)     TRANS

      GRID_SUB%RC%NALLOC=0
      GRID_SUB%IN%NALLOC=0
      GRID_SUB%RL%NALLOC=0

#ifdef MPI
!-----------------------------------------------------------------------
! parallel version sub grid
!-----------------------------------------------------------------------
    NODE_ME=GRID_SUB%COMM%NODE_ME
    IF (LSUB) THEN
      GRID_SUB%RC%NFAST= 1
      GRID_SUB%RC%NROW = GRID_SUB%NGX_rd
      ALLOCATE(GRID_SUB%RC%I2( GRID%RC%NCOL )) ! always save
      ALLOCATE(GRID_SUB%RC%I3( GRID%RC%NCOL ))

      NCOL=GRID%RC%NCOL

      IND=0
      DO NC=1,NCOL
        N2=GRID%RC%I2(NC) ; L2=GRID%LPCTY(N2)
        N3=GRID%RC%I3(NC) ; L3=GRID%LPCTZ(N3)
!     is that one in the sub grid
        CALL SRCH_COL2(L2,L3,GRID_SUB,N2P,N3P)
        IF (N2P/=0 .AND. N3P/=0 ) THEN
          IND=IND+1
          GRID_SUB%RC%I2(IND)=N2P
          GRID_SUB%RC%I3(IND)=N3P
        ENDIF
      ENDDO
      GRID_SUB%RC%NCOL = IND
      IF (GRID_SUB%RC%NCOL >GRID%RC%NCOL) THEN
        WRITE(*,*) 'GEN_GRID_SUB: internal error',GRID_SUB%RC%NCOL
        STOP
      ENDIF
      GRID_SUB%RC%NP    = GRID_SUB%RC%NCOL*GRID_SUB%RC%NROW
      GRID_SUB%RC%NALLOC= GRID_SUB%RC%NCOL*GRID_SUB%RC%NROW
    ELSE
      CALL REC_STDLAY(GRID_SUB)
    ENDIF

      IF (LREAL) THEN
         CALL INTER_STDLAY(GRID_SUB)
         CALL REAL_STDLAY (GRID_SUB)
      ENDIF

#else
!-----------------------------------------------------------------------
! no problems with conventional version
! all grid points are contained within 3d box
!-----------------------------------------------------------------------
      GRID_SUB%RC%NFAST= 1
      GRID_SUB%RC%NCOL = GRID_SUB%NGZ_rd*GRID_SUB%NGY
      GRID_SUB%RC%NROW = GRID_SUB%NGX_rd
      ALLOCATE(GRID_SUB%RC%I2( GRID_SUB%RC%NCOL ))
      ALLOCATE(GRID_SUB%RC%I3( GRID_SUB%RC%NCOL ))
      IND=1
      DO N3=1,GRID_SUB%NGZ_rd
      DO N2=1,GRID_SUB%NGY
        GRID_SUB%RC%I2(IND)=N2
        GRID_SUB%RC%I3(IND)=N3
        IND=IND+1
      ENDDO
      ENDDO
      GRID_SUB%RC%NP    = GRID_SUB%RC%NCOL*GRID_SUB%RC%NROW
      GRID_SUB%RC%NALLOC= GRID_SUB%RC%NCOL*GRID_SUB%RC%NROW

      IF (LREAL) THEN
        CALL REAL_STDLAY (GRID_SUB)
      ENDIF
#endif

!=======================================================================
!
! create transition table for fast index (always x) quite simple
    IF (LSUB) THEN
      NGX=GRID_SUB%NGX

      ALLOCATE (TRANS%IND1(NGX))

      DO N1=1,(NGX/2)+1
        TRANS%IND1(N1)=N1
      ENDDO
      DO N1=(NGX/2)+2,NGX
        TRANS%IND1(N1)=N1-NGX+GRID%NGX
      ENDDO
!
! create transition table for columns
      NCOL=GRID_SUB%RC%NCOL
      ALLOCATE (TRANS%INDCOL(NCOL))

      DO NC=1,NCOL
        N2=GRID_SUB%RC%I2(NC) ; L2=GRID_SUB%LPCTY(N2)
        N3=GRID_SUB%RC%I3(NC) ; L3=GRID_SUB%LPCTZ(N3)
        CALL SRCH_COL1(L2,L3,GRID,NCP)
        TRANS%INDCOL(NC)=NCP
      ENDDO
    ENDIF
      GRID_SUB%MPLWV=MAX(GRID_SUB%RC%NALLOC ,GRID_SUB%IN%NALLOC , GRID_SUB%RL%NALLOC)
      DWRITE 'grid set up return',NODE_ME,GRID_SUB%NGX,GRID_SUB%NGY,GRID_SUB%NGZ,GRID_SUB%MPLWV


      END SUBROUTINE
!
! helper subroutine to search a specific column in a sub grid
! must succeed
!
      SUBROUTINE SRCH_COL1(L2S,L3S,GRID,NC2)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID

      NCOL=GRID%RC%NCOL
      DO NC=1,NCOL
        N2=GRID%RC%I2(NC) ; L2=GRID%LPCTY(N2)
        N3=GRID%RC%I3(NC) ; L3=GRID%LPCTZ(N3)
        IF (L2==L2S .AND. L3==L3S) THEN
          NC2=NC
          RETURN
        ENDIF
      ENDDO
      WRITE(*,*)'SRCH_COL1: internal error',L2S,L3S
      WRITE(*,'(20I6)') GRID%LPCTY
      WRITE(*,*)
      WRITE(*,'(20I6)') GRID%LPCTZ
      STOP
      END SUBROUTINE

!
! helper subroutine to check whether one specific column is in a
! sub grid: on sucess N2P and N3P is set
!
      SUBROUTINE SRCH_COL2(L2,L3,GRID,N2,N3)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID

      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY
        IF (L2==GRID%LPCTY(N2) .AND. L3==GRID%LPCTZ(N3)) RETURN
      ENDDO
      ENDDO
      N2=0
      N3=0

      END SUBROUTINE

!*************************SUBROUTINE REC_STDLAY   **********************
!
! setup reciprocal layout for FFT only parallel version
! columns are distributed among procressors in a round robin fashion
!
!***********************************************************************

      SUBROUTINE REC_STDLAY(GRID)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID
#ifdef MPI
      NCOL_TOT= GRID%NGZ_rd*GRID%NGY
      NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
      ALLOCATE(GRID%RC%I2( NCOL_MAX ))
      ALLOCATE(GRID%RC%I3( NCOL_MAX ))

      GRID%RC%NFAST= 1
      GRID%RC%NCOL = 0
      GRID%RC%NROW = GRID%NGX_rd
! distribute columns among proc
      IND=0
      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY
        NODE_TARGET =MOD(IND,GRID%COMM%NCPU)+1
        IND=IND+1
        IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
            N=GRID%RC%NCOL
            N=N+1
            GRID%RC%NCOL=N
            GRID%RC%I2(N)=N2 ! I2 contains y index
            GRID%RC%I3(N)=N3 ! I3      the z index
       ENDIF
      ENDDO; ENDDO
      IF (GRID%RC%NCOL >NCOL_MAX) THEN
        WRITE(*,*) 'REC_STDLAY: internal error',NCOL_MAX,GRID%RC%NCOL
        STOP
      ENDIF
      GRID%RC%NP    = GRID%RC%NCOL*GRID%RC%NROW
      GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
#else
      WRITE(*,*)'internal ERROR: REC_STDLAY'
      STOP
#endif
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE INTER_STDLAY **********************
!
! setup intermediate layout for FFT only parallel version
! columns are distributed among procressors in a round robin fashion
!
!***********************************************************************

      SUBROUTINE INTER_STDLAY(GRID)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID
#ifdef MPI

      NCOL_TOT= GRID%NGX_rd*GRID%NGZ_rd
      NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
      IF (LPLANE_WISE) THEN
         NCOL_MAX=GRID%NGZ_rd*((GRID%NGX_rd+GRID%COMM%NCPU-1)/GRID%COMM%NCPU)
      ENDIF


      ALLOCATE(GRID%IN%I2( NCOL_MAX ))
      ALLOCATE(GRID%IN%I3( NCOL_MAX ))

      GRID%IN%NFAST= 2
      GRID%IN%NCOL = 0
      GRID%IN%NROW = GRID%NGY
! distribute planes among proc
      IND=0
      DO N3=1,GRID%NGX_rd
      DO N2=1,GRID%NGZ_rd
        IF( LPLANE_WISE ) THEN
          NODE_TARGET=MOD(N3-1,GRID%COMM%NCPU)+1
        ELSE
          NODE_TARGET=MOD(IND,GRID%COMM%NCPU)+1
        ENDIF

        IND=IND+1
        IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
          N=GRID%IN%NCOL
          N=N+1
          GRID%IN%NCOL=N
          GRID%IN%I2(N)=N2 ! I2 contains z index
          GRID%IN%I3(N)=N3 ! I3      the x index
        ENDIF
      ENDDO
      ENDDO
      IF (GRID%IN%NCOL >NCOL_MAX) THEN
        WRITE(*,*) 'INTER_STDLAY: internal error',NCOL_MAX,GRID%IN%NCOL
        STOP
      ENDIF
      GRID%IN%NP    = GRID%IN%NCOL*GRID%IN%NROW
      GRID%IN%NALLOC= GRID%IN%NCOL*GRID%IN%NROW
#else
      WRITE(*,*)'internal ERROR: INTER_STDLAY'
      STOP
#endif
      END SUBROUTINE

!*************************SUBROUTINE REAL_STDLAY ***********************
!
! setup standard layout for real space
! this is quite simple
! all points are distributed among procressors in a round robin fashion
!
!***********************************************************************


      SUBROUTINE REAL_STDLAY(GRID)
      USE prec
      USE mpimy

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID
#ifdef MPI
!-----------------------------------------------------------------------
! version for parallel computers
! (currently zy planes are on one proc, z is fast index)
!-----------------------------------------------------------------------
      NCOL_TOT= GRID%NGY*GRID%NGX
      NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
      IF (LPLANE_WISE) THEN
      NCOL_MAX=GRID%NGY*((GRID%NGX+GRID%COMM%NCPU-1)/GRID%COMM%NCPU)
      ENDIF

      ALLOCATE(GRID%RL%I2( NCOL_MAX ))
      ALLOCATE(GRID%RL%I3( NCOL_MAX ))
      ALLOCATE(GRID%RL%INDEX(0:GRID%NGX-1,0:GRID%NGY-1))

      GRID%RL%NFAST= 3
      GRID%RL%NCOL = 0
      GRID%RL%NROW = GRID%NGZ
      GRID%RL%INDEX= 0
! distribute planes among proc
      IND=0
      DO N3=1,GRID%NGY
      DO N2=1,GRID%NGX
        IF( LPLANE_WISE ) THEN
          NODE_TARGET=MOD(N2-1,GRID%COMM%NCPU)+1
        ELSE
          NODE_TARGET=MOD(IND,GRID%COMM%NCPU)+1
        ENDIF

        IND=IND+1
        IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
          N=GRID%RL%NCOL
          N=N+1
          GRID%RL%INDEX(N2-1,N3-1)=N
          GRID%RL%NCOL=N
          GRID%RL%I2(N)=N2 ! I2 contains x index
          GRID%RL%I3(N)=N3 ! I3      the y index
        ENDIF
      ENDDO
      ENDDO
      IF (GRID%RL%NCOL >NCOL_MAX) THEN
        WRITE(*,*) 'INTER_STDLAY: internal error',NCOL_MAX,GRID%RL%NCOL
        STOP
      ENDIF
      GRID%RL%NP    = GRID%RL%NCOL*GRID%RL%NROW
      IF (GRID%LREAL) THEN
         GRID%RL%NALLOC=NCOL_MAX*GRID%NGZ_rd
      ELSE
         GRID%RL%NALLOC=NCOL_MAX*GRID%NGZ
      ENDIF
#else
!-----------------------------------------------------------------------
! conventional version (really simple)
! include all points x fast
!-----------------------------------------------------------------------
      GRID%RL%NFAST= 1
      GRID%RL%NCOL = GRID%NGZ*GRID%NGY
      GRID%RL%NROW = GRID%NGX

      ALLOCATE(GRID%RL%INDEX(0:GRID%NGY-1,0:GRID%NGZ-1))
      ALLOCATE(GRID%RL%I2( GRID%RL%NCOL ))
      ALLOCATE(GRID%RL%I3( GRID%RL%NCOL ))

      GRID%RL%INDEX=0

      IND=0
      DO N3=1,GRID%NGZ
      DO N2=1,GRID%NGY
        IND=IND+1
        GRID%RL%INDEX(N2-1,N3-1)=IND
        GRID%RL%I2(IND)=N2
        GRID%RL%I3(IND)=N3
      ENDDO
      ENDDO

      GRID%RL%NP    = GRID%RL%NCOL*GRID%RL%NROW
      IF (GRID%LREAL) THEN
         GRID%RL%NALLOC= GRID%RL%NCOL*GRID%RL%NROW/2
      ELSE
         GRID%RL%NALLOC= GRID%RL%NCOL*GRID%RL%NROW
      ENDIF
#endif
      END SUBROUTINE





!************************************************************************
!
! copy real space layout from one GRID_SRC to another GRID_DEST
! it might be that the GRID_DEST    contains more nodes in this
! case remaining nodes in GRID_DEST will have no columns in real space
! for instance consider the following topology (see M_divide in mpi.F)
! wave1           0 (0,0)        1 (0,1)
! wave2           2 (1,0)        3 (1,1)
! wave3           4 (2,0)        5 (2,1)
! wave4           6 (3,0)        7 (3,1)
! GRID_DEST%COMM contains all nodes 0-7, whereas GRID_SRC%COMM is an
! cartesian sub-communicator (in-band-group)
! (which communicate along first row, i.e proc 0 and 1).
! in real space node 0 and 1 will store same columns as GRID_SRC
! all other nodes will have no data in real space
!
!***********************************************************************


      SUBROUTINE  SET_RL_GRID(GRID_DEST,GRID_SRC)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (grid_3d)     GRID_DEST,GRID_SRC
#ifdef MPI
!-----------------------------------------------------------------------
! real space
! see description in mpi.F to see how this routine works
!-----------------------------------------------------------------------
! deallocate old RL space data prior to copy
      DEALLOCATE(GRID_DEST%RL%I2,GRID_DEST%RL%I3,GRID_DEST%RL%INDEX)
!
      IF (GRID_DEST%NGX /= GRID_SRC%NGX .OR. GRID_DEST%NGY /= GRID_SRC%NGY .OR. GRID_DEST%NGZ /= GRID_SRC%NGZ) THEN
         WRITE(0,*) 'ERROR in SET_RL_GRID: GRID_DEST is not conform to GRID_SRC'
         STOP
      ENDIF
      IF (GRID_DEST%COMM%NODE_ME <= GRID_SRC%COMM%NCPU) THEN
        GRID_DEST%RL%NFAST= GRID_SRC%RL%NFAST
        GRID_DEST%RL%NCOL = GRID_SRC%RL%NCOL
        GRID_DEST%RL%NALLOC=GRID_DEST%RL%NCOL*GRID_DEST%NGZ_rd
        GRID_DEST%RL%NROW = GRID_SRC%RL%NROW
        GRID_DEST%RL%NP    =GRID_SRC%RL%NP
        GRID_DEST%RL%INDEX=>GRID_SRC%RL%INDEX
        GRID_DEST%RL%I2   =>GRID_SRC%RL%I2
        GRID_DEST%RL%I3   =>GRID_SRC%RL%I3
      ELSE
        GRID_DEST%RL%NFAST =3
        GRID_DEST%RL%NCOL  =0
        GRID_DEST%RL%NROW  =0
        GRID_DEST%RL%NALLOC=0
        GRID_DEST%RL%NP    =0
        NULLIFY(GRID_DEST%RL%I2); NULLIFY(GRID_DEST%RL%I3)
        NULLIFY(GRID_DEST%RL%INDEX)
      ENDIF
      GRID_DEST%MPLWV=MAX(GRID_DEST%RC%NALLOC ,GRID_DEST%IN%NALLOC , GRID_DEST%RL%NALLOC)
#else
! nothing
#endif
      END SUBROUTINE
      END MODULE


!*************************SUBROUTINE RC_FLIP ***************************
!
! rearranges the storage mode for spin components of charge density:
! 
! given rho_up and rho_down on input the quantities (rho_up+rho_down)
! and (rho_up-rho_down) = total charge and magnetization are returned
! and also the reverse operation is possible if setting LBACK=.TRUE.
! RC_FLIP  operates in reciprocal space (half gird)
! RL_FLIP  operates in real space       (possibly real arrays)
!
! for the non collinear version ISPIN must be set to 4 by the caller
!
! given the total 2x2 density matrix (CHTOT = p00, p11, p01, p10 ) 
! the  charge density (C00) and the magnetisation density 
! (CX, CY, CY) are calculated 
! If LBACK is defined, and given the charge density (C00) and the 
! magnetisation density (CX, CY, CY) the total 2x2 density matrix 
! (CHTOT = p00, p11, p01, p10 ) is calculated
!
!***********************************************************************

!
!  version for operating in reciprocal space (half grid version)
!
      SUBROUTINE RC_FLIP(CHTOT, GRID,ISPIN,LBACK)
      USE prec
      USE mpimy
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      LOGICAL LBACK
      COMPLEX(q) CHTOT(GRID%MPLWV,ISPIN)

      IF (ISPIN==2) THEN

         FAC=1._q
         IF (LBACK) FAC=0.5_q
         DO K=1,GRID%RC%NP
            CQU=CHTOT(K,1)
            CQD=CHTOT(K,2)
            CHTOT(K,1)=FAC*(CQU+CQD)
            CHTOT(K,2)=FAC*(CQU-CQD)
         ENDDO
      ELSE IF ( ISPIN==4 .AND. .NOT. LBACK) THEN
         DO K=1,GRID%RC%NP
            C00=CHTOT(K,1)
            C01=CHTOT(K,2)
            C10=CHTOT(K,3)
            C11=CHTOT(K,4)

            CHTOT(K,1)= C00+C11             
            CHTOT(K,2)= C01+C10             
            CHTOT(K,3)=(C01-C10)*(0._q,1._q)
            CHTOT(K,4)= C00-C11             
         ENDDO
      ELSE IF ( ISPIN==4 .AND. LBACK) THEN
         FAC=0.5_q
         DO K=1,GRID%RC%NP
            C00=CHTOT(K,1)
            CX =CHTOT(K,2)
            CY =CHTOT(K,3)
            CZ =CHTOT(K,4)

            CHTOT(K,1)= (C00+CZ)*FAC           
            CHTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CHTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CHTOT(K,4)= (C00-CZ)*FAC           
         ENDDO
      ENDIF

      END SUBROUTINE

!
!  version for operating in real space (possibly real arrays)
!
      SUBROUTINE RL_FLIP(CHTOTR, GRID,ISPIN,LBACK)
      USE prec
      USE mpimy
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      LOGICAL LBACK

      RGRID CHTOTR(DIMREAL(GRID%MPLWV),ISPIN)

      IF (ISPIN==2) THEN
         FAC=1._q
         IF (LBACK) FAC=0.5_q
         DO K=1,GRID%RL%NP
            CQU=CHTOTR(K,1)
            CQD=CHTOTR(K,2)
            CHTOTR(K,1)=FAC*(CQU+CQD)
            CHTOTR(K,2)=FAC*(CQU-CQD)
         ENDDO
      ELSE IF ( ISPIN==4 .AND. .NOT. LBACK) THEN
         DO K=1,GRID%RL%NP
            C00=CHTOTR(K,1)
            C01=CHTOTR(K,2)
            C10=CHTOTR(K,3)
            C11=CHTOTR(K,4)

            CHTOTR(K,1)= C00+C11             
            CHTOTR(K,2)= C01+C10             
            CHTOTR(K,3)=(C01-C10)*(0._q,1._q)
            CHTOTR(K,4)= C00-C11             
         ENDDO
      ELSE IF ( ISPIN==4 .AND. LBACK) THEN
         FAC=0.5_q
         DO K=1,GRID%RL%NP
            C00=CHTOTR(K,1)
            CX =CHTOTR(K,2)
            CY =CHTOTR(K,3)
            CZ =CHTOTR(K,4)

            CHTOTR(K,1)= (C00+CZ)*FAC           
            CHTOTR(K,2)= (CX-CY*(0._q,1._q))*FAC
            CHTOTR(K,3)= (CX+CY*(0._q,1._q))*FAC
            CHTOTR(K,4)= (C00-CZ)*FAC           
         ENDDO
      ENDIF
      END SUBROUTINE


!************************ SUBROUTINE MRG_GRID_RL ***********************
!
! helper routine to merge a grid-array from all nodes to the
! local node (real space version)
! result is always real and has standard serial layout
!
!***********************************************************************

      SUBROUTINE MRG_GRID_RL(GRID, C, C_LOCAL)
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
#ifdef MPI
      RGRID   C_LOCAL(GRID%RL%NROW,GRID%RL%NCOL)
      REAL(q) C(GRID%NGX,GRID%NGY,GRID%NGZ)
      INTEGER NC,NX,NY,NZ

      C=0
      DO NC=1,GRID%RL%NCOL
        NX= GRID%RL%I2(NC)
        NY= GRID%RL%I3(NC)
        DO NZ=1,GRID%NGZ
          C(NX,NY,NZ)=C_LOCAL(NZ,NC)
        ENDDO
      ENDDO

      CALLMPI( M_sum_d(GRID%COMM, C(1,1,1), GRID%NPLWV))
#else
      RGRID    C_LOCAL(GRID%RL%NP)
      REAL(q)  C(GRID%RL%NP)

      C(1:GRID%RL%NP)=C_LOCAL(1:GRID%RL%NP)
#endif

      END SUBROUTINE

!************************ SUBROUTINE MRG_GRID_RL_PLANE *****************
!
! helper routine to merge a plane of a grid-array from all nodes to the
! local node (real space version)
! result is always real and has standard serial layout
!
!***********************************************************************

      SUBROUTINE MRG_GRID_RL_PLANE(GRID, C, C_LOCAL, NZ)
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
#ifdef MPI
      RGRID   C_LOCAL(GRID%RL%NROW,GRID%RL%NCOL)
      REAL(q) C(GRID%NGX,GRID%NGY)
      INTEGER NC,NX,NY,NZ

      C=0
      DO NC=1,GRID%RL%NCOL
        NX= GRID%RL%I2(NC)
        NY= GRID%RL%I3(NC)
          C(NX,NY)=C_LOCAL(NZ,NC)
      ENDDO

      CALLMPI( M_sum_d(GRID%COMM, C(1,1), GRID%NGX*GRID%NGY))
#else
      INTEGER NZ
      RGRID    C_LOCAL(GRID%NGX*GRID%NGY,NZ)
      REAL(q)  C(GRID%NGX*GRID%NGY)

      C=C_LOCAL(:,NZ)
#endif

      END SUBROUTINE

!************************ SUBROUTINE DIS_GRID_RL ***********************
!
! helper routine to distribute an array from io-node to the
! local nodes (real space version)
! input must be always real and must have  standard serial layout
!
!***********************************************************************

      SUBROUTINE DIS_GRID_RL(GRID, C, C_LOCAL, LBROADCAST )
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      LOGICAL LBROADCAST
#ifdef MPI
      RGRID  C_LOCAL(GRID%RL%NROW,GRID%RL%NCOL)
      REAL(q)   C(GRID%NGX,GRID%NGY,GRID%NGZ)
      INTEGER NC,NX,NY,NZ

      ! broadcast C to all nodes from IONODE
      IF (LBROADCAST) &
        CALLMPI( M_bcast_d(GRID%COMM, C(1,1,1), GRID%NPLWV))
      ! and pick up data
      DO NC=1,GRID%RL%NCOL
        NX= GRID%RL%I2(NC)
        NY= GRID%RL%I3(NC)
        DO NZ=1,GRID%NGZ
          C_LOCAL(NZ,NC)=C(NX,NY,NZ)
        ENDDO
      ENDDO
#else
      RGRID C_LOCAL(GRID%RL%NP)
      REAL(q)  C(GRID%RL%NP)

      C_LOCAL(1:GRID%RL%NP)=C(1:GRID%RL%NP)
#endif

      END SUBROUTINE

!************************ SUBROUTINE DIS_GRID_RL ***********************
!
! helper routine to distribute a plane of a grid based array from io-node 
! to the local nodes (real space version)
! input must be always real and must have  standard serial layout
!
!***********************************************************************

      SUBROUTINE DIS_GRID_RL_PLANE(GRID, C, C_LOCAL, LBROADCAST , NZ  )
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      LOGICAL LBROADCAST
#ifdef MPI
      RGRID  C_LOCAL(GRID%RL%NROW,GRID%RL%NCOL)
      REAL(q)   C(GRID%NGX,GRID%NGY)
      INTEGER NC,NX,NY,NZ

      ! broadcast C to all nodes from IONODE
      IF (LBROADCAST) &
        CALLMPI( M_bcast_d(GRID%COMM, C(1,1), GRID%NGX*GRID%NGY))
      ! and pick up data
      DO NC=1,GRID%RL%NCOL
        NX= GRID%RL%I2(NC)
        NY= GRID%RL%I3(NC)
        C_LOCAL(NZ,NC)=C(NX,NY)
      ENDDO
#else
      INTEGER NZ
      RGRID    C_LOCAL(GRID%NGX*GRID%NGY,NZ)
      REAL(q)  C(GRID%NGX*GRID%NGY)

      C_LOCAL(:,NZ)=C
#endif

      END SUBROUTINE

!************************ SUBROUTINE MRG_GRID_RC ***********************
!
! helper routine to merge a grid-array from all nodes to the
! local node (reciprocal space version)
! result is always real and has standard serial layout
!
!***********************************************************************

      SUBROUTINE MRG_GRID_RC(GRID, C, C_LOCAL)
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
#ifdef MPI
      COMPLEX(q) C_LOCAL(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) C(GRID%NGX_rd,GRID%NGY,GRID%NGZ_rd)
      INTEGER NC,NX,NY,NZ,NPLWV

      C=0
      DO NC=1,GRID%RC%NCOL
        NY= GRID%RC%I2(NC)
        NZ= GRID%RC%I3(NC)
        DO NX=1,GRID%NGX_rd
          C(NX,NY,NZ)=C_LOCAL(NX,NC)
        ENDDO
      ENDDO

      NPLWV=GRID%NGX_rd*GRID%NGY*GRID%NGZ_rd
      CALLMPI( M_sum_z(GRID%COMM, C(1,1,1), NPLWV))
#else
      COMPLEX(q) C_LOCAL(GRID%RC%NP)
      COMPLEX(q) C(GRID%RC%NP)

      C(1:GRID%RC%NP)=C_LOCAL(1:GRID%RC%NP)
#endif

      END SUBROUTINE

!************************ SUBROUTINE DIS_GRID_RC ***********************
!
! helper routine to distribute an array from io-node to the
! local nodes (reciprocal space version)
! input must be always real and must have  standard serial layout
!
!***********************************************************************

      SUBROUTINE DIS_GRID_RC(GRID, C, C_LOCAL, LBROADCAST )
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      LOGICAL LBROADCAST
#ifdef MPI
      COMPLEX(q) C_LOCAL(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) C(GRID%NGX_rd,GRID%NGY,GRID%NGZ_rd)
      INTEGER NC,NX,NY,NZ,NPLWV
      NPLWV=GRID%NGX_rd*GRID%NGY*GRID%NGZ_rd

      ! broadcast C to all nodes from IONODE
      IF (LBROADCAST) &
        CALLMPI( M_bcast_z(GRID%COMM, C(1,1,1), NPLWV))
      ! and pick up data
      DO NC=1,GRID%RC%NCOL
        NY= GRID%RC%I2(NC)
        NZ= GRID%RC%I3(NC)
        DO NX=1,GRID%NGX_rd
          C_LOCAL(NX,NC)=C(NX,NY,NZ)
        ENDDO
      ENDDO
#else
      COMPLEX(q) C_LOCAL(GRID%RC%NP)
      COMPLEX(q) C(GRID%RC%NP)

      C_LOCAL(1:GRID%RC%NP)=C(1:GRID%RC%NP)
#endif

      END SUBROUTINE



!************************ SUBROUTINE RL_ADD  ***************************
!
!  This modul contains a helper routine to perform the operation
!     C = A * SCALE1 + B * SCALE2
!  in REAL(q) SPACE
!  depending on the mode this is done using real or complex
!  arrays
!  C might be equal to A
!
! !! there is no need to optimized these routines                    !!
!
!***********************************************************************

      SUBROUTINE RL_ADD(A,SCALE1,B,SCALE2,C,GRID)
      USE prec
      USE mpimy
      USE mgrid

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      RGRID  C(GRID%RL%NP),A(GRID%RL%NP),B(GRID%RL%NP)

      NP=GRID%RL%NP
!
!   C=A*SCALE1
!
!  )  case SCALE1=0
      IF (SCALE1==0)  THEN
        DO 100 N=1,NP
          C(N)=0
  100   CONTINUE
!  )  case SCALE1=1
      ELSE IF (SCALE1==1)  THEN
        DO 200 N=1,NP
          C(N)=A(N)
  200   CONTINUE
!  )  else
      ELSE
        DO 300 N=1,NP
          C(N)=A(N)*SCALE1
  300   CONTINUE
      ENDIF

!
!   C=C+ B*SCALE2
!
      IF (SCALE2==0) THEN
      ELSE IF (SCALE2==1) THEN
        DO 400 N=1,NP
          C(N)=C(N)+B(N)
  400   CONTINUE
      ELSE
        DO 500 N=1,NP
          C(N)=C(N)+B(N)*SCALE2
  500   CONTINUE
      ENDIF

      RETURN
      END

!************************ SUBROUTINE RC_ADD  ***************************
!
!  This modul contains a helper routine to perform the operation
!     C = A * SCALE1 + B * SCALE2
!  in RECIPROCAL SPACE
!  C might be equal to A
!
!***********************************************************************

      SUBROUTINE RC_ADD(A,SCALE1,B,SCALE2,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) C(GRID%RC%NP),A(GRID%RC%NP),B(GRID%RC%NP)
      NP=GRID%RC%NP
!
!   C=A*SCALE1
!
!  )  case SCALE1=0
      IF (SCALE1==0)  THEN
        DO 100 N=1,NP
          C(N)=0
  100   CONTINUE
!  )  case SCALE1=1
      ELSE IF (SCALE1==1)  THEN
        DO 200 N=1,NP
          C(N)=A(N)
  200   CONTINUE
!  )  else
      ELSE
        DO 300 N=1,NP
          C(N)=A(N)*SCALE1
  300   CONTINUE
      ENDIF

!
!   C=C+ B*SCALE2
!
      IF (SCALE2==0) THEN
      ELSE IF (SCALE2==1) THEN
        DO 400 N=1,NP
          C(N)=C(N)+B(N)
  400   CONTINUE
      ELSE
        DO 500 N=1,NP
          C(N)=C(N)+B(N)*SCALE2
  500   CONTINUE
      ENDIF

      RETURN
      END

!************************ SUBROUTINE RLR_MUL  **************************
!  This modul contains a helper routine to perform the operation
!     C = A * B * SCALE
!  in REAL(q) SPACE, where B is an real array
!                               ----------
!  depending on the mode this is done using real or complex
!  arrays C and A
!   C might be equal to A
!
! !! there is no need to optimized these routines                    !!
!***********************************************************************

      SUBROUTINE RLR_MUL(A,B,SCALE,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      REAL(q)  B(GRID%RL%NP)
      RGRID C(GRID%RL%NP),A(GRID%RL%NP)

      NPLWV=GRID%RL%NP
!
!   C=0
!
!  )  case SCALE=0
      IF (SCALE==0)  THEN
        DO N=1,NPLWV
          C(N)=0
        ENDDO
!  )  case SCALE=1
      ELSE IF (SCALE==1)  THEN
        DO N=1,NPLWV
          C(N)=A(N)*B(N)
        ENDDO
!  )  else
      ELSE
        DO N=1,NPLWV
          C(N)=A(N)*B(N)*SCALE
        ENDDO
      ENDIF

      RETURN
      END





!***********************************************************************
!
!  This subroutine adds an array defined on the small grid
!  GRID_SOFT to  an array defined on the fine grid GRID
!  both arrays are in the half grid mode
!
!***********************************************************************

      SUBROUTINE ADD_GRID(GRID,GRID_SOFT,SOFT_TO_C,CADD,CSUM)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      COMPLEX(q) CSUM(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) CADD(GRID_SOFT%RC%NROW,GRID_SOFT%RC%NCOL)

      col: DO NC=1,GRID_SOFT%RC%NCOL
      row: DO N1=1,GRID_SOFT%RC%NROW
        CSUM(SOFT_TO_C%IND1(N1),SOFT_TO_C%INDCOL(NC))= &
        CSUM(SOFT_TO_C%IND1(N1),SOFT_TO_C%INDCOL(NC))+ CADD(N1,NC)
      ENDDO row
      ENDDO col
      RETURN
      END


!***********************************************************************
!
!  This subroutine copies an array defined on the small grid
!  GRID_SOFT to  an array defined on the fine grid GRID
!  both arrays are in the half grid mode
!
!***********************************************************************

      SUBROUTINE CPB_GRID(GRID,GRID_SOFT,SOFT_TO_C,CSRC,CDEST)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      COMPLEX(q) CDEST(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) CSRC(GRID_SOFT%RC%NROW,GRID_SOFT%RC%NCOL)

      col: DO NC=1,GRID_SOFT%RC%NCOL
      row: DO N1=1,GRID_SOFT%RC%NROW
        CDEST(SOFT_TO_C%IND1(N1),SOFT_TO_C%INDCOL(NC))=CSRC(N1,NC)
      ENDDO row
      ENDDO col
      RETURN
      END

!***********************************************************************
!
!  This subroutine brings an array defined on the fine grid
!  to an array defined on the wavefunctions grid
!  both arays are stored in the half grid mode
!
!***********************************************************************

      SUBROUTINE CP_GRID(GRID,GRID_SOFT,SOFT_TO_C,CSRC,CDEST)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      COMPLEX(q) CSRC(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) CDEST(GRID_SOFT%RC%NROW,GRID_SOFT%RC%NCOL)

      col: DO NC=1,GRID_SOFT%RC%NCOL
      row: DO N1=1,GRID_SOFT%RC%NROW
        CDEST(N1,NC)=CSRC(SOFT_TO_C%IND1(N1),SOFT_TO_C%INDCOL(NC))
      ENDDO row
      ENDDO col
      RETURN
      END

!***********************SUBROUTINE SETUNB ******************************
!
! set contribution of unbalanced lattic-vectors in an array to zero
! if wrap-arround errors are ommited the contribution is zero
!  operates on the reduced half grid
!**********************************************************************

      SUBROUTINE SETUNB(C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID

      COMPLEX(q) C(GRID%RC%NROW,GRID%RC%NCOL)

      N2= (GRID%NGY/2)+1
      N3= (GRID%NGZ/2)+1

      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2) THEN
          DO N1=1,GRID%RC%NROW
            C(N1,NC)=0
          ENDDO
        ENDIF
        IF (GRID%RC%I3(NC)==N3) THEN
          DO N1=1,GRID%RC%NROW
            C(N1,NC)=0
          ENDDO
        ENDIF
      ENDDO

      N1= (GRID%NGX/2)+1
      DO NC=1,GRID%RC%NCOL
       C(N1,NC)=0
      ENDDO

      RETURN
      END SUBROUTINE

!***********************SUBROUTINE SETUNB ******************************
!
! set contribution of unbalanced lattic-vectors in an array to zero
! if wrap-arround errors are ommited the contribution is zero
!  operates on the reduced half grid
!**********************************************************************

      SUBROUTINE SETUNB_COMPAT(C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)     GRID

      COMPLEX(q) C(GRID%RC%NROW,GRID%RC%NCOL)

      IF (.NOT. LCOMPAT) RETURN

      N2= (GRID%NGY/2)+1
      N3= (GRID%NGZ/2)+1

      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2) THEN
          DO N1=1,GRID%RC%NROW
            C(N1,NC)=0
          ENDDO
        ENDIF
        IF (GRID%RC%I3(NC)==N3) THEN
          DO N1=1,GRID%RC%NROW
            C(N1,NC)=0
          ENDDO
        ENDIF
      ENDDO

      N1= (GRID%NGX/2)+1
      DO NC=1,GRID%RC%NCOL
       C(N1,NC)=0
      ENDDO

      RETURN
      END SUBROUTINE

!***********************************************************************
!
!  OPSYNC does nothing
!  just returns
!  this function is called (for instance) from xcgrad.F
!  to avoid that the compiler is too clever
!
!***********************************************************************

      SUBROUTINE OPSYNC(CSRC,CDEST,NPLWV)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)

      RGRID CSRC(NPLWV),CDEST(NPLWV)

      RETURN
      END

!***********************************************************************
!
!  sum all points on the grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUM_RC(STRING,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      COMPLEX(q) C(GRID%RC%NROW,GRID%RC%NCOL)
      CHARACTER*(*) STRING

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      CSUM=0
      DO NC=1,GRID%RC%NCOL
      DO N1=1,GRID%RC%NROW
       CSUM=CSUM+C(N1,NC)
      ENDDO
      ENDDO
      CALLMPI( M_sum_z(GRID%COMM, CSUM,1))
    io_begin
      WRITE(*,'("SUM_RC:",A,2E20.10)') STRING,CSUM
    io_end
      RETURN
      END

!***********************************************************************
!
!  sum all points on the grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUM_IN(STRING,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      COMPLEX(q) C(GRID%IN%NCOL,GRID%IN%NROW)
      CHARACTER*(*) STRING

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      CSUM=0
      DO NC=1,GRID%IN%NCOL
      DO N1=1,GRID%IN%NROW
       CSUM=CSUM+C(NC,N1)
      ENDDO
      ENDDO
      CALLMPI( M_sum_z(GRID%COMM, CSUM,1))
    io_begin
      WRITE(*,'("SUM_IN:",A,2E20.10)') STRING,CSUM
    io_end
      RETURN
      END

!***********************************************************************
!
!  sum all points on the real grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUM_RL(STRING,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      RGRID C(GRID%RL%NROW,GRID%RL%NCOL)
      CHARACTER*(*) STRING

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      CSUM=0
      DO NC=1,GRID%RL%NCOL
      DO N1=1,GRID%RL%NROW
       CSUM=CSUM+C(N1,NC)
      ENDDO
      ENDDO
      CALLMPI( M_sum_z(GRID%COMM, CSUM,1))
    io_begin
      WRITE(*,'("SUM_RL:",A,4E20.10)') STRING,CSUM/GRID%NPLWV,CSUM
    io_end
      RETURN
      END


!*************************SUBROUTINE INIDAT ****************************
!  subroutine to initialize a workarray to some non zero data
!***********************************************************************

      SUBROUTINE INIDAT(NPLWV,CWORK)
      USE prec
      COMPLEX(q)  CWORK(NPLWV)

      DO M=1,NPLWV
        CWORK(M)=(0.1_q,0.2_q)
      ENDDO
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE INIDATR ***************************
!  subroutine to initialize a workarray to some non zero data
!***********************************************************************

      SUBROUTINE INIDATR(NPLWV,CWORK)
      USE prec
      GDEF  CWORK(NPLWV)

      DO M=1,NPLWV
        CWORK(M)=(0.1_q,0.2_q)
      ENDDO
      RETURN
      END SUBROUTINE

!***********************SUBROUTINE WRT_RC_LINE ************************
!
!  dump the contents of an array (defined in reciprocal space)
!  along three lines to an specified unit
!  parallel version currently no support
!
!**********************************************************************

      SUBROUTINE WRT_RC_LINE(IU,GRID,CHDEN)
      USE prec
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (grid_3d) GRID
      COMPLEX(q) CHDEN(GRID%NGX_rd,GRID%NGY,GRID%NGZ_rd)
#ifndef MPI
      NGX=GRID%NGX
      NGY=GRID%NGY
      NGZ=GRID%NGZ

      WRITE(IU,'(2Hx ,10F10.4)') (REAL( CHDEN(NX,1,1) ,KIND=q) ,NX=1,(NGX/2)+1)
      WRITE(IU,'(2Hy ,10F10.4)') (REAL( CHDEN(1,NY,1) ,KIND=q) ,NY=1,(NGY/2)+1)
      WRITE(IU,'(2Hz ,10F10.4)') (REAL( CHDEN(1,1,NZ) ,KIND=q) ,NZ=1,(NGZ/2)+1)
!      WRITE(IU,'(2Hx ,10F10.4)') (CHDEN(NX,1,1) ,NX=1,(NGX/2)+1)
!      WRITE(IU,'(2Hy ,10F10.4)') (CHDEN(1,NY,1) ,NY=1,(NGY/2)+1)
!      WRITE(IU,'(2Hz ,10F10.4)') (CHDEN(1,1,NZ) ,NZ=1,(NGZ/2)+1)
#endif
      RETURN
      END SUBROUTINE
