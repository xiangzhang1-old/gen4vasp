!#define debug
#include "symbol.inc"
!************************************************************************
! RCS:  $Id: wave.F,v 1.6 2002/08/14 13:59:43 kresse Exp $
!
!  this module contains the routines required to setup
!  the distribution of wavefunctions over nodes and all basic routines
!  handling wavedes etc.
!
!***********************************************************************
      MODULE WAVE
      USE prec
      USE mpimy
      INCLUDE "wave.inc"
      CONTAINS

!=======================================================================
!  initialize and descriptor for the wavefunctions
!  mainly allocation
!=======================================================================

      SUBROUTINE ALLOCWDES(WDES,LEXTEND)
      USE prec
      IMPLICIT NONE

      INTEGER NK
      TYPE (wavedes)  WDES
      INTEGER NRPLWV,NKPTS,NCOL
      LOGICAL LEXTEND

      NRPLWV=WDES%NGDIM
      NKPTS =WDES%NKPTS
      NCOL  =WDES%NCOL      

      ALLOCATE( WDES%NPLWKP(NKPTS),WDES%NGVECTOR(NKPTS),WDES%NPLWKP_TOT(NKPTS),WDES%NINDPW(NRPLWV,NKPTS))
      IF (NCOL>0) THEN
        ALLOCATE(WDES%PL_INDEX(NCOL,NKPTS),WDES%PL_COL(NCOL,NKPTS))
      ELSE
        NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
      ENDIF

      IF (LEXTEND) THEN
!-MM- changes to accommodate spin spirals
! original statement
!       ALLOCATE( WDES%DATAKE(NRPLWV,NKPTS), &
!          WDES%IGX(NRPLWV,NKPTS),WDES%IGY(NRPLWV,NKPTS),WDES%IGZ(NRPLWV,NKPTS))
        ALLOCATE(WDES%DATAKE(NRPLWV,NKPTS,2), &
      &    WDES%IGX(NRPLWV,NKPTS),WDES%IGY(NRPLWV,NKPTS),WDES%IGZ(NRPLWV,NKPTS))
!-MM- end of alteration
      ELSE
        NULLIFY(WDES%DATAKE)
        NULLIFY(WDES%IGX); NULLIFY(WDES%IGY); NULLIFY(WDES%IGZ)
        NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
      END IF
      END SUBROUTINE

!=======================================================================
!  deallocate a  descriptor for the wavefunctions
!=======================================================================

      SUBROUTINE DEALLOCWDES(WDES,LEXTEND)
      USE prec
      IMPLICIT NONE
      TYPE (wavedes)  WDES
      LOGICAL LEXTEND

      DEALLOCATE( WDES%NPLWKP,WDES%NGVECTOR,WDES%NPLWKP_TOT,WDES%NINDPW)

      IF (WDES%NCOL>0) THEN
        DEALLOCATE(WDES%PL_INDEX,WDES%PL_COL)
        NULLIFY(WDES%PL_INDEX); NULLIFY(WDES%PL_COL)
      ENDIF
      IF (LEXTEND) THEN
        DEALLOCATE( WDES%DATAKE,WDES%IGX,WDES%IGY,WDES%IGZ)
        NULLIFY(WDES%DATAKE)
        NULLIFY(WDES%IGX); NULLIFY(WDES%IGY); NULLIFY(WDES%IGZ)
      END IF
      END SUBROUTINE


!=======================================================================
!  initialize the projector part of the descriptor for the
!  wavefunctions
!=======================================================================

      SUBROUTINE WDES_SET_NPRO(WDES,T_INFO,P)
      USE prec
      USE  mpimy
      USE  poscar
      USE  pseudo

      TYPE (wavedes)  WDES
      TYPE (type_info) :: T_INFO
      TYPE (potcar)   P(T_INFO%NTYP)
! local varibles
      INTEGER NALLOC,NPRO_TOT,NT,NI,NIS,NODE_TARGET,NPRO, &
             LMMAXC,NIONS,LASTTYP
#ifndef MPI
      WDES%NIONS = T_INFO%NIONS
      WDES%NTYP  = T_INFO%NTYP
      WDES%NITYP =>T_INFO%NITYP
      ALLOCATE(WDES%LMMAX(WDES%NTYP))
      DO NT=1,T_INFO%NTYP
        WDES%LMMAX(NT)=P(NT)%LMMAX
      ENDDO
      WDES%NPRO  =SUM(WDES%LMMAX*WDES%NITYP)
      WDES%NPRO_TOT=SUM(WDES%LMMAX*WDES%NITYP)
      WDES%NPROD =WDES%NPRO
#else
!-----------------------------------------------------------------------
! parallel version
!-----------------------------------------------------------------------
      TYPE (communic) COMM
! first count number of projection operators
      COMM = WDES%COMM_INB
      NPRO_TOT=0
      DO NT=1,T_INFO%NTYP
        NPRO_TOT=NPRO_TOT+P(NT)%LMMAX*T_INFO%NITYP(NT)
      ENDDO
! check implementation of NI_LOCAL and NI_GLOBAL
      DO NI=1,T_INFO%NIONS
        NI_L=NI_LOCAL(NI,COMM)
        IF (NI_L /= 0) THEN
          IF (NI /= NI_GLOBAL(NI_L,COMM)) THEN
            WRITE(*,*)' internal error NI_GLOBAL, NI_LOCAL'
            STOP
          ENDIF
        ENDIF
      ENDDO
! number of projection operators per node
      NALLOC=(T_INFO%NIONS+COMM%NCPU-1)/COMM%NCPU

      WDES%NIONS=0
      WDES%NTYP =0
      WDES%NPRO =0
      WDES%NPRO_TOT=NPRO_TOT
      ALLOCATE(WDES%NITYP(NALLOC))
      ALLOCATE(WDES%LMMAX(NALLOC))
      ALLOCATE(WDES%NPRO_POS(NALLOC))
      WDES%LMMAX=0
      WDES%NITYP=0

      LASTTYP=0
      NIS    =1

      NPRO=0
      DO NT=1,T_INFO%NTYP
        LMMAXC=P(NT)%LMMAX
        IF (LMMAXC==0) GOTO 100
        DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
          ! does this element reside on local node
          IF (NI_LOCAL(NI,COMM) /=0 ) THEN
            WDES%NIONS=WDES%NIONS+1
            WDES%NPRO_POS(WDES%NIONS)=NPRO
            WDES%NPRO =WDES%NPRO +LMMAXC
            IF (NT /= LASTTYP) THEN
              WDES%NTYP=WDES%NTYP+1
              LASTTYP=NT
            ENDIF
            WDES%LMMAX(WDES%NTYP)=LMMAXC
            WDES%NITYP(WDES%NTYP)=WDES%NITYP(WDES%NTYP)+1
          ENDIF
          NPRO= LMMAXC+NPRO
        ENDDO

 100    NIS = NIS+T_INFO%NITYP(NT)
      ENDDO
! check whether everything is right
      NPRO=SUM(WDES%LMMAX*WDES%NITYP)
      IF (NPRO/= WDES%NPRO) THEN
        WRITE(*,*)'internal ERROR 1 in WDES_SET_NPRO:',NPRO,WDES%NPRO
        STOP
      ENDIF
      CALLMPI( M_sum_i(COMM,NPRO ,1))
      IF (NPRO/= NPRO_TOT) THEN
        WRITE(*,*)'internal ERROR 2 in WDES_SET_NPRO:',NPRO,NPRO_TOT
        STOP
      ENDIF
! make NPROD dividable by number of NB_PAR
      WDES%NPROD=((WDES%NPRO+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
! and set it to the maximum value of all processors
      CALLMPI (M_max_i(COMM, WDES%NPROD, 1)) ! (required for  SHMALLOC)
#endif


      WDES%NPRO  =WDES%NPRO      *WDES%NRSPINORS
      WDES%NPRO_TOT=WDES%NPRO_TOT*WDES%NRSPINORS
      WDES%NPROD =WDES%NPROD*WDES%NRSPINORS
      END SUBROUTINE


!=======================================================================
!
! this routine gives the local storage index for
! the non local overlap CQIJ , strength CDIJ matrix elements
! and for the projected wavefunctions
! return is 0 if the element resides on an other processor
! on entry NI is the global index
!=======================================================================

      FUNCTION NI_LOCAL(NI,COMM)
      USE prec
      USE  mpimy
      IMPLICIT NONE
      TYPE (communic)  COMM
      INTEGER NI,NI_LOCAL,NODE_TARGET
#ifdef MPI

      NODE_TARGET=MOD(NI-1, COMM%NCPU)+1
      NI_LOCAL   =   (NI-1)/COMM%NCPU +1

      IF (NODE_TARGET /= COMM%NODE_ME) THEN
        NI_LOCAL=0
      ENDIF
#else
!
! in conventional version all elements are local
!
      NI_LOCAL=NI
#endif

      RETURN
      END FUNCTION
!=======================================================================
!
! this routine gives the global storage index for
! the non local overlap CQIJ ,  strength CDIJ matrix elements
! and for the projected wavefunctions
! return is 0 if the element resides on an other processor
! on entry NI is the local index
!=======================================================================

      FUNCTION NI_GLOBAL(NI,COMM)
      USE prec
      USE  mpimy
      IMPLICIT NONE
      TYPE (communic)  COMM
      INTEGER NI,NI_GLOBAL,NODE_TARGET
#ifdef MPI
      NI_GLOBAL=(NI-1)*COMM%NCPU + COMM%NODE_ME
#else
!
! in conventional version all elements are local
!
      NI_GLOBAL=NI
#endif

      RETURN
      END FUNCTION

!=======================================================================
!  set WDES for one k-point
!  this is quite simple and sometimes necessary
!=======================================================================

      SUBROUTINE CREATE_SINGLE_KPOINT_WDES(WDES_ORIG,WDES,NK)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (wavedes)  WDES,WDES_ORIG

      WDES=WDES_ORIG
      WDES%NKPTS=1

      WDES%NPLWKP=> WDES_ORIG%NPLWKP(NK:NK)
      WDES%NGVECTOR=> WDES_ORIG%NGVECTOR(NK:NK)
      WDES%NPLWKP_TOT=> WDES_ORIG%NPLWKP_TOT(NK:NK)
      WDES%WTKPT => WDES_ORIG%WTKPT (NK:NK)
      WDES%VKPT  => WDES_ORIG%VKPT  (:,NK:NK)
      WDES%NINDPW=> WDES_ORIG%NINDPW(:,NK:NK)
      WDES%IGX   => WDES_ORIG%IGX   (:,NK:NK)
      WDES%IGY   => WDES_ORIG%IGY   (:,NK:NK)
      WDES%IGZ   => WDES_ORIG%IGZ   (:,NK:NK)
!-MM- changes to accommodate spin spirals
! original statement
!     WDES%DATAKE=> WDES_ORIG%DATAKE(:,NK:NK)
      WDES%DATAKE=> WDES_ORIG%DATAKE(:,NK:NK,:)
!-MM- end of alterations
      IF (WDES%NCOL/=0) THEN
        WDES%PL_INDEX => WDES_ORIG%PL_INDEX(:,NK:NK)
        WDES%PL_COL   => WDES_ORIG%PL_COL (:,NK:NK)
      ENDIF

      END SUBROUTINE

!=======================================================================
!  initialize the storage for the wavefunctions
!=======================================================================

      SUBROUTINE ALLOCW(WDES,W,WUP,WDW)
      USE prec
      IMPLICIT NONE
      INTEGER NK
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (wavefun)  WUP,WDW

      INTEGER NRPLWV,NPROD,NKDIM,NBANDS,ISPIN,NB_TOT,NB_PAR,NB_LOW

      NRPLWV=WDES%NRPLWV
      NPROD =WDES%NPROD
      NKDIM =WDES%NKDIM
      NBANDS=WDES%NBANDS
      NB_TOT=WDES%NB_TOT
      NB_LOW=WDES%NB_LOW
      NB_PAR=WDES%NB_PAR
      ISPIN =WDES%ISPIN

      ALLOCATE(W%CPTWFP(NRPLWV,NBANDS,NKDIM,ISPIN), &
               W%CPROJ (NPROD, NBANDS,NKDIM,ISPIN), &
               W%CELTOT(NB_TOT,NKDIM,ISPIN), &
               W%FERTOT(NB_TOT,NKDIM,ISPIN))

      W%CPTWFP=0
      W%CPROJ =0
      W%CELTOT=0
      W%FERTOT=0
      W%FERWE => W%FERTOT(NB_LOW:NB_TOT:NB_PAR,:,:)
      W%CELEN => W%CELTOT(NB_LOW:NB_TOT:NB_PAR,:,:)

      W%OVER_BAND=.FALSE.
!     W%WDES  => WDES

      WUP%CELTOT=> W%CELTOT(:,:,1)
      WUP%FERTOT=> W%FERTOT(:,:,1)
      WUP%CELEN => W%CELEN(:,:,1)
      WUP%FERWE => W%FERWE(:,:,1)
      WUP%CPTWFP=> W%CPTWFP(:,:,:,1)
      WUP%CPROJ => W%CPROJ(:,:,:,1)      
      WUP%OVER_BAND=.FALSE.
!     WUP%WDES  => WDES

      IF (WDES%ISPIN==2) THEN
      WDW%CELTOT=> W%CELTOT(:,:,2)
      WDW%FERTOT=> W%FERTOT(:,:,2)
      WDW%CELEN => W%CELEN(:,:,2)
      WDW%FERWE => W%FERWE(:,:,2)
      WDW%CPTWFP=> W%CPTWFP(:,:,:,2)
      WDW%CPROJ => W%CPROJ(:,:,:,2)
      WDW%OVER_BAND=.FALSE.
!     WDW%WDES  => WDES
      ENDIF

      END SUBROUTINE

!=======================================================================
!  initialize and descriptor for one wavefunction  (wavedes1) from
!  an descriptor of an array of wavefunctions  (wavedes)
!  (kpoint index must be supplied)
!=======================================================================
      SUBROUTINE SETWDES(WDES,WDES1,NK)
      USE prec
      IMPLICIT NONE
      INTEGER NK
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      WDES1%RSPIN= WDES%RSPIN
      WDES1%LNONCOLLINEAR=WDES%LNONCOLLINEAR
      WDES1%NRSPINORS=WDES%NRSPINORS
      WDES1%NRPLWV=WDES%NRPLWV
      WDES1%NGDIM=WDES%NGDIM
      WDES1%NPROD= WDES%NPROD
      WDES1%NBANDS=WDES%NBANDS
      WDES1%NB_TOT=WDES%NB_TOT
      WDES1%NB_PAR=WDES%NB_PAR
      WDES1%NB_LOW=WDES%NB_LOW
      WDES1%NPL  = WDES%NPLWKP(NK)
      WDES1%NGVECTOR= WDES%NGVECTOR(NK)
      WDES1%NPL_TOT= WDES%NPLWKP_TOT(NK)
      WDES1%NPRO = WDES%NPRO
      WDES1%NPRO_TOT = WDES%NPRO_TOT
      WDES1%NIONS= WDES%NIONS
      WDES1%NTYP = WDES%NTYP
      WDES1%NITYP=>WDES%NITYP
      WDES1%LMMAX=>WDES%LMMAX
      WDES1%NPRO_POS=>WDES%NPRO_POS
      WDES1%NINDPW=>WDES%NINDPW(:,NK)
      WDES1%IGX  =>WDES%IGX(:,NK)
      WDES1%IGY  =>WDES%IGY(:,NK)
      WDES1%IGZ  =>WDES%IGZ(:,NK)
!-MM- changes to accommodate spin spirals
! original statement
!     WDES1%DATAKE=>WDES%DATAKE(:,NK)
      WDES1%DATAKE=>WDES%DATAKE(:,NK,:)
      WDES1%LSPIRAL=WDES%LSPIRAL
      WDES1%LZEROZ=WDES%LZEROZ
      WDES1%QSPIRAL=WDES%QSPIRAL
!-MM- end of alteration
      WDES1%NK    =NK
      WDES1%NCOL  =WDES%NCOL
      IF (WDES%NCOL/=0) THEN
        WDES1%PL_INDEX => WDES%PL_INDEX(:,NK)
        WDES1%PL_COL   => WDES%PL_COL (:,NK)
      ENDIF
! can not initialize here
      WDES1%RINPL =1._q
      WDES1%NPLWV =0
      WDES1%NPLWVL=0
      WDES1%COMM        => WDES%COMM
      WDES1%COMM_INTER  => WDES%COMM_INTER
      WDES1%COMM_INB    => WDES%COMM_INB

      END SUBROUTINE

!=======================================================================
!  initialize the optional datas required for real space calculations
!  in a descriptor for one single wavefunction (wavedes1)
!=======================================================================

      SUBROUTINE SETWGRID(WDES1,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d)  GRID
      TYPE (wavedes1) WDES1

      WDES1%RINPL =1._q/GRID%NPLWV
      WDES1%NPLWV =GRID%NPLWV    ! total number of points in real space
      WDES1%NPLWVL=GRID%RL%NP    ! local number of points in real space
      WDES1%MPLWV =GRID%MPLWV    ! dimension of arrays in real space

      END SUBROUTINE

!=======================================================================
!  create storage for one wavefunction W
!  optionally real space array is allocated
!=======================================================================

      SUBROUTINE NEWWAV(W,WDES,MPLWV,ALLOC_REAL)
      USE prec
      IMPLICIT NONE
      TYPE (wavefun1) W
      TYPE (wavedes)  WDES
      LOGICAL ALLOC_REAL
      INTEGER MPLWV

      IF (ALLOC_REAL) THEN
        ALLOCATE(W%CPTWFP(WDES%NRPLWV),W%CPROJ(WDES%NPROD),W%CR(MPLWV))
      ELSE
        ALLOCATE(W%CPTWFP(WDES%NRPLWV),W%CPROJ(WDES%NPROD))
      NULLIFY(W%CR)
      ENDIF
      END SUBROUTINE

!=======================================================================
!  destroy storage for one wavefunction W
!  optionally real space array is deallocated
!=======================================================================

      SUBROUTINE DELWAV(W,DEALLOC_REAL)
      USE prec
      IMPLICIT NONE
      TYPE (wavefun1) W
      LOGICAL DEALLOC_REAL

      IF (DEALLOC_REAL.AND. ASSOCIATED(W%CR)) THEN
      DEALLOCATE(W%CPTWFP,W%CPROJ,W%CR)
      ELSE
      DEALLOCATE(W%CPTWFP,W%CPROJ)
      ENDIF
      END SUBROUTINE

!=======================================================================
!  set one singe wavefunction (W1) from an array of wavefunctions
!  local band index and k point must be supplied
!=======================================================================

      SUBROUTINE SETWAV(W,W1,NB,NK)
      USE prec
      IMPLICIT NONE
      INTEGER NK,NB
      TYPE (wavefun)  W
      TYPE (wavefun1) W1

      W1%CPTWFP=>W%CPTWFP(:,NB,NK)
      W1%CPROJ =>W%CPROJ(:,NB,NK)
      W1%FERWE =W%FERWE(NB,NK)
      W1%CELEN =W%CELEN(NB,NK)

      END SUBROUTINE

!=======================================================================
!  set one singe wavefunction (W1) from an array of wavefunctions
!  local band index and k point must be supplied
!=======================================================================

      SUBROUTINE SETWAV_(W,W1,NB,NK,ISP)
      USE prec
      IMPLICIT NONE
      INTEGER NK,NB,ISP
      TYPE (wavespin) W
      TYPE (wavefun1) W1

      W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
      W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
      W1%FERWE =W%FERWE(NB,NK,ISP)
      W1%CELEN =W%CELEN(NB,NK,ISP)

      END SUBROUTINE

!=======================================================================
!  set wavefunction (wavefun) from an spin array of wavefunctions
!  spin must be supplied
!=======================================================================

      SUBROUTINE SETW_SPIN(W,W1,ISPIN)
      USE prec
      IMPLICIT NONE
      INTEGER ISPIN
      TYPE (wavespin) W
      TYPE (wavefun)  W1

      W1%CPTWFP=>W%CPTWFP(:,:,:,ISPIN)
      W1%CPROJ =>W%CPROJ (:,:,:,ISPIN)
      W1%FERWE =>W%FERWE (:,:,ISPIN)
      W1%CELEN =>W%CELEN (:,:,ISPIN)
      W1%FERTOT=>W%FERTOT(:,:,ISPIN)
      W1%CELTOT=>W%CELTOT(:,:,ISPIN)

      END SUBROUTINE

!************************* SUBROUTINE WVREAL ***************************
!
! this subroutine makes a wavefunction real
! it is required for the gamma point only mode
! to avoid that small non real components develop
!***********************************************************************

      SUBROUTINE WVREAL(WDES,GRID,W)
      USE prec
      USE  mgrid
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (grid_3d)  GRID
#ifdef gammareal
!
! search for position of G=0 zero element
!
      NK  =1
      NIND=0
      NPL =WDES%NGVECTOR(NK)

      DO N=1,NPL
        IND=WDES%NINDPW(N,NK)
        N1= MOD((IND-1),GRID%RC%NROW)+1
        NC= (IND-1)/GRID%RC%NROW+1
        N2= GRID%RC%I2(NC)
        N3= GRID%RC%I3(NC)
        IF (N1==1 .AND. N2==1 .AND. N3==1) THEN
          NIND=N
          EXIT
        ENDIF
      ENDDO

      IF (NIND /= 0) THEN
       DO ISP=1,WDES%ISPIN
       DO NB=1,WDES%NBANDS
         W%CPTWFP(NIND,NB,NK,ISP)= REAL( W%CPTWFP(NIND,NB,NK,ISP) ,KIND=q)
       ENDDO
       ENDDO
       IF (WDES%LNONCOLLINEAR) THEN
          DO ISP=1,WDES%ISPIN
          DO NB=1,WDES%NBANDS
             W%CPTWFP(NIND+NPL,NB,NK,ISP)= REAL( W%CPTWFP(NIND+NPL,NB,NK,ISP) ,KIND=q)
          ENDDO
          ENDDO
       ENDIF
      ENDIF

#endif
      RETURN
      END SUBROUTINE

!=======================================================================
!
! NB_LOCAL returns the local storage index of a band
! if bands are distributed over processors
!
!=======================================================================

      FUNCTION NB_LOCAL(NB,WDES1)
      USE prec
      IMPLICIT NONE
      INTEGER NB,NB_LOCAL
      TYPE (wavedes1)    WDES1

      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
        NB_LOCAL=1+(NB-1)/WDES1%NB_PAR
      ELSE
        NB_LOCAL=0
      ENDIF

      END FUNCTION

!***************************SUBROUTINE WFINIT***************************
!
! this subroutine initializes the wavefunction array W
! it use always a random number generator
! to initialize the coefficients
!
!***********************************************************************

#define serial_compat
#define random

      SUBROUTINE WFINIT(GRID,WDES,W, ENINI,INIWAV)
      USE prec

      USE mpimy
      USE mgrid
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1

      REAL(q) G05CAF
! work arrays
      COMPLEX (q),ALLOCATABLE ::  CW(:)
      COMPLEX (q) :: YY
      INTEGER :: ISPINOR, I, NK

      W%CELTOT=0

      W%CPTWFP=0

      NALLOC = MAXVAL(WDES%NPLWKP_TOT)
      ALLOCATE(CW(NALLOC))

      spin:   DO I=1,WDES%ISPIN
      YR=RANE()  ! what is that ? it is here for compatibility
      kpoint: DO NK=1,WDES%NKPTS

      CALL SETWDES(WDES,WDES1,NK)

#if defined(serial_compat) && defined(MPI)
      band:   DO NB_GLOBAL=1,WDES%NB_TOT
      NB=NB_LOCAL(NB_GLOBAL,WDES1)
!=======================================================================
!     initialize for serial compatible version
!=======================================================================
! use only for testing it is not entirely correct (ENINI is not supported)
       CW=0  ! just in case initialize CW
       NPL_GLOBAL=WDES%NPLWKP_TOT(NK)/WDES%NRSPINORS

       DO ISPINOR=0,WDES%NRSPINORS-1 ! to be included later, Testing only
       DO M=1,NPL_GLOBAL
          YY=RANE()
!  at the gamma it is somethimes better to use random phase factor
!  (!!! but if the cell has inversion symmtry it is a bad choice !!!)
#if defined(random) || defined(gammareal)
          IF (M/=1) YY=CMPLX(REAL(YY,q),0.2*RANE()-0.1)
#endif
          MM=M+NPL_GLOBAL*ISPINOR
          CW(MM)=YY
       ENDDO
       ENDDO
       CALL DIS_PW_BAND(WDES1,NB_GLOBAL,CW,W%CPTWFP(1,1,NK,I))
       IF (NB==0) CYCLE band  ! not on local node, than cycle

       NPL=WDES%NGVECTOR(NK)

       DO ISPINOR=0,WDES%NRSPINORS-1
       DO M=1,NPL
          MM=M+NPL*ISPINOR
!-MM- changes to accommodate spin spirals
! original statement
!          IF(WDES%DATAKE(M,NK)>=ENINI) THEN
          IF(WDES%DATAKE(M,NK,ISPINOR+1)>=ENINI) THEN
             W%CPTWFP(MM,NB,NK,I)=0
          ENDIF
!-MM- changes to accommodate spin spirals
! original statement
!         WW=WDES%DATAKE(M,NK)
          WW=WDES%DATAKE(M,NK,ISPINOR+1)      
!-MM- end of alteration
          IF(WW<=0.000001_q) WW=0.1_q
          W%CPTWFP(MM,NB,NK,I)=W%CPTWFP(MM,NB,NK,I)/WW
       ENDDO
       ENDDO

!=======================================================================
!     initialize components of wavefunction with random numbers
!     weighted with the inverse of the kinetic energy of the plane
!     wave (random phase factor added gK)
!=======================================================================
#else
       NPL=WDES%NGVECTOR(NK)

       NSTEP=1
!       IF (WDES%LNONCOLLINEAR)  NSTEP=2
       ISPINOR=0

       band1:   DO NB=1,WDES%NBANDS,NSTEP

        spinor: DO ISPINOR=0,WDES%NRSPINORS-1 ! to be included later, Testing only
        DO M=1,NPL
           YY=RANE()
!  at the gamma it is somethimes better to use random phase factor
!  (!!! but if the cell has inversion symmtry it is a bad choice !!!)
#if defined(random) || defined(gammareal)
           IF (M/=1) YY=CMPLX(REAL(YY,q),0.2*RANE()-0.1)
#endif

!-MM- changes to accommodate spin spirals
! original statement
!          IF(WDES%DATAKE(M,NK)>=ENINI) THEN
           IF(WDES%DATAKE(M,NK,ISPINOR+1)>=ENINI) THEN
              YY=0
           ENDIF
!-MM- end of alteration
!-MM- changes to accommodate spin spirals
! original statement
!          WW=WDES%DATAKE(M,NK)
           WW=WDES%DATAKE(M,NK,ISPINOR+1)          
!-MM- end of alteration
           IF(WW<=0.000001_q) WW=0.1_q
           YY=YY/WW
           W%CPTWFP(M+ NPL*ISPINOR, NB, NK, I)=YY
       ENDDO
       ENDDO spinor

       IF (NSTEP==2)  THEN   ! to be removed later, Testing only
          DO M=1,NPL
             W%CPTWFP(M+ NPL, NB+1, NK, I) =  W%CPTWFP(M, NB, NK, I)
          ENDDO
       ENDIF
       ENDDO band1

       band:   DO NB=1,WDES%NBANDS
#endif
!=======================================================================
! calculate magnitude squared of wavefunction
!=======================================================================
        WFMAG=0
        DO M=1,WDES%NPLWKP(NK)
           CCC=W%CPTWFP(M,NB,NK,I)
           WFMAG=WFMAG+CCC*CONJG(CCC)
        ENDDO
        CALLMPI( M_sum_d(WDES%COMM_INB, WFMAG, 1))
!=======================================================================
! check that it is nonzero
!=======================================================================
        IF (WFMAG<=0.000001_q) THEN
          WRITE(6,10)
 10       FORMAT('ERROR: WFINIT: wavefunctions linearily dependent at', &
              ' random-number initialization ')
          STOP
        ENDIF
!=======================================================================
!     normalize the wavefunction
!     and set CELEN to kinetic energy
!=======================================================================
        WFMINV=1._q/SQRT(WFMAG)
        DO M=1,WDES%NPLWKP(NK)
          W%CPTWFP(M,NB,NK,I)=W%CPTWFP(M,NB,NK,I)*WFMINV
        ENDDO

        SUM_=0
        SUM2=0

        DO ISPINOR=0,WDES%NRSPINORS-1
        DO M=1,NPL
           MM=M+NPL*ISPINOR
!-MM- changes to accommodate spin spirals
! original statement
!         SUM_=SUM_+W%CPTWFP(MM,NB,NK,I)*CONJG(W%CPTWFP(MM,NB,NK,I))*WDES%DATAKE(M,NK)
          SUM_=SUM_+W%CPTWFP(MM,NB,NK,I)*CONJG(W%CPTWFP(MM,NB,NK,I))*WDES%DATAKE(M,NK,ISPINOR+1)         
!-MM- end of alteration
           SUM2=SUM2+W%CPTWFP(MM,NB,NK,I)*CONJG(W%CPTWFP(MM,NB,NK,I))
        ENDDO
        ENDDO
        CALLMPI( M_sum_d(WDES%COMM_INB, SUM_, 1))
        CALLMPI( M_sum_d(WDES%COMM_INB, SUM2, 1))

        W%CELEN (NB,NK,I)=SUM_
      ENDDO band
      ENDDO kpoint
      ENDDO spin

      CALL MRG_CEL(WDES,W)
#ifdef debug
      #ifdef MPI
      IF (WDES%COMM%NODE_ME ==  WDES%COMM%IONODE) THEN
        WRITE(*,'(I4,F14.10)') (I,REAL(W%CELTOT(I,1,1)),I=1,WDES%NB_TOT)
      ENDIF
      #else
        WRITE(*,'(I4,F14.10)') (I,REAL(W%CELTOT(I,1,1)),I=1,WDES%NB_TOT)
      #endif
#endif

      DEALLOCATE(CW)
      RETURN
      END SUBROUTINE

      END MODULE

!*************************SUBROUTINE GEN_LAYOUT**************************
!
! subroutine GENL_LAYOUT performs a number of tasks:
! ) determines the layout (distribution) of the columns on parallel
!      computers
! also the following  arrays are allocated:
!     GRID%RC%I2 , GRID%RC%I3, GRID%RL%I2, GRID%RL%I3, GRID%RL%INDEX
!     WDES%NPLWKP  WDES%NINDPW
! for LSETUP=.TRUE. the kinetic energy arrays and the G-vector array
! are additionally allocated
!     WDES%DATAKE, WDES%IGX,Y,Z
!
! the data layout is based on the initial reciprocal lattice vectors
! stored in BI
!
!***********************************************************************

      SUBROUTINE GEN_LAYOUT(GRID,WDES, B,BI,IU6,LSETUP)
      USE prec
      USE mgrid
      USE wave
      USE constant
      USE base
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (wavedes)     WDES
      DIMENSION B(3,3),BI(3,3) ! current lattice, and initial lattice
      LOGICAL LSETUP,LUP

      LOGICAL, ALLOCATABLE :: LUSE_IN(:)
      INTEGER, ALLOCATABLE :: USED_ROWS(:),IND2(:),IND3(:)
      INTEGER, ALLOCATABLE :: REDISTRIBUTION_INDEX(:)
      COMMON /WAVCUT/   IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX

      GRID%NGX_rd=xw(GRID%NGX)
      GRID%NGY_rd=GRID%NGY
      GRID%NGZ_rd=zw(GRID%NGZ)

    ! wavefunctions are allways complex in the direct grid in VASP
    ! hence LREAL is set to .FALSE.
      GRID%LREAL=.FALSE.

      GRID%RL%NALLOC=0
      GRID%RC%NALLOC=0
      GRID%IN%NALLOC=0
      WDES%NCOL     =0
!=======================================================================
! determine the layout
! i.e. all required columns
! or (y,z) pairs which are required for the 3d-FFT grid
!=======================================================================
#ifdef MPI
!-----------------------------------------------------------------------
! parallel version
!-----------------------------------------------------------------------
! reciprocal space first
!-----------------------------------------------------------------------
!
! determine all required grid points (usually a circle around center)
!
      ALLOCATE(USED_ROWS(GRID%NGY*GRID%NGZ_rd), &
               IND2(GRID%NGY*GRID%NGZ_rd),IND3(GRID%NGY*GRID%NGZ_rd), &
               REDISTRIBUTION_INDEX(GRID%NGY*GRID%NGZ_rd), &
               LUSE_IN(GRID%NGZ_rd))
! following arrays are required to setup the reciprocal layout
      USED_ROWS=0     ! gives the total number of points in each column
                      ! contained within the cutoff sphere
      NCOL_TOT=0      ! total number of columns
      IND2=0          ! y index of each used column
      IND3=0          ! z index of each used column
! LUSE_IN array is used to find out which z-planes must be contained
! in the intermediate representation
      LUSE_IN=.FALSE.

      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY

        IUSED_MAX=0   ! maxmimal number of rows for this column

        DO NK=1,WDES%NKPTS
          G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
          G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          IUSED=0     ! number of rows for this column and k-point
          DO N1=1,GRID%NGX_rd
            G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
            GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
            GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
            GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI
            ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

            ! exclude some components for gamma-only version (C(G)=C*(-G))
           #ifdef wNGXhalf
            IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
            IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
           #endif
           #ifdef wNGZhalf
            IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
            IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
           #endif
            IF(ENERGI<WDES%ENMAX) THEN
              IUSED=IUSED+1
              LUSE_IN(N3)=.TRUE.
            ENDIF
          ENDDO
          IUSED_MAX=MAX(IUSED,IUSED_MAX)
        ENDDO
        ! if this column is used, enter it in required array
        IF (IUSED_MAX>0) THEN
          NCOL_TOT=NCOL_TOT+1        ! increase column counter
          USED_ROWS(NCOL_TOT)=IUSED  ! set number of rows in this column
          IND2(NCOL_TOT)=N2          ! set y index of this column
          IND3(NCOL_TOT)=N3          ! set z index of this column
        ENDIF
      ENDDO; ENDDO
! setup redistribution array
! by sorting according to the number of rows in each columns
! we try to get better load balancing
      DO I=1,NCOL_TOT
        REDISTRIBUTION_INDEX(I)=I
      ENDDO
! do not sort 1 column (it must hold q=000) !
      CALL SORT_REDIS(NCOL_TOT,USED_ROWS(1),REDISTRIBUTION_INDEX(1))

! distribute among processors

      NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU
      ALLOCATE(GRID%RC%I2( NCOL_MAX ))
      ALLOCATE(GRID%RC%I3( NCOL_MAX ))

      GRID%RC%NFAST= 1
      GRID%RC%NROW = GRID%NGX_rd
      GRID%RC%NCOL = 0

      IND=0
!
! step through all columns and distribute them onto proc.
! in the manner 1 ... NCPU - NCPU ... 1 - 1 ... NCPU - etc.
!
      NODE_TARGET=0  ! NODE onto which column has to go
      LUP=.TRUE.     ! determines whether NODE_TARGET is increased or decreased

      DO NCOL=1,NCOL_TOT
        IND_REDIS=REDISTRIBUTION_INDEX(NCOL)
        N2=IND2(IND_REDIS)
        N3=IND3(IND_REDIS)
        IF (LUP) THEN
          IF (NODE_TARGET == GRID%COMM%NCPU) THEN
            LUP=.FALSE.
          ELSE
            NODE_TARGET=NODE_TARGET+1
          ENDIF
        ELSE
          IF (NODE_TARGET == 1) THEN
            LUP=.TRUE.
          ELSE
            NODE_TARGET=NODE_TARGET-1
          ENDIF
        ENDIF

        IND_ON_CPU=(NCOL-1)/GRID%COMM%NCPU+1

        IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
           GRID%RC%NCOL=GRID%RC%NCOL+1
           IF (IND_ON_CPU /= GRID%RC%NCOL) THEN
              WRITE(*,*)'GENSP: internal error(1) ',GRID%COMM%NODE_ME,IND_ON_CPU,GRID%RC%NCOL
              STOP
           ENDIF
           GRID%RC%I2(IND_ON_CPU)=N2
           GRID%RC%I3(IND_ON_CPU)=N3
        ENDIF
      ENDDO
      GRID%RC%NP    = GRID%RC%NCOL*GRID%RC%NROW
      GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
      WDES%NCOL     = GRID%RC%NCOL
!
! set up intermediate grid
! this usually contains a certain number of z-planes around center
!
      LUSE_IN=.TRUE.

      NCOL_TOT=0
      DO N3=1,GRID%NGX_rd
      DO N2=1,GRID%NGZ_rd
        IF (LUSE_IN(N2)) NCOL_TOT=NCOL_TOT+1
      ENDDO
      ENDDO

      NCOL_MAX=(NCOL_TOT+GRID%COMM%NCPU-1)/GRID%COMM%NCPU

      IF (LPLANE_WISE) THEN
      NCOL_MAX=GRID%NGZ*((GRID%NGX+GRID%COMM%NCPU-1)/GRID%COMM%NCPU)
      ENDIF

      ALLOCATE(GRID%IN%I2( NCOL_MAX ))
      ALLOCATE(GRID%IN%I3( NCOL_MAX ))

      GRID%IN%NFAST= 2
      GRID%IN%NCOL = 0
      GRID%IN%NROW = GRID%NGY

      IND=0
      DO N3=1,GRID%NGX_rd
      DO N2=1,GRID%NGZ_rd
        IF (LUSE_IN(N2)) THEN
          IF (LPLANE_WISE) THEN
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
        ENDIF
      ENDDO
      ENDDO
      IF (GRID%IN%NCOL >NCOL_MAX) THEN
        WRITE(*,*) 'GEN_LAYOUT: internal error',NCOL_MAX,GRID%IN%NCOL
        STOP
      ENDIF
      GRID%IN%NP    = GRID%IN%NCOL*GRID%IN%NROW
      GRID%IN%NALLOC= GRID%IN%NCOL*GRID%IN%NROW
!
! set up real space grid, use all points
! in this case GRID will be conformable to GRID_SOFT in real space
!
      CALL REAL_STDLAY(GRID)

      DEALLOCATE(USED_ROWS,IND2,IND3,REDISTRIBUTION_INDEX,LUSE_IN)
#else
!-----------------------------------------------------------------------
! set reciprocal space and real space layout for non parallel computers
! always x-first (or x-fast) layout
! all grid points are used for FFT
!-----------------------------------------------------------------------
      GRID%RC%NFAST= 1
      GRID%RC%NCOL = GRID%NGZ_rd*GRID%NGY
      GRID%RC%NROW = GRID%NGX_rd
      GRID%RC%NP   = GRID%RC%NCOL*GRID%RC%NROW
      GRID%RC%NALLOC= GRID%RC%NCOL*GRID%RC%NROW
      ALLOCATE(GRID%RC%I2( GRID%RC%NCOL ))
      ALLOCATE(GRID%RC%I3( GRID%RC%NCOL ))
      IND=1
      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY
        GRID%RC%I2(IND)=N2
        GRID%RC%I3(IND)=N3
        IND=IND+1
      ENDDO
      ENDDO

      CALL REAL_STDLAY(GRID)
#endif
!=======================================================================
! count number of plane wave coefficients
! and allocate required arrays
!=======================================================================
      NRPLWV=0

      DO NK=1,WDES%NKPTS
        IND=0
        DO NC=1,GRID%RC%NCOL
        N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
        N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
        DO N1=1,GRID%RC%NROW

        G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))

        GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
        GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
        GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

        ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

! exclude some components for gamma-only version (C(G)=C*(-G))
#ifdef wNGXhalf
        IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
        IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
#endif
#ifdef wNGZhalf
        IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
        IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
#endif
       IF(ENERGI<WDES%ENMAX) THEN
          IND=IND+1
        ENDIF
        ENDDO
        ENDDO
        NRPLWV=MAX(NRPLWV,IND)
      ENDDO
! make WDES%NRPLWV dividable by NB_PAR
      WDES%NRPLWV=((NRPLWV+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
      WDES%NGDIM=WDES%NRPLWV

      WDES%NRPLWV = WDES%NRPLWV*WDES%NRSPINORS

!    CALL  MAKE_STRIDE(WDES%NRPLWV)

      GRID%MPLWV=MAX(GRID%RC%NALLOC ,GRID%IN%NALLOC , GRID%RL%NALLOC)
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      DWRITE 'gen_layout',NODE_ME,GRID%RC%NALLOC ,GRID%IN%NALLOC , GRID%RL%NALLOC
      CALL ALLOCWDES(WDES,LSETUP)
      WDES%NPLWKP=0
      WDES%NGVECTOR=0
     
      IF      (WDES%ISPIN==1  .AND. .NOT. WDES%LNONCOLLINEAR ) THEN
        WDES%NCDIJ=1 
      ELSE IF (WDES%ISPIN==2) THEN
        WDES%NCDIJ=2
      ELSE IF (WDES%ISPIN==1  .AND. WDES%LNONCOLLINEAR ) THEN
        WDES%NCDIJ=4 
      ELSE
        WRITE(*,*) 'internal error: can not set NCDIJ'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE

!=======================================================================
! sorts RA in descending order, and rearanges an index array RB
! seems to be a quicksort, by I am not sure
! subroutine writen by Florian Kirchhof
!=======================================================================

      SUBROUTINE SORT_REDIS(N,RA,RB)
      INTEGER RA(N),RB(N)
      INTEGER RRA,RRB

      IF (N==0) RETURN

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).GT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.GT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END SUBROUTINE

!=======================================================================
! sorts RA in ascending order, and rearanges an index array RB
! seems to be a quicksort, by I am not sure
! subroutine writen by Florian Kirchhof
!=======================================================================

      SUBROUTINE SORT_REDIS_ASC(N,RA,RB)
      INTEGER RA(N),RB(N)
      INTEGER RRA,RRB

      IF (N==0) RETURN

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END SUBROUTINE

!*************************SUBROUTINE COUNT_ROWS ************************
!
! this subroutine counts the total number of plane waves contained
! within the cutoff sphere up to (but excluding) a certain column
! this array is required to find out which index a certain column would
! have in the serial version
! mind the index is 0 based
! the total number of plane waves is returned in NUSED
! USED_POINTS returns the number of plane wave coefficients
! up to column N1,N3
!
!***********************************************************************

      SUBROUTINE  COUNT_ROWS(GRID, WDES, BI, NK, USED_POINTS, NUSED)
      USE prec
      USE mgrid
      USE wave
      USE constant
      USE base
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (wavedes)     WDES
      DIMENSION BI(3,3)

      INTEGER :: USED_POINTS(GRID%NGY,GRID%NGZ)

      NUSED=0
      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY
        G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
        G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
        USED_POINTS(N2,N3)=NUSED

         row: DO N1=1,GRID%NGX_rd

          G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
          GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
          GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
          GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

          ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
         ! exclude some components for gamma-only version (C(G)=C*(-G))
         #ifdef wNGXhalf
          IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
          IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
         #endif
         #ifdef wNGZhalf
          IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
          IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
         #endif
          IF(ENERGI<WDES%ENMAX) THEN
            NUSED=NUSED+1
          ENDIF
         ENDDO row
      ENDDO
      ENDDO
      END SUBROUTINE

!*************************SUBROUTINE REPAD_INDEX_ARRAY  ****************
!
! this subroutine calculates two index array that allow
! to "restore" a plane wave array from an old cutoff and lattice
! to a new one
! this operation works only if the wavefunctions are stored in
! the serial layout (not parallel)
! 
! DO I=INDMAX
!   CW(IND(I))=CWI(INDI(I))
! ENDDO
!***********************************************************************

      SUBROUTINE  REPAD_INDEX_ARRAY(GRID, VKPT, VKPTI, B,  BI, ENMAX, ENMAXI, & 
                NP, NPI, IND, INDI, INDMAX, IFAIL )
      USE prec
      USE mgrid
      USE constant
      USE base
      IMPLICIT NONE

      TYPE (grid_3d)     GRID  ! grid descriptor
      REAL(q) :: VKPT(3)       ! new k-point
      REAL(q) :: VKPTI(3)      ! old k-point
      REAL(q) ::  B (3,3)      ! new reciprocal lattice constant
      REAL(q) ::  BI(3,3)      ! old reciprocal lattice constant
      REAL(q) ::  ENMAX        ! new cutoff
      REAL(q) ::  ENMAXI       ! old cutoff
      INTEGER ::  NP           ! number of plane wave coefficients old
      INTEGER ::  NPI          ! number of plane wave coefficients new
                               ! MIND: NP and NPI must be set by the caller
      INTEGER ::  IND(MAX(NP,NPI))  ! index array new
      INTEGER ::  INDI(MAX(NP,NPI)) ! index array old
      INTEGER ::  INDMAX       ! on return maximum index
      INTEGER ::  IFAIL        ! 0  NP and NPI were correct
                               ! 1  NP was incorrect, 2 NPI was incorrect
! local
      INTEGER NP_,NPI_,N1,N2,N3
      REAL(q) :: G1,G2,G3, G1I,G2I,G3I, GIX,GIY,GIZ, GX,GY,GZ, ENERGI, ENERG

      IFAIL = 0

      NP_ =0
      NPI_=0
      INDMAX=0

      DO N3=1,GRID%NGZ_rd
      DO N2=1,GRID%NGY
        G3=(GRID%LPCTZ(N3)+VKPT(3))
        G2=(GRID%LPCTY(N2)+VKPT(2))
        G3I=(GRID%LPCTZ(N3)+VKPTI(3))
        G2I=(GRID%LPCTY(N2)+VKPTI(2))

         row: DO N1=1,GRID%NGX_rd

          G1=(GRID%LPCTX(N1)+VKPT(1))
          G1I=(GRID%LPCTX(N1)+VKPTI(1))

          GIX= (G1I*BI(1,1)+G2I*BI(1,2)+G3I*BI(1,3)) *TPI
          GIY= (G1I*BI(2,1)+G2I*BI(2,2)+G3I*BI(2,3)) *TPI
          GIZ= (G1I*BI(3,1)+G2I*BI(3,2)+G3I*BI(3,3)) *TPI

          GX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
          GY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
          GZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI

          ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
          ENERG =HSQDTM*( (GX**2)+ (GY**2)+ (GZ**2))

         ! exclude some components for gamma-only version (C(G)=C*(-G))
         #ifdef wNGXhalf
          IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
          IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE
         #endif
         #ifdef wNGZhalf
          IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
          IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
         #endif
          IF (ENERG <ENMAX) THEN
             NP_=NP_+1    ! increase index for new array
          ENDIF
          IF (ENERGI<ENMAXI) THEN
             NPI_=NPI_+1  ! increase index for old array
          ENDIF
          IF (ENERG<ENMAX .AND. ENERGI<ENMAXI) THEN
             INDMAX= MIN( MIN(INDMAX+1 ,NP ), NPI) 
                          ! increase index, and avoid overrun
             IND(INDMAX) =NP_
             INDI(INDMAX)=NPI_
          ENDIF
         ENDDO row
      ENDDO
      ENDDO

      IF  (NP_ /= NP) THEN
         NP=NP_
         IFAIL=1
      ENDIF
      IF  (NPI_ /= NPI) THEN
         NPI=NPI_
         IFAIL=1
      ENDIF

      END SUBROUTINE

      SUBROUTINE REPAD_WITH_INDEX_ARRAY(INDMAX,IND,INDI, CW, CWI)
      USE prec
      IMPLICIT NONE

      INTEGER INDMAX
      INTEGER ::  IND(INDMAX)  ! index array new
      INTEGER ::  INDI(INDMAX) ! index array old
      COMPLEX(q) :: CW(*),CWI(*)
  ! local
      INTEGER I

      DO I=1,INDMAX
         CW(IND(I))=CWI(INDI(I))
      ENDDO

      END SUBROUTINE

!*************************SUBROUTINE GEN_INDEX ************************
!
! subroutine GEN_INDEX calculates the following arrays:
! ) the indexing array NINDPW for copying the plane wave coefficients 
!   from the continuous array CW to the column wise layout used for
!   the 3d-FFT
! for LSETUP=.TRUE., additionally the following array are set up:
! ) the kinetic energies of the plane wave basis states are computed
! ) the G vector corresponding to each plane wave basis state is stored
!
! ) in the parallel version, the arrays PL_INDEX and PL_COL
!   are set up and stored
!     PL_INDEX(NC,NK) stores the position of a column at
!                     which data is stored in the serial version
!     PL_COL(NC,NK)   number of data in this column
!
! the data layout is based on the initial reciprocal lattice vectors
! stored in BI
!
!***********************************************************************


      SUBROUTINE GEN_INDEX(GRID,WDES, B,BI, IU6,IU0,LSETUP)
      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (wavedes)     WDES
      DIMENSION B(3,3),BI(3,3) ! current lattice, and initial lattice
      LOGICAL LSETUP
      INTEGER, ALLOCATABLE :: USED_POINTS(:,:)

#ifdef MPI
      NODE_ME= WDES%COMM%NODE_ME
      IONODE = WDES%COMM%IONODE
#endif
!=======================================================================
! now setup the required quantities
!=======================================================================
      TESTMX=0.0_q

      IXMAX=0
      IYMAX=0
      IZMAX=0
      IXMIN=0
      IYMIN=0
      IZMIN=0
      
      ALLOCATE(USED_POINTS(GRID%NGY,GRID%NGZ))

      kpoint: DO NK=1,WDES%NKPTS
        NLBOXI=0
        IND=1
        CALL COUNT_ROWS(GRID,WDES,BI,NK, USED_POINTS,NUSED)

        IF (WDES%LNONCOLLINEAR) THEN
           NUSED=NUSED*WDES%NRSPINORS
        ENDIF

        col: DO NC=1,GRID%RC%NCOL
        N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
        N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
        IN_THIS_ROW=0

        row: DO N1=1,GRID%RC%NROW
        NLBOXI=NLBOXI+1

        G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
!-MM- changes to accommodate spin spirals
! original statements (the three lines directly below)
        GX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
        GY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
        GZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI
        ENERG =HSQDTM*((GX**2)+(GY**2)+(GZ**2))
      ! kinetic energy of plane wave components of spin up part of the spinor
        GX= ((G1-WDES%QSPIRAL(1)/2)*B(1,1)+(G2-WDES%QSPIRAL(2)/2)*B(1,2)+(G3-WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
        GY= ((G1-WDES%QSPIRAL(1)/2)*B(2,1)+(G2-WDES%QSPIRAL(2)/2)*B(2,2)+(G3-WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
        GZ= ((G1-WDES%QSPIRAL(1)/2)*B(3,1)+(G2-WDES%QSPIRAL(2)/2)*B(3,2)+(G3-WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
        ENERGUP=HSQDTM*((GX**2)+(GY**2)+(GZ**2))
      ! kinetic energy of plane wave components of spin up part of the spinor
        GX= ((G1+WDES%QSPIRAL(1)/2)*B(1,1)+(G2+WDES%QSPIRAL(2)/2)*B(1,2)+(G3+WDES%QSPIRAL(3)/2)*B(1,3)) *TPI
        GY= ((G1+WDES%QSPIRAL(1)/2)*B(2,1)+(G2+WDES%QSPIRAL(2)/2)*B(2,2)+(G3+WDES%QSPIRAL(3)/2)*B(2,3)) *TPI
        GZ= ((G1+WDES%QSPIRAL(1)/2)*B(3,1)+(G2+WDES%QSPIRAL(2)/2)*B(3,2)+(G3+WDES%QSPIRAL(3)/2)*B(3,3)) *TPI
        ENERGDN=HSQDTM*((GX**2)+(GY**2)+(GZ**2))
!-MM- end of alterations
        GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
        GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
        GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI

        ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
        TESTMX=MAX(TESTMX,ENERGI)
       !
       ! exclude some components for gamma-only version (C(G)=C*(-G))
       #ifdef wNGXhalf
        IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) CYCLE row
        IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) CYCLE row
       #endif
       #ifdef wNGZhalf
        IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE row
        IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE row
       #endif
      ! check to see if the kinetic energy of the plane wave is less than
      ! ENMAX in which case the plane wave is included in the set of basis
      ! states for this k point
        IF(ENERGI<WDES%ENMAX) THEN
!       IF ((ENERGUP<WDES%ENMAX).AND.(ENERGDN<WDES%ENMAX)) THEN
          IN_THIS_ROW=IN_THIS_ROW+1

          IXMAX=MAX(IXMAX,GRID%LPCTX(N1))
          IYMAX=MAX(IYMAX,GRID%LPCTY(N2))
          IZMAX=MAX(IZMAX,GRID%LPCTZ(N3))
          IXMIN=MIN(IXMIN,GRID%LPCTX(N1))
          IYMIN=MIN(IYMIN,GRID%LPCTY(N2))
          IZMIN=MIN(IZMIN,GRID%LPCTZ(N3))

          IF (LSETUP) THEN
            WDES%IGX(IND,NK)=GRID%LPCTX(N1)
            WDES%IGY(IND,NK)=GRID%LPCTY(N2)
            WDES%IGZ(IND,NK)=GRID%LPCTZ(N3)
!-MM- changes to accommodate spin spirals
! original statement
!           WDES%DATAKE(IND,NK)=ENERG
            WDES%DATAKE(IND,NK,1)=ENERGUP
            WDES%DATAKE(IND,NK,2)=ENERGDN
!-MM- end of alterations
          ENDIF
          WDES%NINDPW(IND,NK)=NLBOXI
          IND=IND+1
        ENDIF
        ENDDO row
        IF (WDES%NCOL /= 0) THEN
          WDES%PL_INDEX(NC,NK)=USED_POINTS(N2,N3)
          WDES%PL_COL  (NC,NK)=IN_THIS_ROW
        ENDIF
        ENDDO col
!=======================================================================
! check to see if there are less than NRPLWV basis states at this kpoint
! if not stop
!=======================================================================
        IND=IND-1

        ! at this point IND is set to the number of plane wave coefficients
        ! for the current k-point
        IND=IND*WDES%NRSPINORS

        IF(WDES%NRPLWV < IND) THEN
           WRITE(*,*)'internal ERROR: GEN_INDEX: number of plane waves is too large', &
           IND,WDES%NRPLWV
          STOP
        ENDIF
        IF (WDES%NPLWKP(NK)/=0 .AND. WDES%NPLWKP(NK)/=IND) THEN
          WRITE(*,*) 'GEN_INDEX: number of plane waves is incorrect', &
                     ' propably incorrect WAVECAR read in'
          STOP
        ENDIF
        WDES%NPLWKP(NK)=IND
        WDES%NPLWKP_TOT(NK)=IND
        WDES%NGVECTOR(NK)=WDES%NPLWKP(NK)/WDES%NRSPINORS

        CALLMPI( M_sum_i(GRID%COMM,WDES%NPLWKP_TOT(NK)  ,1))

        IF (WDES%NPLWKP_TOT(NK) /= NUSED) THEN
          WRITE(*,*)'internal ERROR 2: GEN_INDEX:',WDES%NPLWKP_TOT(NK),NUSED
          STOP
        ENDIF

        IF (IU6>=0) WRITE(IU6,10)NK,WDES%VKPT(1:3,NK),WDES%NPLWKP_TOT(NK)
      ENDDO kpoint

      DEALLOCATE(USED_POINTS)
!=======================================================================
! write maximum index for each direction and give optimal values for
! NGX NGY and NGZ
!=======================================================================
  10  FORMAT(' k-point ',I2,' :  ',3F6.4,'  plane waves: ',I6)

      NPLMAX=0
      DO NK=1,WDES%NKPTS
        NPLMAX=MAX( WDES%NPLWKP_TOT(NK),NPLMAX)
      ENDDO

      NPLMAX_LOC=0
      NPLMIN_LOC=-NPLMAX
      DO NK=1,WDES%NKPTS
        NPLMAX_LOC=MAX( WDES%NPLWKP(NK),NPLMAX_LOC)
        NPLMIN_LOC=MAX(-WDES%NPLWKP(NK),NPLMIN_LOC)
      ENDDO

      CALLMPI( M_max_i(GRID%COMM,NPLMAX_LOC ,1))
      CALLMPI( M_max_i(GRID%COMM,NPLMIN_LOC ,1))
      NPLMIN_LOC=-NPLMIN_LOC

      IXMIN=-IXMIN
      IYMIN=-IYMIN
      IZMIN=-IZMIN
      CALLMPI( M_max_i(GRID%COMM,IXMAX  ,1))
      CALLMPI( M_max_i(GRID%COMM,IYMAX  ,1))
      CALLMPI( M_max_i(GRID%COMM,IZMAX  ,1))
      CALLMPI( M_max_i(GRID%COMM,IXMIN  ,1))
      CALLMPI( M_max_i(GRID%COMM,IYMIN  ,1))
      CALLMPI( M_max_i(GRID%COMM,IZMIN  ,1))
      IXMIN=-IXMIN
      IYMIN=-IYMIN
      IZMIN=-IZMIN

   IF (IU6>=0) THEN

     #ifdef MPI
      WRITE(IU6,21) NPLMAX_LOC,NPLMIN_LOC
 21   FORMAT(/' maximum and minimum number of plane-waves per node : ',2I6)
     #endif

      WRITE(IU6,20) NPLMAX,IXMAX,IYMAX,IZMAX,IXMIN,IYMIN,IZMIN
  20  FORMAT(/' maximum number of plane-waves: ',I6/ &
     &        ' maximal index in each direction: ',/ &
     &        '   IXMAX=',I3,'   IYMAX=',I3,'   IZMAX=',I3/ &
     &        '   IXMIN=',I3,'   IYMIN=',I3,'   IZMIN=',I3/)

      IWARN=0
      IF (IXMIN==0) IXMIN=-IXMAX
      IF ((IXMAX-IXMIN)*2+1>=GRID%NGX) THEN
        WRITE(IU6,30)'NGX',(IXMAX-IXMIN)*2+2
        IWARN=1
      ELSE
        WRITE(IU6,31)'NGX',(IXMAX-IXMIN)*2+2
      ENDIF

      IF (IYMIN==0) IYMIN=-IYMAX
      IF ((IYMAX-IYMIN)*2+1>=GRID%NGY) THEN
        WRITE(IU6,30)'NGY',(IYMAX-IYMIN)*2+2
        IWARN=1
      ELSE
        WRITE(IU6,31)'NGY',(IYMAX-IYMIN)*2+2
      ENDIF

      IF (IZMIN==0) IZMIN=-IZMAX
      IF ((IZMAX-IZMIN)*2+1>=GRID%NGZ) THEN
        WRITE(IU6,30)'NGZ',(IZMAX-IZMIN)*2+2
        IWARN=1
      ELSE
        WRITE(IU6,31)'NGZ',(IZMAX-IZMIN)*2+2
      ENDIF

      IF (IWARN==1 .AND. IU0>=0 ) &
     &WRITE(IU0,*)'WARNING: wrap around errors must be expected'

  30  FORMAT(' WARNING: wrap around error must be expected', &
     &       ' set ',A3,' to ',I3)
  31  FORMAT(' ',A3,' is ok and might be reduce to ',I3)
   ENDIF

      DWRITE 'gen_index done',NODE_ME,NPLMAX

      RETURN
      END SUBROUTINE
      
      
!-MM- Added to restart spin spiral calculations from a WAVECAR
!     obtained at a different q-vector or using a different
!     value for ENINI
!      
      SUBROUTINE CLEANWAV(WDES,W,ENINI)

      USE prec
      USE constant
      USE wave
      
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1

      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS

#ifdef MPI
      CALL SETWDES(WDES,WDES1,NK)
      band: DO NB_GLOBAL=1,WDES%NB_TOT
      NB=NB_LOCAL(NB_GLOBAL,WDES1)       
        
      IF (NB==0) CYCLE band
          
      NPL=WDES%NGVECTOR(NK)

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1

         DO M=1,NPL
            IF(WDES%DATAKE(M,NK,ISPINOR+1)>=ENINI) W%CPTWFP(M+NPL*ISPINOR,NB,NK,I)=0
         ENDDO
         
      ENDDO spinor
      ENDDO band
#else
      NPL=WDES%NGVECTOR(NK)

      band:   DO NB=1,WDES%NBANDS
      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
       
         DO M=1,NPL
            IF(WDES%DATAKE(M,NK,ISPINOR+1)>=ENINI) W%CPTWFP(M+NPL*ISPINOR,NB,NK,I)=0
         ENDDO

      ENDDO spinor
      ENDDO band
#endif            
      ENDDO kpoint
      ENDDO spin
      
      RETURN
      END SUBROUTINE
!-MM- end of addition
