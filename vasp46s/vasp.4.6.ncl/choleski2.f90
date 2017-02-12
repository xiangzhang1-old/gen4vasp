!#define dotiming
!#define debug
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






      MODULE choleski
      USE prec
      CONTAINS
!************************ SUBROUTINE ORTHCH ****************************
! RCS:  $Id: choleski2.F,v 1.5 2002/04/16 07:28:38 kresse Exp $
!
! this subroutine orthonormalizes a set of complex (wave-)functions
! using a Choleski-decomposition of the overlap matrix (O = L L^H)
! in conjunction with inversion of the result of this decomposition
! (U --> U^-1). If we have arbitrary functions {|cfin_i>} on input,
! we have first to build up the overlap matrix OVL_ij=<cfin_i|cfin_j>
! then we have to decompose it (OVL_ij=L_ik U_kj), have to invert
! U_ij and then to form the output set |cfout_i>=U^-1_ji|cfin_j>. As
! (1._q,0._q) can show easily it holds: <cfout_i|cfout_j>=delta_ij !! Voila!
!
!***********************************************************************

      SUBROUTINE ORTHCH(WDES,W, LOVERL,LMDIM,CQIJ,NBLK)
      USE prec
      USE scala
      USE dfast
      USE wave
      USE wave_mpi
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      PARAMETER (NSTRIPD=16)

      LOGICAL   LOVERL
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),ALLOCATABLE,TARGET :: CPROW(:,:),COVL(:,:)
! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:)
      COMPLEX(q)   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:)
      LOGICAL DO_REDIS
      TYPE (REDIS_PW_CTR),POINTER :: H_PW

      NCPU=1
      NODE_ME=0
      IONODE=0
!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      IF (NCPU /= 1) THEN

        DO_REDIS=.TRUE.
        NRPLWV_RED=WDES%NRPLWV/NCPU
        NPROD_RED =WDES%NPROD /NCPU

      ELSE

        DO_REDIS=.FALSE.
        NRPLWV_RED=WDES%NRPLWV
        NPROD_RED =WDES%NPROD

      ENDIF
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

      ! set NSTRIP between [2 and 32]
      NSTRIP=MAX(MIN(NSTRIPD,32/NCPU,NBANDS),1)

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS),COVL(NB_TOT,NB_TOT))

!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
!=======================================================================
      CALL SETWDES(WDES,WDES1,NK)
!   get pointers for redistributed wavefunctions
!   I can not guarantee that this works with all f90 compilers
!   please see comments in wave_mpi.F
      IF (NCPU /= 1) THEN
        CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CPTWFP(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
      ELSE
        CW_RED    => W%CPTWFP(:,:,NK,ISP)
        CPROJ_RED => W%CPROJ(:,:,NK,ISP)
        CPROW_RED => CPROW(:,:)
      ENDIF

!   set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO

      

      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0
!=======================================================================
!  calculate overlap matrix (only upper triangle is needed):
!=======================================================================
      IF (DO_REDIS .AND. LASYNC) THEN
         CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
         DO NPOS=1,NSTRIP
           CALL REDIS_PW_START(WDES, W%CPTWFP(1,NPOS,NK,ISP), NPOS, H_PW)
        ENDDO
      ENDIF

      CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW(1,1))
    ! redistribute everything

      
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
        
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        
        IF (.NOT. LASYNC)  CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
        
      ENDIF

      DO N=1,NB_TOT
      DO I=1,NB_TOT
        COVL(I,N)=(0._q,0._q)
      ENDDO; ENDDO
!
! there is a strange bug in the PII optimized blas DGEMM, which seems
! to access in certain instances data beyond the last array element
! if a matrix is multiplied with a vector (second call to ORTH1)
! to work around this I avoid calling ORTH1 with NB_TOT-NPOS+1=1
      DO NPOS=1,NBANDS-NSTRIP,NSTRIP
        IF (DO_REDIS .AND. LASYNC) THEN
        DO NPOS_=NPOS,NPOS+NSTRIP-1
          CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,NPOS_,NK,ISP), NPOS_, H_PW)
          IF (NPOS_+NSTRIP<=NBANDS) &
          CALL REDIS_PW_START(WDES, W%CPTWFP(1,NPOS_+NSTRIP,NK,ISP), NPOS_+NSTRIP, H_PW)
        ENDDO
        ENDIF

        NPOS_RED  =(NPOS-1)*NCPU+1
        NSTRIP_RED=NSTRIP*NCPU

        CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
          CPROW_RED(1,NPOS_RED),NB_TOT,NBLK, &
          NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
      ENDDO

      IF (DO_REDIS .AND. LASYNC) THEN
      DO NPOS_=NPOS,NBANDS
          CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,NPOS_,NK,ISP), NPOS_, H_PW)
      ENDDO
      ENDIF
      
      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=(NBANDS-NPOS+1)*NCPU

      CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
        CPROW_RED(1,NPOS_RED),NB_TOT,NBLK, &
        NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))

      IF (DO_REDIS .AND. LASYNC) CALL REDIS_PW_DEALLOC(H_PW)

      
      
      

!=======================================================================
! Choleski-decomposition of the overlap matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZTRTRI (inversion):
!=======================================================================
      IF (LscaLAPACK .AND. LscaLU ) THEN
         CALL pPOTRF_TRTRI(WDES%COMM, COVL(1,1),WDES%NB_TOT)

         
      ELSE
         INFO=0
         CALL ZPOTRF &
          & ('U',NB_TOT,COVL(1,1),NB_TOT,INFO)
         IF (INFO/=0) THEN
            WRITE(*,*) 'LAPACK: Routine ZPOTRF failed!',INFO
            STOP
         ENDIF
         CALL ZTRTRI &
          & ('U','N',NB_TOT,COVL(1,1),NB_TOT,INFO)
         IF (INFO/=0) THEN
            WRITE(*,*) 'LAPACK: Routine ZTRTRI failed!',INFO
            STOP
         ENDIF
      ENDIF

      
!=======================================================================
!  construct the orthogonal set:
!=======================================================================

      CALL LINCOM('U',CW_RED(1,1),CPROJ_RED(1,1),COVL(1,1), &
             NB_TOT,NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, &
             NBLK,CW_RED(1,1),CPROJ_RED(1,1))
      

     !  back redistribution
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        IF (LASYNC) THEN
           W%OVER_BAND=.TRUE.
        ELSE
           CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
        ENDIF
      ENDIF
      
!=======================================================================
      ENDDO kpoints
      ENDDO spin
!=======================================================================
      DEALLOCATE(CPROW,COVL)

      RETURN
      END SUBROUTINE
      END MODULE
