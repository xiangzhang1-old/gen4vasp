#include "symbol.inc"
#ifndef CACHE_SIZE
#define CACHE_SIZE 8000
#endif
#ifndef MINLOOP
#define MINLOOP 1
#endif
!-----------------------------------------------------------------------
! RCS:  $Id: fftmpi.F,v 1.4 2003/06/27 13:22:18 kresse Exp kresse $
!
!  ADVICE FOR USERS OF "SCALAR MACHINES" (RISC-WORKSTATIONS):
!  ==========================================================
!
!  if you know the cache size please define it on previous line
!  in complex double pression words
!  or somewhere in the preprocessor call
!  if you have no idea at all set to zero
!  in this case each plane will be transformed seperately
!  only for extremely large FFT-grids this might slow the speed
!  some typical CACHE_SIZES
!  IBM 590: 32768
!  SUN    : 16000
!  T3D    : 8000
!
!  You can use SGI SCSL FFT routines (at least for complex complex):
!#define FFT_SCSL
!  these routines provide about 10% performance improvement on Altix
!=======================================================================
!   3-d parallel fast fourier transformation
!   fft-kernels from Juergen Furthmueller (using Cray like FPASSM)
!   communication routines and set of communication routines
!   in fftmpi_map.F written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!=======================================================================

      SUBROUTINE FFTMAKEPLAN(A,GRID)
      USE prec
      USE mgrid
      TYPE (grid_3d) GRID
      REAL(q) A(*)
      RETURN
      END SUBROUTINE

      SUBROUTINE FFTBAS(A,GRID,ISIGN)
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      REAL(q) A(*)

      INTEGER,PARAMETER   :: NMAXM=4096  ! propably sufficient forever
      REAL(q),SAVE:: TRIG1(2*NMAXM),TRIG2(2*NMAXM),TRIG3(2*NMAXM)
#ifdef FFT_SCSL
      REAL(q),SAVE:: SCSL_TRIG1(2*NMAXM),SCSL_TRIG2(2*NMAXM),SCSL_TRIG3(2*NMAXM)
#endif
      INTEGER,SAVE:: IFAC(19,3),NXO=0,NYO=0,NZO=0

      LOGICAL, EXTERNAL :: FACERR

#if defined (use_malloc)
      REAL(q)            :: work(*)
      POINTER               (work_p, work)
#else
      REAL(q),POINTER,SAVE    ::  WORK(:)
#endif
#if defined (use_symmetric_heap)
      INTEGER            :: sym_len
      COMPLEX (q)        :: rcvbuf(*), sndbuf(*)
      POINTER               (rcvbuf_p, rcvbuf), (sndbuf_p, sndbuf)
      COMMON /sym_heap/     rcvbuf_p, sndbuf_p, sym_len
#else
#ifndef T3D_SMA
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
#else
      COMMON /SHM/ MALLOC_DONE, PBUF
      POINTER ( PBUF, BUF ); REAL(q) :: BUF(*)

      POINTER ( RCV_P, RCVBUF )
      POINTER ( SND_P, SNDBUF )
      COMPLEX(q) :: RCVBUF(GRID%MPLWV),SNDBUF(GRID%MPLWV)
      RCV_P = PBUF
      SND_P = PBUF+GRID%MPLWV*2
#endif
#endif

!=======================================================================
! initialization
! check size of work arrays TRIG1,2,3 and allocate WORK
!=======================================================================
      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)
      NMAX=MAX(NX,NY)
      NMAX=MAX(NMAX,NZ)

      IF (NX>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NY>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NZ>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF

      NALLOC=MAX(CACHE_SIZE, MAX(8,4*MINLOOP)* MAX(NX,NY,NZ))
      IF (CACHE_SIZE==0) NALLOC=4*MAX(NX*GRID%RC%NCOL, NY*GRID%IN%NCOL, NZ*GRID%RL%NCOL)

#if defined (use_malloc)
      work_p = malloc(q * nalloc)
#else
      CALL SMART_ALLOCATE_REAL(WORK,NALLOC)
#endif
#if ! defined (use_symmetric_heap)
#ifndef T3D_SMA
      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)
#else
      MALLOC=GRID%MPLWV*4
      IF ( ISHM_CHECK(MALLOC)==0) THEN
        WRITE(*,*)'internal error (1)  in FFTBAS:',MALLOC,MALLOC_DONE
        STOP
      ENDIF
#endif
#endif

!
! Initialize FFT if necessary (changes of mesh size, first call)
!
      IF ((NX/=NXO).OR.(NY/=NYO).OR.(NZ/=NZO)) THEN

#if defined (use_symmetric_heap)
         IF (4 * grid%mplwv .GT. sym_len) THEN
            IF (sym_len .NE. 0) THEN
               CALL shpdeallc (sndbuf_p, info, -1)
               CALL shpdeallc (rcvbuf_p, info, -1)
            ENDIF
            sym_len = q * (4 * grid%mplwv)
            CALLMPI (M_max_i( GRID%COMM, sym_len, 1))
            CALL shpalloc (sndbuf_p, sym_len+1024, info, -1)
            CALL shpalloc (rcvbuf_p, sym_len+1024, info, -1)
         ENDIF
#endif
         IF ((FACERR(NX)).OR.(FACERR(NY)).OR.(FACERR(NZ))) THEN
            WRITE(*,*) 'INIT FFT3D: FACERR'
            STOP
         ENDIF
#ifdef FFT_SCSL
         CALL ZZFFTM(0, NX, 1, 1.0_q, DUMMY, 1, DUMMY, 1, SCSL_TRIG1(1), WORK(1), 1)
         CALL ZZFFTMR(0, NY, 1, 1.0_q, DUMMY, 1, DUMMY, 1, SCSL_TRIG2(1), WORK(1), 1)
         CALL ZZFFTM(0, NZ, 1, 1.0_q, DUMMY, 1, DUMMY, 1, SCSL_TRIG3(1), WORK(1), 1)
#else
         CALL CFTTAB(-NX,IFAC(1,1),TRIG1(1))
         CALL CFTTAB(-NY,IFAC(1,2),TRIG2(1))
         CALL CFTTAB(-NZ,IFAC(1,3),TRIG3(1))
#endif
         ! Remember last mesh size
         NXO=NX
         NYO=NY
         NZO=NZ
      END IF

      IDX=NX
      IDY=NY
      IDZ=NZ
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
       IF (ISIGN ==1) THEN
! transformation along first dimension:
#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NX, GRID%RC%NCOL, 1.0_q, A(1), NX, A(1), NX, &
                     SCSL_TRIG1(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG1(1),IFAC(1,1),2,2*IDX,NX,GRID%RC%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

! transformation along second dimension:
#ifdef FFT_SCSL
         CALL ZZFFTMR (ISIGN, NY, GRID%IN%NCOL, 1.0_q, A(1), GRID%IN%NCOL, &
                     A(1), GRID%IN%NCOL, SCSL_TRIG2(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG2(1),IFAC(1,2),2*GRID%IN%NCOL,2,NY,GRID%IN%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_FORWARD(A(1), GRID%RL%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

! transformation along third dimension:

#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NZ, GRID%RL%NCOL, 1.0_q, A(1), NZ, &
                     A(1), NZ, SCSL_TRIG3(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG3(1),IFAC(1,3),2,2*IDZ,NZ,GRID%RL%NCOL,ISIGN,CACHE_SIZE)
#endif
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
       ELSE
! transformation along third dimension:
#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NZ, GRID%RL%NCOL, 1.0_q, A(1), NZ, &
                     A(1), NZ, SCSL_TRIG3(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG3(1),IFAC(1,3),2,2*IDZ,NZ,GRID%RL%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

! transformation along second dimension:
#ifdef FFT_SCSL
         CALL ZZFFTMR (ISIGN, NY, GRID%IN%NCOL, 1.0_q, A(1), GRID%IN%NCOL, &
                     A(1), GRID%IN%NCOL, SCSL_TRIG2(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG2(1),IFAC(1,2),2*GRID%IN%NCOL,2,NY,GRID%IN%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

! transformation along first dimension:

#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NX, GRID%RC%NCOL, 1.0_q, A(1), NX, A(1), NX, &
                     SCSL_TRIG1(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG1(1),IFAC(1,1),2,2*IDX,NX,GRID%RC%NCOL,ISIGN,CACHE_SIZE)
#endif
      ENDIF

#if defined (use_malloc)
      CALL free (work_p)
#endif

      RETURN
      END SUBROUTINE

!=======================================================================
!   3-d parallel fast fourier transformation
!   fft-kernels from Juergen Furthmueller (using Cray like FPASSM)
!   communication routines and set of communication routines
!   in fftmpi_map.F written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!=======================================================================

      SUBROUTINE FFTBRC(A,GRID,ISIGN)
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      REAL(q) A(*)

      INTEGER,PARAMETER   :: NMAXM=4096  ! propably sufficient forever
      REAL(q),SAVE:: TRIG1(2*NMAXM),TRIG2(2*NMAXM),TRIG3(2*NMAXM)
#ifdef FFT_SCSL
      REAL(q),SAVE:: SCSL_TRIG1(2*NMAXM),SCSL_TRIG2(2*NMAXM),SCSL_TRIG3(2*NMAXM)
#endif
      INTEGER,SAVE:: IFAC(19,3),NXO=0,NYO=0,NZO=0

      LOGICAL, EXTERNAL :: FACERR
#if defined (use_malloc)
      REAL(q)            :: work(*)
      POINTER               (work_p, work)
#else
      REAL(q),POINTER,SAVE    ::  WORK(:)
#endif
#if defined (use_symmetric_heap)
      INTEGER            :: sym_len
      COMPLEX (q)        :: rcvbuf(*), sndbuf(*)
      POINTER               (rcvbuf_p, rcvbuf), (sndbuf_p, sndbuf)
      COMMON /sym_heap/     rcvbuf_p, sndbuf_p, sym_len
#else
#ifndef T3D_SMA
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
#else
      COMMON /SHM/ MALLOC_DONE, PBUF
      POINTER ( PBUF, BUF ); REAL(q) :: BUF(*)

      POINTER ( RCV_P, RCVBUF )
      POINTER ( SND_P, SNDBUF )
      COMPLEX(q) :: RCVBUF(GRID%MPLWV),SNDBUF(GRID%MPLWV)
      RCV_P = PBUF
      SND_P = PBUF+GRID%MPLWV*2
#endif
#endif
!=======================================================================
! initialization
! check size of work arrays TRIG1,2,3 and allocate WORK
!=======================================================================

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)
      NMAX=MAX(NX,NY)
      NMAX=MAX(NMAX,NZ)

      IF (NX>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NY>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NZ>NMAXM) THEN
         WRITE(*,*) ' FFTMPI: Increase NMAXM to ',NMAX
         STOP
      ENDIF

      NALLOC=MAX(CACHE_SIZE, MAX(8,4*MINLOOP)* MAX(NX,NY,NZ))
      IF (CACHE_SIZE==0) NALLOC=4*MAX(NX*GRID%RC%NCOL, NY*GRID%IN%NCOL, NZ*GRID%RL%NCOL)

#if defined (use_malloc)
      work_p = malloc(q * nalloc)
#else
      CALL SMART_ALLOCATE_REAL(WORK,NALLOC)
#endif
#if ! defined (use_symmetric_heap)
#ifndef T3D_SMA
      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)
#else
      MALLOC=GRID%MPLWV*4
      IF ( ISHM_CHECK(MALLOC)==0) THEN
        WRITE(*,*)'internal error (1)  in FFTBAS:',MALLOC,MALLOC_DONE
        STOP
      ENDIF
#endif
#endif
!
! Initialize FFT if necessary (changes of mesh size, first call)
!
      IF ((NX/=NXO).OR.(NY/=NYO).OR.(NZ/=NZO)) THEN

#if defined (use_symmetric_heap)
!         WRITE(*,*)'fftmpi: PE = ', GRID%COMM%NODE_ME, 'MPLWV = ', grid%mplwv
         IF (4 * grid%mplwv .GT. sym_len) THEN
            IF (sym_len .NE. 0) THEN
               CALL shpdeallc (sndbuf_p, info, 0)
               CALL shpdeallc (rcvbuf_p, info, 0)
            END IF
            sym_len = q * (4 * grid%mplwv)
            CALLMPI (M_max_i( GRID%COMM, sym_len, 1))
            CALL shpalloc (sndbuf_p, sym_len+1024, info, -1) 
            CALL shpalloc (rcvbuf_p, sym_len+1024, info, -1) 
         END IF
#endif

         IF ((FACERR(NX)).OR.(FACERR(NY)).OR.(FACERR(NZ))) THEN
            WRITE(*,*) 'INIT FFT3D: FACERR'
            STOP
         ENDIF
#ifdef FFT_SCSL
         CALL ZZFFTM(0, NX, 1, 1.0_q, DUMMY, 1, DUMMY, 1, SCSL_TRIG1(1), WORK(1), 1)
         CALL ZZFFTMR(0, NY, 1, 1.0_q, DUMMY, 1, DUMMY, 1, SCSL_TRIG2(1), WORK(1), 1)
#else
         CALL CFTTAB(-NX,IFAC(1,1),TRIG1(1))
         CALL CFTTAB(-NY,IFAC(1,2),TRIG2(1))
#endif
         CALL RFTTAB(-NZ,IFAC(1,3),TRIG3(1))
         ! Remember last mesh size
         NXO=NX
         NYO=NY
         NZO=NZ
      END IF

      IDX=NX
      IDY=NY
      IDZ=NZ
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
       IF (ISIGN ==1) THEN
! transformation along first dimension:
#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NX, GRID%RC%NCOL, 1.0_q, A(1), NX, A(1), NX, &
                     SCSL_TRIG1(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG1(1),IFAC(1,1),2,2*IDX,NX,GRID%RC%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
#ifdef FFT_SCSL
         CALL ZZFFTMR (ISIGN, NY, GRID%IN%NCOL, 1.0_q, A(1), GRID%IN%NCOL, &
                     A(1), GRID%IN%NCOL, SCSL_TRIG2(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG2(1),IFAC(1,2),2*GRID%IN%NCOL,2,NY,GRID%IN%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_FORWARD(A(1), GRID%RL%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

! transformation along third dimension (complex to real)
         ITRANS=2 ! select complex storage mode on entry, real on exit
         CALL RHFFTM(A(1),WORK(1),TRIG3(1), &
                     IFAC(1,3),1,IDZ+2,NZ,GRID%RL%NCOL,ITRANS,CACHE_SIZE)

!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
       ELSE
! transformation along third dimension:
         ITRANS=-2 ! select complex storage mode on exit
         CALL RHFFTM(A(1),WORK(1),TRIG3(1), &
                     IFAC(1,3),1,IDZ+2,NZ,GRID%RL%NCOL,ITRANS,CACHE_SIZE)
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)

! transformation along second dimension:
#ifdef FFT_SCSL
         CALL ZZFFTMR (ISIGN, NY, GRID%IN%NCOL, 1.0_q, A(1), GRID%IN%NCOL, &
                     A(1), GRID%IN%NCOL, SCSL_TRIG2(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG2(1),IFAC(1,2),2*GRID%IN%NCOL,2,NY,GRID%IN%NCOL,ISIGN,CACHE_SIZE)
#endif
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)

! transformation along first dimension:
#ifdef FFT_SCSL
         CALL ZZFFTM (ISIGN, NX, GRID%RC%NCOL, 1.0_q, A(1), NX, A(1), NX, &
                     SCSL_TRIG1(1),WORK(1),1)
#else
         CALL CFFTML(A(1),A(2),WORK(1), &
                     TRIG1(1),IFAC(1,1),2,2*IDX,NX,GRID%RC%NCOL,ISIGN,CACHE_SIZE)
#endif

      ENDIF

#if defined (use_malloc)
      CALL free (work_p)
#endif

      RETURN
      END SUBROUTINE

!
!===============================================================================
!     3-d fast fourier transformation for wavefunctions
!  this routine must be called only via FFTWAV and FFTEXT
!  if this is not the case results are unpredictable
!  any scaling should be done on the reduced linear grid in the routines
!  FFTWAV and FFTEXT, this saves time
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)    complex to real
!     -1  r->q   vq= sum(r) vr exp(-iqr)    real to complex
!   MIND: result is a complex array !
!===============================================================================

      SUBROUTINE FFT3D(C,GRID,ISN)
      USE prec
      USE  mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)  GRID
#ifdef wNGZhalf
!-------------------------------------------------------------------------------
!  real mode version wWGZhalf
!  in real space the first dimension in VASP is NGZ (COMPLEX data !!)
!  but the FFT requires NGZ+2 (REAL data)
!  therefore some data movement is required
!-------------------------------------------------------------------------------
      REAL(q) C(1)

      NZ=GRID%NGPTAR(3)
!     q->r FFT
      IF (ISN==1) THEN
        CALL FFTBRC(C,GRID,ISN)
!       concat  z-lines (go from stride NZ+2 to NZ)
!DIR$ IVDEP
!OCL NOVREC
        DO IL=GRID%RL%NCOL-1,0,-1
          NDEST=IL* NZ*2
          NSRC =IL*(NZ+2)
!DIR$ IVDEP
!OCL NOVREC
          DO  NZZ=NZ,1,-1
            C(NDEST+NZZ*2-1)=C(NSRC+NZZ)
            C(NDEST+NZZ*2)=0
          ENDDO
        ENDDO
      ELSE
!     r->q FFT
!       z-lines (go from complex stride NZ to real stride NZ+2)
!DIR$ IVDEP
!OCL NOVREC
        DO IL=0,GRID%RL%NCOL-1
        NSRC =IL* NZ*2
        NDEST=IL*(NZ+2)
!DIR$ IVDEP
!OCL NOVREC
          DO NZZ=1,NZ
            C(NDEST+NZZ)=C(NSRC+NZZ*2-1)
          ENDDO
        ENDDO
        CALL FFTBRC(C,GRID,ISN)
      ENDIF
      RETURN
#elif defined(wNGXhalf)
      NOT IMPLEMENTED
#else
!-------------------------------------------------------------------------------
!  complex version
!-------------------------------------------------------------------------------
      COMPLEX(q) C(0:1)
      CALL FFTBAS(C,GRID,ISN)
      RETURN
#endif
      END SUBROUTINE


!************************* SUBROUTINE FFTINI ***************************
!
!  if necessary this routine performes initialization
!  for FFTWAV and FFTEXT
!  usually this is only necessary for the Gamma point only
!  1-kpoint version
!
!   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!               from the FFT grid (FFTEXT)
!   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!               the grid
!***********************************************************************

      SUBROUTINE  FFTINI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)  GRID
      DIMENSION NPLWKP(NKPTS)
      DIMENSION NINDPW(NRPLW,NKPTS)

#if defined(wNGZhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA

      IF (NKPTS>1) THEN
        WRITE(*,*)'FFT3D: real version works only for 1 k-point'
        STOP
      ENDIF

      NK=1
      NPL=NPLWKP(NK)

      ALLOCATE(FFTSCA(NPL,2))

      DO N=1,NPL
        IND=NINDPW(N,NK)
        N1= MOD((IND-1),GRID%RC%NROW)+1
        NC= (IND-1)/GRID%RC%NROW+1
        N2= GRID%RC%I2(NC)
        N3= GRID%RC%I3(NC)

        FACTM=SQRT(2._q)
        IF (N1==1 .AND. N2==1 .AND. N3==1) FACTM=1
        FFTSCA(N,1)= FACTM
        FFTSCA(N,2)= 1/FACTM
        ! this statment seems to be strange but is required
        ! because for z==0 only half of the FFT components are set
        IF (N3==1) FFTSCA(N,2)=FACTM
      ENDDO
#endif
      RETURN
      END SUBROUTINE

!************************* SUBROUTINE FFTWAV ***************************
!
!  this subroutine transforms a wavefunction C defined  within  the
!  cutoff-sphere to real space CR
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGZ != 0 (wNGZhalf)
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be removed before the FFT transformation !
! (scaling with   FFTSCA(M,2))
!
!***********************************************************************

      SUBROUTINE FFTWAV(NPL,NINDPW,CR,C,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      TYPE (grid_3d)     GRID
      DIMENSION C(NPL),CR(GRID%NPLWV)
      DIMENSION NINDPW(NPL)
#if defined(wNGZhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA
#endif

!DIR$ IVDEP
      DO M=1,GRID%RL%NCOL*GRID%NGZ
        CR(M)=(0.0_q,0.0_q)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
#if defined(wNGZhalf)
        CR(NINDPW(M))=C(M)*FFTSCA(M,2)
#else
        CR(NINDPW(M))=C(M)
#endif
      ENDDO
      CALL FFT3D(CR,GRID,1)
      RETURN
      END SUBROUTINE

!************************* SUBROUTINE FFTEXT ***************************
!
! this subroutine performes a FFT to reciprocal space and extracts data
! from the FFT-mesh
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGX != 0 (wNGXhalf) or NGZ != 0 (wNGZhalf)
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be applied after the FFT transformation !
!  (scaling with   FFTSCA(M))
!
!***********************************************************************
      SUBROUTINE FFTEXT(NPL,NINDPW,CR,C,GRID,LADD)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      DIMENSION C(NPL),CR(GRID%NPLWV)
      DIMENSION NINDPW(NPL)
      LOGICAL   LADD
#if defined(wNGZhalf) ||  defined(wNGXhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA
#endif
      CALL FFT3D(CR,GRID,-1)

      IF (LADD) THEN
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
#if defined(wNGZhalf)
          C(M)=C(M)+CR(NINDPW(M))*FFTSCA(M,1)
#else
          C(M)=C(M)+CR(NINDPW(M))
#endif
        ENDDO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
#if defined(wNGZhalf)
          C(M)=CR(NINDPW(M))*FFTSCA(M,1)
#else
          C(M)=CR(NINDPW(M))
#endif
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE


!===============================================================================
!    3-d fast fourier transform (possibly real to complex and vice versa)
!    for chardensities and potentials
!     +1  q->r   vr= sum(q) vq exp(+iqr)    (might be complex to real)
!     -1  r->q   vq= sum(r) vr exp(-iqr)    (might be real to complex)
!
!===============================================================================

      SUBROUTINE FFT3RC(C,GRID,ISN)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)   GRID
      REAL(q) C(1)

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

!-------------------------------------------------------------------------------
!  complex version
!-------------------------------------------------------------------------------
      IF (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd ) THEN

         CALL FFTBAS(C,GRID,ISN)

!-------------------------------------------------------------------------------
!  in real space the first dimension in VASP is NGZ (REAL data)
!  but the FFT required NGZ+2 (real data)
!  therefore some data movement is required
!-------------------------------------------------------------------------------
      ELSE IF (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd ) THEN

         NZ=GRID%NGPTAR(3)
         !     q->r FFT
         IF (ISN==1) THEN
            CALL FFTBRC(C,GRID,ISN)
            !       concat  z-lines (go from stride NZ+2 to NZ)
!DIR$ IVDEP
!OCL NOVREC
            DO IL=0,GRID%RL%NCOL-1
               NDEST=IL* NZ
               NSRC =IL*(NZ+2)
!DIR$ IVDEP
!OCL NOVREC
               DO NZZ=1,NZ
                  C(NDEST+NZZ)=C(NSRC+NZZ)
               ENDDO
            ENDDO
         ELSE
            !     r->q FFT
            !       x-lines (go from stride NZ to NZ+2)
!DIR$ IVDEP
!OCL NOVREC
            DO IL=GRID%RL%NCOL-1,0,-1
               NSRC =IL*NZ
               NDEST=IL*(NZ+2)
!DIR$ IVDEP
!OCL NOVREC
               DO NZZ=NZ,1,-1
                  C(NDEST+NZZ)=C(NSRC+NZZ)
               ENDDO
            ENDDO
            CALL FFTBRC(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  presently not supported
!-------------------------------------------------------------------------------
      ELSE
         WRITE(0,*) 'ERROR in FFT3RC: this version does not support the required half grid mode'
         WRITE(0,*) NX, NY, NZ
         WRITE(0,*) GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd
         STOP
      ENDIF
      RETURN
      END SUBROUTINE


!=======================================================================
!   this routine returns the next correct setting for the
!   three dimensional FFT
!=======================================================================

      SUBROUTINE FFTCHK(NFFT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NFFT(3)
      LOGICAL FFTCH1

      DO 100 IND=1,3
  200 CONTINUE
        IF (FFTCH1(NFFT(IND))) GOTO 100
        NFFT(IND)=NFFT(IND)+1
        GOTO 200
  100 CONTINUE
      END

      LOGICAL FUNCTION FFTCH1(NIN)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (NFACT=4)
      DIMENSION IFACT(NFACT),NCOUNT(NFACT)
      DATA      IFACT /2,3,5,7/
      N=NIN
      DO 100 I=1,NFACT
        NCOUNT(I)=0
  120   NEXT=N/IFACT(I)
        IF (NEXT*IFACT(I)==N) THEN
          N=NEXT
          NCOUNT(I)=NCOUNT(I)+1
          GOTO 120
        ENDIF
  100 CONTINUE
      IF (N==1 .AND. (NCOUNT(1)/=0)) &
     &  THEN
        FFTCH1=.TRUE.
      ELSE
        FFTCH1=.FALSE.
      ENDIF
      RETURN
      END

