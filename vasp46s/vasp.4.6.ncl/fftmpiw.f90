!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define 1             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part         --------------------------------
!
!   charge density: half grid mode X direction
!
!
!   charge density real
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






!-----------------------------------------------------------------------
! RCS:  $Id: fftmpi.F,v 1.3 2002/08/14 13:59:38 kresse Exp $
!
!   3-d parallel fast fourier transformation using fftw
!   written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
! 
!   the FFTBAS_PLAN routine performs both the complex to complex, and
!   complex to real FFT
!
!   the FFTBAS routine is the calling interface for the 
!    complex, complex FFT 
!   whereas FFTBRC is the calling interface for complex to real FFT
!  
!=======================================================================

!
!  this subroutine calls the FFTBAS_PLAN routine  with FFTW_ESTIMATE
!
      SUBROUTINE FFTBAS(A,GRID,ISIGN)
      USE prec
      USE mgrid
      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft
      include 'fftw3.f'

      CALL FFTBAS_PLAN(A,GRID,ISIGN,FFTW_ESTIMATE)

      END SUBROUTINE


      SUBROUTINE FFTBAS_PLAN(A,GRID,ISIGN,IPLAN)
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft
      INTEGER IPLAN   !  make a plan (/=FFTW_ESTIMATE)
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
      INTEGER :: planx, plany, planz

      include 'fftw3.f'
!=======================================================================
! initialization
!=======================================================================
      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)

      IDX=NX
      IDY=NY
      IDZ=NZ

      IF (ISIGN==1) THEN
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_BACKWARD, IPLAN)
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_BACKWARD, IPLAN)
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_c2r(planz, 1, NZ , GRID%RL%NCOL, &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             A(1), NZ, 1, IDZ+2 , &
                             IPLAN)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_BACKWARD, IPLAN)
         ENDIF
      ELSE
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_r2c(planz, 1, NZ , GRID%RL%NCOL, &
                             A(1), NZ, 1, IDZ+2 , &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             IPLAN)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_FORWARD, IPLAN)
         ENDIF
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_FORWARD, IPLAN)
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_FORWARD, IPLAN)
      ENDIF
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
       IF (ISIGN ==1 .AND. IPLAN==FFTW_ESTIMATE) THEN
! transformation along first dimension:
         CALL dfftw_execute(planx)
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
         CALL dfftw_execute(plany)
         CALL MAP_FORWARD(A(1), GRID%RL%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
         CALL dfftw_execute(planz)
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
       ELSE IF(IPLAN==FFTW_ESTIMATE) THEN
! transformation along third dimension:
         CALL  dfftw_execute(planz)
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
         CALL  dfftw_execute(plany)
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
         CALL dfftw_execute(planx)
      ENDIF

      call dfftw_destroy_plan(planx)
      call dfftw_destroy_plan(plany)
      call dfftw_destroy_plan(planz)

      RETURN
      END SUBROUTINE


!=======================================================================
!   3-d parallel real to complex fast fourier transformation using 
!   fftw-kernels
!   communication routines and set of communication routines
!   in fftmpi_map.F written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!=======================================================================

!
!  this subroutine calls the FFTBAS_PLAN routine  with FFTW_ESTIMATE
!
      SUBROUTINE FFTBRC(A,GRID,ISIGN)
      USE prec
      USE mgrid
      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft
      include 'fftw3.f'

      CALL FFTBAS_PLAN(A,GRID,ISIGN,FFTW_ESTIMATE)

      END SUBROUTINE

!===============================================================================
!
!   requires a plan therefore this new calling interface is provided
!  which calles the FFTBAS and FFTBRC routine for generating these
!  plans
!
!===============================================================================

      SUBROUTINE FFTMAKEPLAN(A,GRID)
      USE prec
      USE mgrid
      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft

      include 'fftw3.f'

      CALL FFTBAS_PLAN(A,GRID, 1 ,FFTW_EXHAUSTIVE)
      CALL FFTBAS_PLAN(A,GRID,-1,FFTW_EXHAUSTIVE)

      END SUBROUTINE


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
!-------------------------------------------------------------------------------
!  complex version
!-------------------------------------------------------------------------------
      COMPLEX(q) C(0:1)
      CALL FFTBAS(C,GRID,ISN)
      RETURN
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

      DO M=1,GRID%RL%NCOL*GRID%NGZ
        CR(M)=(0.0_q,0.0_q)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
        CR(NINDPW(M))=C(M)
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
      CALL FFT3D(CR,GRID,-1)

      IF (LADD) THEN
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
          C(M)=C(M)+CR(NINDPW(M))
        ENDDO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
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

