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





!===============================================================================
! RCS:  $Id: fft3dsimple.F,v 1.2 2002/08/14 13:59:38 kresse Exp $
!
! this modul implements the following FFT routines (which are called by VASP)
!   FFTWAV
!   FFTEXT
!   FFT3RC
!
! the basic fft routines
!   FFTBAS complex <-> complex
!   FFTBRC complex <-> real
! are missing
! usually this modul should be included in the main fft fortran file
! using the statement
!  #include "fft3dsimple.F"
! the main fft  file should contain the FFTBAS and (optionally FFTBRC) routine;
! but it is also possible to copy this file, and to implement FFTBAS
! (and FFTBRC) directly in the copy.
! On some machines it might be a good idea to optimize all routines
! in this modul, according to the functionallity of the basic FFT routine
! available on the specific machine
! (especially the gamma point only version which is fastest if an
!  real to complex FFT exists requires some optimization)
!
!===============================================================================
!===============================================================================
! template for  basic complex 3-d fast fourier transformation routine
! should give you some idea how the basic complex 3-d fft should be implemented
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!===============================================================================


      SUBROUTINE FFTBAS_(C,N,ISN)
      USE prec
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(0:1)
      DIMENSION N(3)
      WRITE(*,*)"FFTBAS: is not implemented"
      STOP
      RETURN
      END

!===============================================================================
! template for  basic real to complex 3-d fast fourier transformation routine
!   C is used as input and output array
!   in real space C is defined as
!    REAL(q)    C(1:N(1)+2  ,1:N(2),1:N(3)
!   in reciprocal space C is defined as
!    COMPLEX(q) C(1:N(1)/2+1,1:N(2),1:N(3)
!   this conforms to the ESSL and to the CRAY routines
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!===============================================================================

      SUBROUTINE FFTBRC_(C,N,ISN)
      USE prec
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(0:1)
      DIMENSION N(3)
      WRITE(*,*)"FFTBRC: is not implemented"
      STOP
      RETURN
      END

!===============================================================================
!     3-d fast fourier transformation for wavefunctions
!  this routine must be called only via FFTWAV and FFTEXT
!  if this is not the case results are unpredictable
!  any scaling should be done on the reduced linear grid in the routines
!  FFTWAV and FFTEXT, this save time
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)    complex to real
!     -1  r->q   vq= sum(r) vr exp(-iqr)    real to complex
!
!===============================================================================

      SUBROUTINE FFT3D(C,GRID,ISN)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)  GRID
!-------------------------------------------------------------------------------
!  complex version
!-------------------------------------------------------------------------------
      COMPLEX(q) C(0:1)
      CALL FFTBAS(C,GRID%NGPTAR(1),ISN)
      RETURN
      END

!************************* SUBROUTINE FFTINI ***************************
!
!  if necessary this routine performes initialization
!  for FFTWAV and FFTEXT
!  usually this is only necessary for the Gamma point only
!  single  k-point version
!
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
      END


!************************* SUBROUTINE FFTWAV ***************************
!  this subroutine transforms a wavefunction C defined  within  the
!  cutoff-sphere to real space CR
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGX != 0 (wNGXhalf) or NGZ != 0 (wNGZhalf)
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be removed before the FFT transformation !
! (scaling with   FFTSCA(M,2))
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

      DO M=1,GRID%NPLWV
        CR(M)=(0.0_q,0.0_q)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
        CR(NINDPW(M))=C(M)
      ENDDO
! here you can make what you want - safety first ...
      CALL FFT3D(CR,GRID,1)
! ... or 'fastness' first (but often it is not   so   much faster ...):
!      CALL FFTQ2Q(CR,GRID%NGPTAR(1),1)
      RETURN
      END

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
      END


!===============================================================================
!
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

         CALL FFTBAS(C,GRID%NGPTAR,ISN)

!-------------------------------------------------------------------------------
!  real to complex version NGXhalf
!  in real space the first dimension in VASP is NGX (REAL data)
!  but the FFT required NGX+2 (real data)
!  therefore some data movement is required
!-------------------------------------------------------------------------------
      ELSE IF (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd ) THEN


         INC3X=(NX/2+1)*NY
         INC2X= NX/2+1
         !     q->r FFT
         IF (ISN==1) THEN
            CALL FFTBRC(C,GRID%NGPTAR(1),ISN)
            !       concat  x-lines (go from stride NX+2 to NX)
!DIR$ IVDEP
!OCL NOVREC
            DO IL=0,NY*NZ-1
               NDEST=IL*NX
               NSRC =IL*(NX+2)
               !DIR$ IVDEP
               !OCL NOVREC
               DO NXX=1,NX
                  C(NDEST+NXX)=C(NSRC+NXX)
               ENDDO
            ENDDO
         ELSE
            !     r->q FFT
            !     x-lines (go from stride NX to NX+2)
            !
!DIR$ IVDEP
!OCL NOVREC
            DO IL=NY*NZ-1,0,-1
               NSRC =IL*NX
               NDEST=IL*(NX+2)
!DIR$ IVDEP
!OCL NOVREC
               DO NXX=NX,1,-1
                  C(NDEST+NXX)=C(NSRC+NXX)
               ENDDO
            ENDDO
            CALL FFTBRC(C,GRID%NGPTAR(1),ISN)
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

      END SUBROUTINE FFT3RC

!************************ SUBROUTINE MULZ    ***************************
!
!  this subroutine multiplies the Z!=0 components by a factor FACT
!  or divides
!
!***********************************************************************

      SUBROUTINE MULZ(C,NGX,NGY,NGZ,FACT)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(0:NGX-1,0:NGY-1,0:NGZ-1)

      DO N3=1,NGZ/2-1
       DO N2=0,NGY-1
        DO N1=0,NGX-1
          C(N1,N2,N3)= C(N1,N2,N3)*FACT
      ENDDO
      ENDDO 
      ENDDO

      RETURN
      END






!*************************************************************************
!*        Fast Fourier Transform for VASP using the FFTW package.        *
!*                                                                       *
!* This routine is just to give a good interface between VASP and FFTW.  *
!* FFTW3D gives an identical interface to the FFTBAS routine, apart from *
!* that the complex array is placed last instead of first. This is just  *
!* to prevent alignment errors. I.e. if the c-array would be too large   *
!* it is simply cut off with this arrangement (the other way around might*
!* cause trouble on the Macintosh.)                                      *
!* The Fortran calls simply makes use of the FFTW library written in     *
!* C/C++ (see also http://www.fftw.org). The result is considerably      *
!* faster then the fft3dfurth/fftbas system for all array sizes I have   *
!* tested. For small arrays (n<32) you get an improvement of several     *
!* 100%, and for larger you still get an improvement in the range of 50%.*
!* Those who like to use the FFTW should get the package from the url    *
!* above, compile it and link with the library. For the Mac the Absoft   *
!* C/C++ compiler makes a good job even if you have to create your own   *
!* your own makefile.                                                    *
!* It should also be noted that FFTW is just as conservative about mem-  *
!* ory use as FFTBAS, possibly even more so.                             *
!*                                                                       *
!*                                   Ph.D. Stefan Mankefors Mars 3 2000  *
!*                                         sem@fy.chalmers.se            *
!*                                                                       *
!* NOTE: When you compile the library you have to configure the Fortran  *
!* conventions used, i.e. 'upper case', 'lowercase' etc. This is done at *
!* the end of the 'config.h' file where you define *(1._q,0._q)* of the the conv-*
!* ventions used for Fortran. ('config.h' comes with all F77 conventions *
!* undefined.) Please note that this might be machine dependent.         *
!*                                                                       *
!* NOTE2: The real to complex FFT might also be exchanged in a similar   *
!* way. I have not tested this, however, and the gain in time is very    *
!* slim since this type of FFT is used very little by the VASP - at least*
!* as far as the benchmark runs goes (I do not have the experience yet to*
!* claim anything else.) Hence it is a question of an additional library *
!* against a gain of perhaps 1%. I have choosen not to use FFTW for this.*
!* Please observe that this means that fft3dlib still is needed in the   *
!* make process of VASP.                                                 *
!*************************************************************************

       subroutine FFTBAS(c,grid,isign)

       use prec

       implicit none

       include 'fftw3.f'

       Complex(q) c(*), cdummy
       integer grid(3), isign, plan
       integer i,j,k, idummy

       if (isign.le.0) then
        call dfftw_plan_dft_3d(plan,grid(1),grid(2),grid(3),&
                           c, c, & 
                           FFTW_FORWARD, FFTW_ESTIMATE)
       else
        call dfftw_plan_dft_3d(plan,grid(1),grid(2),grid(3),&
                           c, c, & 
                           FFTW_BACKWARD, FFTW_ESTIMATE)
       endif
  
       call dfftw_execute(plan)
       call dfftw_destroy_plan(plan)

       return

       end subroutine

       subroutine FFTMAKEPLAN(c,grid)

       use prec

       implicit none

       include 'fftw3.f'

       integer grid(3),plan
       Complex(q) c(*), cdummy

       call dfftw_plan_dft_3d(plan,grid(1),grid(2),grid(3),&
                           c, c, & 
                           FFTW_FORWARD, FFTW_EXHAUSTIVE)
       call dfftw_destroy_plan(plan)
       call dfftw_plan_dft_3d(plan,grid(1),grid(2),grid(3),&
                           c, c, & 
                           FFTW_BACKWARD, FFTW_EXHAUSTIVE)
       call dfftw_destroy_plan(plan)
       end subroutine


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
      LOGICAL, EXTERNAL :: FFTCHK_FURTH

      FFTCH1=FFTCHK_FURTH(NIN)
      END FUNCTION

      MODULE fft_private
      USE prec
      REAL(q),POINTER,SAVE ::  WORK(:)
      END MODULE


!=======================================================================
!   generic   3-d fast fourier transformation
!   written by Jueregen Furthmueller
!   performes the 3-d real to complex FFT
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)    complex to real
!     -1  r->q   vq= sum(r) vr exp(-iqr)    real to complex
!
!=======================================================================
      SUBROUTINE FFTBRC(A,MF,ISIGN)
      USE prec
      USE fft_private
      USE smart_allocate
      IMPLICIT REAL(q) (A-H,O-Z)

      COMPLEX(q) A(*)
      INTEGER, PARAMETER :: NMAXM=512 ! propably sufficient forever
      DIMENSION TRIGX(2*NMAXM),TRIGY(2*NMAXM),TRIGZ(2*NMAXM)
      DIMENSION IFAC(19,3),MF(3)
      SAVE TRIGX,TRIGY,TRIGZ,IFAC,NXO,NYO,NZO
      DATA NXO /0/, NYO /0/, NZO /0/

      NX=MF(1)
      NY=MF(2)
      NZ=MF(3)
      NMAX=MAX(NX,NY)
      NMAX=MAX(NMAX,NZ)

      NALLOC=MAX(12000, 4*1*MAX(NX,NY,NZ))
      CALL SMART_ALLOCATE_REAL(WORK,NALLOC)

      IF (NX>NMAXM) THEN
         WRITE(*,*) ' FFT3DFURTH: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NY>NMAXM) THEN
         WRITE(*,*) ' FFT3DFURTH: Increase NMAXM to ',NMAX
         STOP
      ENDIF
      IF (NZ>NMAXM) THEN
         WRITE(*,*) ' FFT3DFURTH: Increase NMAXM to ',NMAX
         STOP
      ENDIF
! Initialize FFT if necessary (changes of mesh size, first call)!
      IF ((NX/=NXO).OR.(NY/=NYO).OR.(NZ/=NZO)) THEN
         CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &               TRIGX,TRIGY,TRIGZ,0,0,IERR,12000)
         IF (IERR/=0) THEN
            WRITE(*,*) 'INIT FFT3D: IERR =',IERR
            STOP
         ENDIF
! Remember last mesh sizes!
         NXO=NX
         NYO=NY
         NZO=NZ
      END IF
! Do the transformation!
      CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &            TRIGX,TRIGY,TRIGZ,ISIGN,-ISIGN,IERR,12000)
      IF (IERR/=0) THEN
         WRITE(*,*) 'FFT3D: IERR =',IERR
         STOP
      ENDIF

      RETURN
      END
