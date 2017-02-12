#include "symbol.inc"
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
#ifdef wNGZhalf
!-------------------------------------------------------------------------------
!  real mode version wWGZhalf
!  this version requires only a complex 3d-FFT
!  the inversion through the midpoint is avoided by a trick:
!   the elements with z!=0 have to be multiplied by a factor 2
!   before calling FFT3D
!   in FFT3D the FFT is performed and the result is made real
!   the multiplication with the scaling factor is performed in FFTWAV
!-------------------------------------------------------------------------------
      COMPLEX(q) C(0:GRID%NGPTAR(1)*GRID%NGPTAR(2)*GRID%NGPTAR(3)-1)
      CALL FFTBAS(C,GRID%NGPTAR(1),ISN)
      IF (ISN==1) THEN
        DO M=0,GRID%NGPTAR(1)*GRID%NGPTAR(2)*GRID%NGPTAR(3)-1
            C(M)= REAL( C(M) ,KIND=q)
        ENDDO
      ENDIF
#elif defined(wNGXhalf)
!-------------------------------------------------------------------------------
!  real mode version wWGXhalf
!  this version requires a real 3d-FFT
!-------------------------------------------------------------------------------
      REAL(q) C(1)

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)
      INC3X=(NX/2+1)*NY
      INC2X= NX/2+1

      IF (ISN==1) THEN
        CALL FFTBRC(C,GRID%NGPTAR(1),ISN)
!       concat  x-lines (go from stride NX+2 to NX)
!DIR$ IVDEP
!OCL NOVREC
        DO IL=NY*NZ-1,0,-1
        NDEST=IL*NX*2
        NSRC =IL*(NX+2)
!DIR$ IVDEP
!OCL NOVREC
          DO  NXX=NX,1,-1
            C(NDEST+NXX*2-1)=C(NSRC+NXX)
            C(NDEST+NXX*2)=0
          ENDDO
        ENDDO
      ELSE
!     r->q FFT
!       x-lines (go from complex stride NX to real stride NX+2)
!DIR$ IVDEP
!OCL NOVREC
        DO IL=0,NY*NZ-1
        NSRC =IL*NX*2
        NDEST=IL*(NX+2)
!DIR$ IVDEP
!OCL NOVREC
          DO NXX=1,NX
            C(NDEST+NXX)=C(NSRC+NXX*2-1)
          ENDDO
        ENDDO
        CALL FFTBRC(C,GRID%NGPTAR(1),ISN)
      ENDIF
      RETURN
#else
!-------------------------------------------------------------------------------
!  complex version
!-------------------------------------------------------------------------------
      COMPLEX(q) C(0:1)
      CALL FFTBAS(C,GRID%NGPTAR(1),ISN)
      RETURN
#endif
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

#if defined(wNGZhalf) ||  defined(wNGXhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA

      IF (NKPTS>1) THEN
        WRITE(*,*)'FFT3D: real version works only for 1 k-point'
        STOP
      ENDIF

      NK=1
      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)
      NPL=NPLWKP(NK)

      ALLOCATE(FFTSCA(NPL,2))

      DO N=1,NPL
        IND=NINDPW(N,NKPTS)
        N1= MOD((IND-1),xw(NX))+1
        FACTM=SQRT(2._q)
        IF (N==1) FACTM=1
        FFTSCA(N,1)= FACTM
#ifdef  wNGZhalf
        FFTSCA(N,2)= FACTM
#else
        FFTSCA(N,2)= 1/FACTM
        IF (N1==1) FFTSCA(N,2)= FACTM
#endif
      ENDDO

      RETURN
!-------------------------------------------------------------------------------
!  complex version
!  no initialization necessary
!-------------------------------------------------------------------------------
#else
      RETURN
#endif
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
#if defined(wNGZhalf) ||  defined(wNGXhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA
#endif

      DO M=1,GRID%NPLWV
        CR(M)=(0.0_q,0.0_q)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
#if defined(wNGZhalf) ||  defined(wNGXhalf)
        CR(NINDPW(M))=C(M)*FFTSCA(M,2)
#else
        CR(NINDPW(M))=C(M)
#endif
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
#if defined(wNGZhalf) ||  defined(wNGXhalf)
      REAL(q), POINTER :: FFTSCA(:,:)
      COMMON /FFTSCA/ FFTSCA
#endif

      CALL FFT3D(CR,GRID,-1)

      IF (LADD) THEN
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
#if defined(wNGZhalf) ||  defined(wNGXhalf)
          C(M)=C(M)+CR(NINDPW(M))*FFTSCA(M,1)
#else
          C(M)=C(M)+CR(NINDPW(M))
#endif
        ENDDO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
        DO M=1,NPL
#if defined(wNGZhalf) ||  defined(wNGXhalf)
          C(M)=CR(NINDPW(M))*FFTSCA(M,1)
#else
          C(M)=CR(NINDPW(M))
#endif
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






