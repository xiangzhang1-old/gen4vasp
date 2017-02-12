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





!************************* SUBROUTINE FFT_RC_SCALE *********************
! RCS:  $Id: charge.F,v 1.5 2002/08/14 13:59:37 kresse Exp $
!
! subroutine transforms a real space chargedensity to
! reciprocal space applying the real to complex FFT transformation
!***********************************************************************

      SUBROUTINE FFT_RC_SCALE(CHDENR,CHDEN,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID

      COMPLEX(q) CHDEN(GRID%RC%NP)
      COMPLEX(q)      CHDENR(GRID%RL%NP)

      RINPLW=1.0_q/GRID%NPLWV
      CALL  RL_ADD(CHDENR,RINPLW,CHDENR,0.0_q,CHDEN,GRID)
      CALL FFT3RC(CHDEN,GRID,-1)
      RETURN
      END SUBROUTINE


      MODULE charge
      USE prec

      CONTAINS

!***********************************************************************
!
! subroutine CHSP constructs the electronic charge density according
! to the current wavefunctions and fermi-weights
!
!***********************************************************************

      SUBROUTINE SOFT_CHARGE(GRID,GRID_SOFT,W,WDES, CHDEN)
      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE wave_mpi
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (grid_3d)     GRID,GRID_SOFT
      COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
! work arrays
      COMPLEX(q),   ALLOCATABLE :: CDWORK(:,:)
      COMPLEX(q), ALLOCATABLE,TARGET :: CPTDUM(:)
      INTEGER ISPINOR
      INTEGER, PARAMETER :: NSTRIPD=2
      TYPE (REDIS_PW_CTR),POINTER :: H_PW

      ! MPLWV is the allocation in complex words
      ! hence if CDWORK is REAL (1._q,0._q) needs to double the allocation 
      ALLOCATE(CDWORK(GRID%MPLWV*2,WDES%NCDIJ),CPTDUM(GRID%MPLWV*WDES%NRSPINORS))

      IF (W%OVER_BAND) THEN
         NCPU=1
         NSTRIP=MIN(NSTRIPD,WDES%NBANDS)
         CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
      ENDIF

      CDWORK=0
!=======================================================================
! loop over k-points and bands
!=======================================================================
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS

      IF (W%OVER_BAND) THEN
         DO N=1,NSTRIP
            CALL REDIS_PW_START(WDES, W%CPTWFP(1,N,NK,ISP), N, H_PW)
         ENDDO
      ENDIF
      
      band: DO N=1,WDES%NBANDS

      IF (W%OVER_BAND) THEN
         CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,N,NK,ISP), N, H_PW)
         IF (N+NSTRIP<=WDES%NBANDS) &
         CALL REDIS_PW_START(WDES, W%CPTWFP(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
      ENDIF

      WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
      IF (WEIGHT==0) CYCLE

      NPL=WDES%NGVECTOR(NK)
!=======================================================================
! fourier-transformation of wave-function
! sum up  real space charge density
!=======================================================================
      DO ISPINOR=0,WDES%NRSPINORS-1
         CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM(1+ISPINOR*GRID%MPLWV),W%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
      ENDDO
      spinor: DO ISPINOR=0,WDES%NRSPINORS-1 
      DO ISPINOR_=0,WDES%NRSPINORS-1 
            DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              MM_=M+ISPINOR_*GRID%MPLWV
              CDWORK(M,ISP+ISPINOR_+2*ISPINOR)=CDWORK(M,ISP+ISPINOR_+2*ISPINOR)+CPTDUM(MM)*CONJG(CPTDUM(MM_))*WEIGHT
         ENDDO
      ENDDO
      ENDDO spinor
      ENDDO band
      ENDDO kpoints
      ENDDO spin

      IF (W%OVER_BAND) THEN
         W%OVER_BAND=.FALSE.
         CALL REDIS_PW_DEALLOC(H_PW)
      ENDIF
!=======================================================================
! Fourier-Transformation of charge-density using GRID_SOFT
! only input data from first in-band-group is used, and all 
! are involved in the FFT, final result is distributed among 
! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
!=======================================================================
      DO I=1,WDES%NCDIJ
! now merge the chargedensity from all 
         
         CALL FFT_RC_SCALE(CDWORK(1,I),CHDEN(1,I),GRID_SOFT)
! set the charge-density of unbalanced lattic-vectors to 0
         CALL SETUNB(CHDEN(1,I),GRID_SOFT)
      ENDDO

      DEALLOCATE(CDWORK,CPTDUM)

      RETURN
      END SUBROUTINE



!*************************SUBROUTINE RHOAT0 ****************************
!
!  This routine calculates the term G^2 in the Taylor expansion
!  of the atomic chargedensities multiplies this term with
!  the number of atoms per type and the electronic field constant
!
!***********************************************************************

      SUBROUTINE RHOAT0(P,T_INFO, BETATO,OMEGA)
      USE prec
      USE pseudo
      USE poscar
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)

      BETATO=0._q
      DO 200 NT=1,T_INFO%NTYP
        IF (P(NT)%PSPRHO(1)==0) GOTO 200
        IMAX=4
        DQ=P(NT)%PSGMAX/NPSPTS
        B1=0
        B2=0
        DO 100 I=1,IMAX
          A=(P(NT)%PSPRHO(I+1)-P(NT)%PSPRHO(1))/(DQ*I)**2
          B1=B1+A
          B2=B2+A*A
  100   CONTINUE
        B1=B1/IMAX
        B2=B2/IMAX
        BETATO=BETATO+B1*T_INFO%NITYP(NT)*EDEPS/OMEGA

  200 CONTINUE
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE RHOATO ****************************
!  This routine calculates the chargedensity and its derivatives
!  corresponding to overlapping  atoms
!  LFOUR
!   .FALSE. set up charge density in reciprocal space
!   .TRUE.  set up charge density in real space
!  LPAR
!   .FALSE. use PSPRHO
!   .TRUE.  use PSPCOR
!***********************************************************************

      SUBROUTINE RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CHDER)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%RC%NP),CHDER(GRIDC%RC%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! local variables
      REAL(q), POINTER :: PRHO(:)

      CHTOT=0
      CHDER=0
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
       IF (LPAR) THEN
         PRHO=>P(NT)%PSPCOR
       ELSE
         PRHO=>P(NT)%PSPRHO
       ENDIF
       IF (.NOT.ASSOCIATED(PRHO)) GOTO 200
!=======================================================================
! calculate the scaling factor ARGSC that converts the magnitude of a
! reciprocal lattice vector to the correponding position in the
! pseudopotential arrays
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          CHTOT(N)=CHTOT(N)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N,NT)
          CHDER(N)=CHDER(N)+(T1+REM*(2*T2+REM*3*T3))*ARGSC*CSTRF(N,NT)
        ELSE IF (G==0) THEN
          CHTOT(N)=CHTOT(N)+PRHO(1)*CSTRF(N,NT)
          CHDER(N)=0
        ENDIF
      ENDDO
  200 CONTINUE
      ENDDO typ
!=======================================================================
! set the charge-density of unbalanced lattice-vectors to 0
! and transform the charge-density to real space
!=======================================================================
      CALL SETUNB(CHTOT(1),GRIDC)
      CALL SETUNB(CHDER(1),GRIDC)

      IF (LFOUR) THEN
        CALL FFT3RC(CHTOT,GRIDC,1)
      ENDIF

      RETURN
      END SUBROUTINE
!
!  small subroutine which allocates CWORK explicitly
!  (we could use OPTIONAL, but I am not sure whether this is
!   ok on vector computers)
      SUBROUTINE RHOATO_WORK(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE


      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%MPLWV)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! work arrays
      COMPLEX(q),ALLOCATABLE :: CWORK(:)

      ALLOCATE(CWORK(1:GRIDC%MPLWV))
      CALL RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CWORK)
      DEALLOCATE(CWORK)
      

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE RHOATR ****************************
!  subroutine to calculate the charge-density (or partial-core
!  charge-densities) from overlapping atoms in real space
!  and store it in an real array
!***********************************************************************

      SUBROUTINE RHOPAR(GRIDC,T_INFO,B,P,CSTRF,DENS)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)        DENS(GRIDC%RL%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
! work arrays
      COMPLEX(q),ALLOCATABLE::    CWORK1(:),CWORK2(:)

      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
      CALL RHOATO(.TRUE.,.TRUE.,GRIDC,T_INFO,B,P,CSTRF,CWORK1,CWORK2)
      CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,DENS,GRIDC)

      DEALLOCATE(CWORK1,CWORK2)

      RETURN
      END SUBROUTINE



!*************************SUBROUTINE MAGNET ****************************
!
!  This routine calculates the magnetization density and its derivatives
!  corresponding to overlapping  atoms
!  LFOUR  set up charge density in reciprocal space and store in CHTOT
!
!***********************************************************************

      SUBROUTINE MRHOATO(LFOUR,GRIDC,T_INFO,B,P, CHTOT, ISPIN)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      INTEGER ISPIN
      COMPLEX(q)   CHTOT(GRIDC%MPLWV,ISPIN)
      REAL(q)      B(3,3)
      LOGICAL   LFOUR
      INTEGER ISP
! local variables
      INTEGER NT,NADDR,NI,N,N3,N2,N1,NA(T_INFO%NIONS),NIS,NC
      REAL(q) V1,V2,V3,V4,T0,T1,T2,T3,REM,ARGSC,PSGMA2, &
              GX,GY,GZ,G,ARG,AMAG,FNORM
! work arrays
      COMPLEX(q),ALLOCATABLE:: CSTRF(:)

      IF (ISPIN==0) RETURN

      ALLOCATE(CSTRF(GRIDC%RC%NP))

      CHTOT=0
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
       IF (P(NT)%PSPRHO(1)==0) GOTO 200

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3
      FNORM=1._q/P(NT)%ZVALF

 ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

    ! set up phase factor for this ion
      CALL STUFAK_ONE(GRIDC,1,T_INFO%POSION(1,NI),CSTRF)

    ! AMAG is the magnitude of the  magnetization for this ion
    ! and for the present direction
spin: DO ISP=1,ISPIN

      AMAG=T_INFO%ATOMOM(ISP+(NI-1)*ISPIN)*FNORM

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
        ! calculate the magnitude of the reciprocal lattice vector
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
        ! convert the magnitude of the reciprocal latice vector to a position
        ! in the charge-dens. array  and interpolate the atomic-chargedensity
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=P(NT)%PSPRHO(NADDR-1)
          V2=P(NT)%PSPRHO(NADDR)
          V3=P(NT)%PSPRHO(NADDR+1)
          V4=P(NT)%PSPRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          ! voila: phase factor  * charge * AMAG
          CHTOT(N,ISP)=CHTOT(N,ISP)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N)*AMAG
        ELSE IF (G==0) THEN
          CHTOT(N,ISP)=CHTOT(N,ISP)+P(NT)%PSPRHO(1)*CSTRF(N)*AMAG
        ENDIF
      ENDDO
      ENDDO spin
      ENDDO ion
  200 NIS=NIS+T_INFO%NITYP(NT)
!=======================================================================
      ENDDO typ
!=======================================================================

! set the charge-density of unbalanced lattic-vectors to 0
! and transform the charge-density to real space

      DO ISP=1,ISPIN
         CALL SETUNB(CHTOT(1,ISP),GRIDC)
         IF (LFOUR) THEN
            CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
         ENDIF
      ENDDO

      DEALLOCATE( CSTRF )

      RETURN
      END SUBROUTINE
      END MODULE

!*************************SUBROUTINE RHO0  *****************************
! this subroutine calculates the total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      FUNCTION RHO0(GRID, CHTOT)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO_SUM=CHTOT(NC,N1)
        ENDIF
      ENDDO

      
      RHO0= RHO_SUM
      RETURN
      END FUNCTION

!*************************SUBROUTINE RHO0  *****************************
! this subroutine calculates the total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      SUBROUTINE GET_RHO0(GRID, CHTOT, RHO0)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)
      N2= 1
      N3= 1
      N1= 1
      RHO=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO=CHTOT(NC,N1)
        ENDIF
      ENDDO

      
      RHO0= RHO
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE SET_RHO0  *************************
! this subroutine sets the  total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      SUBROUTINE SET_RHO0(GRID, CHTOT, RHO_SOLL)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          CHTOT(NC,N1)=RHO_SOLL
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE





!*************************SUBROUTINE TRUNC_HIGH_FREQU  *****************
!
! this subroutine truncates the high frequency components
! of an array by imposing a spherical cutoff
! the spherical cutoff is chosen in such a way that the entire sphere
! fits into the parallelepided spanned by the plane wave basis set
! this subroutine is required for GGA calculations to maintain the
! full symmetry of the exchange correlation potential
!
!***********************************************************************

      MODULE compat_gga
        LOGICAL GGA_COMPAT
      END MODULE compat_gga

      SUBROUTINE TRUNC_HIGH_FREQU(LATT_CUR, GRID, C)
      USE prec
      USE mgrid
      USE lattice
      USE compat_gga
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt)    LATT_CUR
      COMPLEX(q)     C(GRID%MPLWV)
      REAL(q) GX, GY, GZ, G2, GMIN2
      INTEGER I, N1, N2, N3, NC, NZERO
!
! if you want to revert to the behaviour of vasp before vasp.4.6.16
! comment in the RETURN statment below

      IF (GGA_COMPAT) RETURN

      GX=GRID%NGX/2/LATT_CUR%ANORM(1)
      GY=GRID%NGY/2/LATT_CUR%ANORM(2)
      GZ=GRID%NGZ/2/LATT_CUR%ANORM(3)

      GMIN2=MIN(GX,GY,GZ)**2
      NZERO=0

      DO I=1,GRID%RC%NP
         N1= MOD((I-1),GRID%RC%NROW) +1
         NC= (I-1)/GRID%RC%NROW+1
         N2= GRID%RC%I2(NC)
         N3= GRID%RC%I3(NC)
         GX=(GRID%LPCTX(N1)*LATT_CUR%B(1,1)+GRID%LPCTY(N2)*LATT_CUR%B(1,2)+GRID%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRID%LPCTX(N1)*LATT_CUR%B(2,1)+GRID%LPCTY(N2)*LATT_CUR%B(2,2)+GRID%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRID%LPCTX(N1)*LATT_CUR%B(3,1)+GRID%LPCTY(N2)*LATT_CUR%B(3,2)+GRID%LPCTZ(N3)*LATT_CUR%B(3,3))
         G2=GX*GX+GY*GY+GZ*GZ

         IF (G2>GMIN2) THEN
            NZERO=NZERO+1
            C(I)=0
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE





!****************** subroutine GGA_COMPAT_MODE ***********************
!
! If GGA_COMPAT is .TRUE. the vasp.4.4-4.6 behavior is used whereas
! for GGA_COMPAT .FALSE. the new corrected version is used
!
! GGA_COMPAT = .TRUE.   is the default for vasp.4.6
! GGA_COMPAT = LCOMPAT  is the default for vasp.5.0
!
!***********************************************************************

      SUBROUTINE GGA_COMPAT_MODE(IU5, IU0, LCOMPAT)
        USE prec
        USE compat_gga
        IMPLICIT NONE 
        
        INTEGER               :: IU5,IU0    ! input unit
        LOGICAL :: LCOMPAT

        LOGICAL :: LOPEN,LDUM
        INTEGER :: IDUM, N, IERR
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER*1 :: CHARAC
        
        LOPEN=.FALSE.
        OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')

        GGA_COMPAT=.TRUE.
        CALL RDATAB(LOPEN,'INCAR',IU5,'GGA_COMPAT','=','#',';','L', &
             &  IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,N,1,IERR)

        IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'Error reading item ''GGA_COMPAT'' from file INCAR.'
           ENDIF
           GGA_COMPAT=.FALSE.
        ENDIF

        CALL XML_INCAR('GGA_COMPAT','L',IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,N)

      END SUBROUTINE GGA_COMPAT_MODE

      SUBROUTINE XML_WRITE_GGA_COMPAT_MODE
        USE prec
        USE compat_gga
        IMPLICIT NONE
        INTEGER :: IDUM
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER*1 :: CHARAC

        CALL XML_INCAR('GGA_COMPAT','L',IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,1)

      END SUBROUTINE XML_WRITE_GGA_COMPAT_MODE
