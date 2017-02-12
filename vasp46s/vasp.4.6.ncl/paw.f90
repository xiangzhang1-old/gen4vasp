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





!*******************************************************************
! RCS:  $Id: paw.F,v 1.14 2003/06/27 13:22:22 kresse Exp kresse $
!
!  PAW-module
!  implements the top level PAW functions
!  most of the functionality is implemented in the MODULS radial
!
!  all routine written by Georg Kresse
!  (even the symmetrization routines :)
!*******************************************************************
  MODULE paw
    USE prec
    USE pseudo
      ! LMAX_MIX specifies the maximal L component mixed in the
      ! broyden mixer, the number can be customized by the users
      ! set to a large value if all components should be mixed,
      ! it seems the L=2 is sufficient, and acutally even faster
      ! than mixing all L
      INTEGER  :: SET_LMAX_MIX_TO=2

      INTEGER LMAX_MIX   ! this is the maximum L mixed by the Broyden mixer
      INTEGER IONS_LOCAL ! here we store the number of ions treated locally
                         ! in PAW
      LOGICAL, ALLOCATABLE :: DO_LOCAL(:)
      REAL(q), ALLOCATABLE :: METRIC(:)
      LOGICAL MIMIC_US
      REAL(q)       DOUBLEC_PS_ATOM,DOUBLEC_AE_ATOM

    CONTAINS

!*******************************************************************
!
!  start up procedure for PAW
!  checks internal consistency of the QPAW
!  and sets up the compensation charges
!  on the radial grid, and spline coefficients which are used
!  to interpolate compensation charges in us.F
!
!*******************************************************************

      SUBROUTINE SET_AUG(NTYPD, P, IU6, LEXCH, LEXCHG, LMAX_CALC, LMETAGGA, LCOMPAT)
        USE radial
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER NTYPD
        INTEGER IU6
        INTEGER LEXCH,LEXCHG
        TYPE (potcar),TARGET :: P(NTYPD)
        TYPE (rgrid),POINTER :: R
        INTEGER LMAX_CALC  ! if -1 mimic US PP (see below)
      ! local variables
        INTEGER NTYP,CHANNELS,LMAX,L,LP,N,I
        INTEGER, PARAMETER :: NQ=2
        REAL(q)  QQ(NQ)
        REAL(q)  A(NQ),B,ALPHA
        REAL(q)  SUM,QR,BJ,STEP,X,DEXC,DEXCM
        REAL(q),PARAMETER ::  TH=1E-6_q
        LOGICAL LMETAGGA, LCOMPAT


      ! if LMAX_CALC == -1 the PAW will run in a special mode
      ! resulting in essentially US-PP like behaviour
      ! i.e. all terms are linearized around the atomic reference configuration

        IF (LMAX_CALC ==-1) THEN
           MIMIC_US=.TRUE.
        ELSE
           MIMIC_US=.FALSE.
        ENDIF

        LMAX_MIX=0

        typ: DO NTYP=1,NTYPD
           IF (.NOT. ASSOCIATED(P(NTYP)%QPAW)) CYCLE

           ! maximal L mixed in Broyden mixer
           LMAX_MIX=MAX(LMAX_MIX,P(NTYP)%LMAX_CALC)

           R => P(NTYP)%R
           CALL  RAD_ALIGN(R)        ! reallign RMAX with grid

           CALL RAD_CHECK_QPAW( R, P(NTYP)%LMAX, &
                P(NTYP)%WAE, P(NTYP)%WPS, P(NTYP)%QPAW, P(NTYP)%QTOT , P(NTYP)%LPS)
!OBengone modify start
!          CALL PHI_DOT_PHI(R, NTYP, P(NTYP)%LMAX, P(NTYP)%WAE, P(NTYP)%WPS, P(NTYP)%LPS) 
           IF ((USELDApU().OR.LCALC_ORBITAL_MOMENT()).AND.INTEGRALS_LDApU()) THEN
              CALL OVERLAP_AE(R,NTYP,NTYPD,P(NTYP)%LMAX,P(NTYP)%WAE,P(NTYP)%LPS)
           ENDIF
!OBengone modify end

           CHANNELS=P(NTYP)%LMAX
           DO L=1, CHANNELS
              DO LP=1, P(NTYP)%LMAX
                 IF (P(NTYP)%LPS(L)==P(NTYP)%LPS(LP)) THEN
                    P(NTYP)%QION(L,LP)=P(NTYP)%QPAW(L,LP,0)
                 ENDIF
              ENDDO
           ENDDO

           LMAX=0
           DO I=1,CHANNELS
              LMAX=MAX( P(NTYP)%LPS(I),LMAX )
           ENDDO

           LMAX=LMAX*2               ! maximum l in augmentation charges
           ALLOCATE(P(NTYP)%QDEP(NPSRNL,5,0:LMAX), &
                    P(NTYP)%AUG (R%NMAX,0:LMAX) )

!           IF (IU6>=0) WRITE(IU6,1) NTYP
1          FORMAT(' L augmenation charges for type=',I4)
2          FORMAT(I2,' q=',2F10.6,'  a=',2F10.6)

           ll: DO L=0,LMAX

        ! find q values
              CALL AUG_SETQ(L,R,QQ,A,LCOMPAT)
!              IF (IU6>=0) WRITE(IU6,2) L,QQ,A

        ! setup augmentation charge on radial grid  rho(r) r^2

              DO N=1,R%NMAX
                 SUM=0
                 IF (R%R(N) <= R%RMAX) THEN
                    DO I=1,NQ
                       QR=QQ(I)*R%R(N)
                       CALL SBESSEL( QR, BJ, L)
                       SUM=SUM+BJ*A(I)*R%R(N)*R%R(N)
                    ENDDO
                 ENDIF
                 P(NTYP)%AUG(N,L)=SUM
              ENDDO

        ! setup spline for augmentation charge
              ! the spline ends at PSDMAX*(NPSRNL-1)/NPSRNL see SETDEP
              STEP= R%RMAX/(NPSRNL-1)
              P(NTYP)%PSDMAX=NPSRNL*STEP

              DO N=1,NPSRNL
                 X=STEP*(N-1)
                 SUM=0
                 DO I=1,NQ
                    QR=QQ(I)*X
                    CALL SBESSEL( QR, BJ, L)
                    SUM=SUM+BJ*A(I)
                 ENDDO
                 P(NTYP)%QDEP(N,1,L) = X
                 P(NTYP)%QDEP(N,2,L) = SUM
              ENDDO
              ! derivative at startpoint
              X=STEP/100
              SUM=0
              DO I=1,NQ
                 QR=QQ(I)*X
                 CALL SBESSEL( QR, BJ, L)
                 SUM=SUM+BJ*A(I)
              ENDDO
              SUM=(SUM-P(NTYP)%QDEP(1,2,L))/X
              CALL SPLCOF(P(NTYP)%QDEP(1,1,L),NPSRNL,NPSRNL,SUM)
           ENDDO ll

        ! finally calculate the exchange energy of the core charge
        ! the contribution from outside the core radius (R%NMAX)
        ! was read from the POTCAR file (entry DEXC and stored
        ! in P%DEXCCORE by the POTCAR reader

           CALL RAD_CORE_XC( P(NTYP)%R, LEXCH, LEXCHG, P(NTYP)%RHOAE, &
                DEXC)
        ! the same for metaGGA
           IF (LMETAGGA) THEN
              CALL RAD_CORE_META_XC( P(NTYP)%R, P(NTYP)%RHOAE, P(NTYP)%TAUAE, DEXCM)
           ELSE
              DEXCM=0
           ENDIF

           IF (P(NTYP)%DEXCCORE == -1) THEN
              WRITE(0,*) P(NTYP)%DEXCCORE, ' set to 0';  P(NTYP)%DEXCCORE=0
              P(NTYP)%DEXCCOREM=0
	   ELSE
              P(NTYP)%DEXCCOREM=P(NTYP)%DEXCCORE+DEXCM
              P(NTYP)%DEXCCORE=P(NTYP)%DEXCCORE+DEXC
!              WRITE(0,*)'core exchange correlation GGA, MGGA',P(NTYP)%DEXCCORE,P(NTYP)%DEXCCOREM
           ENDIF


        ENDDO typ

 ! mix only L=SET_LMAX_MIX_TO components in Broyden mixer
        LMAX_MIX=MIN(LMAX_MIX,SET_LMAX_MIX_TO)

 ! if US-PP are mimiced no need to mix any onsite components
        IF (MIMIC_US) LMAX_MIX=-1

      END SUBROUTINE SET_AUG

!*******************************************************************
!
! SET_RHO_PAW_ELEMENTS
! calculates the number of elements of the augmentation
! occupancies which must be mixed on the local node
! also allocates and sets the mask array DO_LOCAL which determines
! which ions are treated localy
!
!*******************************************************************

    SUBROUTINE SET_RHO_PAW_ELEMENTS(WDES, P , T_INFO, LOVERL, ELEMENTS )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      USE asa
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      TYPE (wavedes)   WDES
      LOGICAL  LOVERL

    ! local variables
      TYPE (potcar),POINTER::  PP
      INTEGER NT,NI,NIP,LMAX_TABLE,NODE_TARGET,LM,LMP,LL,LLP,LMIN,LMAX
      INTEGER LMAIN,ELEMENTS,CH1,CH2,NI_PAW_LOCAL
      REAL(q) :: AMETRIC
!=======================================================================
! quick return and allocation of work space
!=======================================================================
      ELEMENTS=0

      IF (.NOT.LOVERL) RETURN

      ! allocate the DO_LOCAL array
      ALLOCATE (DO_LOCAL(T_INFO%NIONS))

      DO_LOCAL=.FALSE.
      IONS_LOCAL=0
!=======================================================================
! cycle all ions and calculate number of L=0 channels
!=======================================================================
      ELEMENTS=0

      ion: DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         NT =T_INFO%ITYP(NI)
    ! not on local node, cycle
    ! no PAW for this ion, cycle
         IF ( NIP/=0 .AND. ASSOCIATED(P(NT)%QPAW) ) THEN
            IONS_LOCAL=IONS_LOCAL+1

            ! distribute ions in a round robin fashion
            ! between  that share (1._q,0._q) DIJ
            DO_LOCAL(NI)= .TRUE.

            IF (DO_LOCAL(NI) ) THEN
            ! how many elements are there
               DO CH1=1,P(NT)%LMAX
                  DO CH2=CH1,P(NT)%LMAX
                     ! quantum numbers l and lp of these two channels
                     LL =P(NT)%LPS(CH1)
                     LLP=P(NT)%LPS(CH2)
                     ! Lmin and Lmax
                     LMIN=ABS(LL-LLP) ; LMAX=MIN(LMAX_MIX,P(NT)%LMAX_CALC,ABS(LL+LLP))
                     DO LMAIN=LMIN,LMAX,2
                        ELEMENTS=ELEMENTS+LMAIN*2+1
                        ELEMENTS=ELEMENTS+LMAIN*2+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO ion
!=======================================================================
! finally calculate the metric tensor of each component that is mixed
! locally
!=======================================================================
      IF (.NOT. ALLOCATED( METRIC) ) THEN
         ALLOCATE( METRIC(ELEMENTS))
      ENDIF

      NI_PAW_LOCAL=0
      DO NI=1,T_INFO%NIONS
         IF (.NOT. DO_LOCAL(NI)) CYCLE

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)
         DO CH1=1,P(NT)%LMAX
            DO CH2=CH1,P(NT)%LMAX
               ! quantum numbers l and lp of these two channels
               LL =P(NT)%LPS(CH1)
               LLP=P(NT)%LPS(CH2)
               ! Lmin and Lmax
               LMIN=ABS(LL-LLP) ; LMAX=MIN(LMAX_MIX,P(NT)%LMAX_CALC,ABS(LL+LLP))
               DO LMAIN=LMIN,LMAX,2
                 ! metric is defined as
                 ! M = ABS( rho(ae)^2 - (rho(ps)+rho(comp))^2)
                 CALL RAD_METRIC( PP%R, PP%AUG(:,LMAIN), PP%QPAW(CH1,CH2,LMAIN), &
                  PP%WAE(:,CH1), PP%WAE(:,CH2), &
                  PP%WPS(:,CH1), PP%WPS(:,CH2), AMETRIC)
                 ! generally larger L qauntum numbers are less important
                 AMETRIC=SQRT(MIN(MAX(ABS(AMETRIC),1E-6_q),1E-2_q))
                 ! if you want no metric, it will be slower ...
                 METRIC( NI_PAW_LOCAL+1:NI_PAW_LOCAL+LMAIN*2+1 ) = AMETRIC
                 METRIC( NI_PAW_LOCAL+1+ELEMENTS/2:NI_PAW_LOCAL+LMAIN*2+1+ELEMENTS/2 ) = AMETRIC
                 NI_PAW_LOCAL=NI_PAW_LOCAL+LMAIN*2+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO

    END SUBROUTINE SET_RHO_PAW_ELEMENTS

!*******************************************************************
!
! SET_RHO_PAW calculates the elements of the augmentation
! occupancies CRHODE that must be mixed on the local node and
! stores them in RHOLM
!
!  calling convention for this must be changed to
!  (charge, m) instead of (up, down)
!
!*******************************************************************

    SUBROUTINE SET_RHO_PAW(WDES, P , T_INFO, LOVERL, &
         ISPIN, LMDIM, CRHODE , RHOLM_STORE )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::      P(T_INFO%NTYP)
      TYPE (wavedes)     WDES
      INTEGER LMDIM,ISPIN
      COMPLEX(q)  CRHODE(LMDIM,LMDIM,WDES%NIONS,ISPIN)
      LOGICAL  LOVERL
      REAL(q) RHOLM_STORE(:,:)

    ! local variables
      TYPE (potcar),POINTER:: PP
      INTEGER NT,NI,NIP,ISP,ITMP
      INTEGER ISIZE,IBASE,IADD
      INTEGER, EXTERNAL :: MAXL1
      REAL(q) RHOLM(LMDIM*LMDIM)

!=======================================================================
! quick return if possible
!=======================================================================
      IF (.NOT.LOVERL .OR. IONS_LOCAL == 0 ) RETURN

      IF (MIMIC_US) RETURN
!=======================================================================
! cycle all ions and set required elements
!=======================================================================
      IBASE=1
      ISIZE=UBOUND(RHOLM_STORE,1) ! runtimecheck for insufficient allocation

      ion: DO NI=1,T_INFO%NIONS
         IF (.NOT. DO_LOCAL(NI)) CYCLE ion

         NIP=NI_LOCAL(NI, WDES%COMM_INB) ! local storage index
         IF (NIP==0) THEN
            WRITE(0,*) 'SET_RHO_PAW: internal error: ion not local'
            STOP
         ENDIF

         NT=T_INFO%ITYP(NI)

         PP=> P(NT)
         DO ISP=1,ISPIN
    ! transform CRHODE (lm,lpmp) to llp,LM
            CALL TRANS_RHOLM( CRHODE(:,:,NIP,ISP), RHOLM, PP )
    ! retrieve elements which are mixed
            CALL STORE_RHOLM( RHOLM, RHOLM_STORE(IBASE:,ISP),  &
                    METRIC(IBASE:), IADD, PP, .FALSE., ITMP )
            CALL TRANS_RHOLM_IM( CRHODE(:,:,NIP,ISP), RHOLM, PP )
            CALL STORE_RHOLM( RHOLM, RHOLM_STORE(IBASE+ISIZE/2:,ISP),  &
                    METRIC(IBASE+ISIZE/2:), IADD, PP, .FALSE., ITMP )
         ENDDO
         IBASE=IBASE+IADD
         IF (IBASE > ISIZE+1) THEN
            WRITE(0,*) 'internal error SET_RHO_PAW: insufficient space'
            STOP
         ENDIF
      ENDDO ion
      ! WRITE(*,'("RHOOUT",6F10.6)') RHOLM_STORE

    END SUBROUTINE SET_RHO_PAW

!*******************************************************************
!
! WRT_RHO_PAW write the PAW occupancies to a file specified by
! IU
!
!*******************************************************************

    SUBROUTINE WRT_RHO_PAW(P, T_INFO, LOVERL, RHOLM_STORE, COMM, IU )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      IMPLICIT NONE

      INTEGER IU               ! io unit
      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      LOGICAL  LOVERL          ! overlap matrix used ?
      REAL(q) RHOLM_STORE(:)   ! storage for the channel occupancies
      TYPE(communic) :: COMM

    ! local variables
      TYPE (potcar),POINTER :: PP
      INTEGER NT, NI, I
      INTEGER IBASE, IADD, NELEMENTS, LYMAX, LMMAX
      INTEGER, EXTERNAL :: MAXL_AUG
      INTEGER NODE_ME, IONODE

      REAL(q), ALLOCATABLE::  BUFFER(:)


!=======================================================================
! quick return if possible
!=======================================================================
      IF (.NOT.LOVERL .OR. MIMIC_US ) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      ALLOCATE( BUFFER(LMMAX*LMMAX))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
      IBASE=1

      ion: DO NI=1,T_INFO%NIONS
         BUFFER=0
         IADD  =0

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)

         NELEMENTS =0 
         IF (DO_LOCAL(NI)) THEN
            CALL RETRIEVE_RHOLM( BUFFER, RHOLM_STORE(IBASE:), &
                 METRIC(IBASE:), IADD, PP, .TRUE., NELEMENTS)

            IF (NELEMENTS > LMMAX*LMMAX) THEN
               WRITE(*,*)'internal ERROR: WRT_RHO_PAW running out of buffer'
               STOP
            ENDIF
            
            IBASE=IBASE+IADD
         ENDIF

         
         
         
         
         WRITE(IU,'("augmentation occupancies",2I4)') NI, NELEMENTS
         WRITE(IU,'(5E15.7)') (BUFFER(I),I=1,NELEMENTS)
         

      ENDDO ion

      DEALLOCATE(BUFFER)
    END SUBROUTINE WRT_RHO_PAW

!*******************************************************************
!
! RD_RHO_PAW read the PAW occupancies from a file specified by
! IU
! only components up to LMAXPAW are read in
! the remaining components are required by the routine
! SET_DD_PAW
! the currentl logic is that these are determined from CRHODE
! which is initialised according to atomic occupancies
!
!*******************************************************************

    SUBROUTINE RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM_STORE, COMM, IU, IERR )
      USE pseudo
      USE poscar
      USE wave
      USE constant
      IMPLICIT NONE

      INTEGER IU               ! io unit
      INTEGER IERR             ! error status on return
                               ! 0 ok, 1 error occured
      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      LOGICAL  LOVERL          ! overlap matrix used ?
      REAL(q) RHOLM_STORE(:)   ! storage for the channel occupancies
      TYPE(communic) :: COMM

    ! local variables
      TYPE (potcar),POINTER :: PP
      INTEGER NT, NI, I
      INTEGER IBASE, IADD, NELEMENTS, LYMAX, LMMAX
      INTEGER, EXTERNAL :: MAXL_AUG
      INTEGER NI_READ, NELEMENTS_READ
      INTEGER NODE_ME, IONODE
      REAL(q) :: RELEMENTS

      REAL(q), ALLOCATABLE::  BUFFER(:)
      CHARACTER*(80) CH


!=======================================================================
! quick return if possible
!=======================================================================
      IF (.NOT.LOVERL .OR. MIMIC_US ) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      ALLOCATE( BUFFER(LMMAX*LMMAX))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
      IBASE=1

      ion: DO NI=1,T_INFO%NIONS
         BUFFER=0
         IADD  =0

         NT=T_INFO%ITYP(NI)
         PP=> P(NT)

         IERR=0


         
         NI_READ=0
         NELEMENTS_READ=0

         READ(IU,'(24X,2I4)',IOSTAT=IERR) &
              NI_READ, NELEMENTS_READ

         IF (NELEMENTS_READ > LMMAX*LMMAX) THEN
            WRITE(*,*)'internal ERROR: RD_RHO_PAW running out of buffer'
            STOP
         ENDIF

         IF (IERR == 0 ) READ(IU,*,IOSTAT=IERR) (BUFFER(I),I=1,NELEMENTS_READ)
         IF (NI_READ /= NI .OR. IERR/=0 ) THEN
            IERR=1
         ENDIF
         

         
         

         NELEMENTS=0
         IF (DO_LOCAL(NI)) THEN
            CALL STORE_RHOLM( BUFFER, RHOLM_STORE(IBASE:), &
                 METRIC(IBASE:), IADD, PP, .TRUE. , NELEMENTS)
            IBASE=IBASE+IADD
         ENDIF

         
         

         
         IF (NELEMENTS /= NELEMENTS_READ)  THEN
           IERR=1
         ENDIF
         

         
         
         IF (IERR/=0) THEN
             WRITE(*,*) 'RD_RHO_PAW: ion', NI,'data corrupt'
            DEALLOCATE (BUFFER)
            RETURN
         ENDIF
         
      ENDDO ion
      DEALLOCATE(BUFFER)

    END SUBROUTINE RD_RHO_PAW

!***********************************************************************
!                                                                      *
!   Routine AUGSYM symmetrizes arrays like RHOLM or DLM                *
!   or other arrays with a similar structure                           *
!   The procedure is a quite                                           *
!   simple: Rotate the input array and add it to the array at the      *
!   rotated atomic position (for all space group operations ...).      *
!   Finally divide the result by the number of symmetry operations!    *
!                                                                      *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!   see routine FSYM in symlib.F                                       *
!                                                                      *
!      MAP(NAX,1,48,NP) contains a table connecting the rotated and    *
!                   the unrotated positions (stored is the index of    *
!                   the rotated position)                              *
!      S(3,3,48) contains the INTEGER rotation matrices.               *
!      NR contains the number of given space group operations.         *
!      NP contains the number of "primitive" translations in the cell. *
!                                                                      *
!      NTYP is the number of atomic species.                           *
!      NITYP(NSP) contains the number of atoms per species.            *
!                                                                      *
!      A(:,I)  contains the direct (real space) lattice vectors.       *
!      B(:,I)  contains the reciprocal lattice vectors.                *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      Array MAT contains the symmetrized matrix on output (input is   *
!      overwritten!).                                                  *
!                                                                      *
!                                                                      *
!***********************************************************************


      SUBROUTINE AUGSYM(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
            NIOND,NR,NP,MAP,MAGROT,S,A,B,ISP)
      USE prec
      USE pseudo
      IMPLICIT NONE

! parameters required from the symmetry package
      INTEGER NIOND      ! dimension for number of ions
      INTEGER NIONS      ! number of ions
      INTEGER NP         ! number of primitive translations
      INTEGER MAP(NIOND,48,NP)
      INTEGER NR         ! number of rotations
      INTEGER NTYP       ! number of species
      INTEGER NITYP(NTYP)! number of atoms per species
      INTEGER ISP        ! spin component
      INTEGER S(3,3,48)
      REAL(q) A(3,3),B(3,3),MAGROT(48,NP)
!
      INTEGER LMDIM      ! first dimension of MAT
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS)
      TYPE (potcar) P(NTYP)

! local varibales
      COMPLEX(q),ALLOCATABLE :: TMP(:,:,:),ROTMAT(:,:,:)
      REAL(q),ALLOCATABLE :: SL(:,:,:)
      INTEGER NROT,LMAX,LDIM,ITRANS,IA,IAP,MMAX,IROT,IS,ISTART
      INTEGER,EXTERNAL :: MAXL,MAXL1
      REAL(q) SCALE
      INTEGER L,LP,NI
!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------



      LDIM=MAXL(NTYP,P  )         ! maximum l quantum number
      MMAX=(2*LDIM+1)             ! maximum m quantum number
      ALLOCATE (TMP(LMDIM,LMDIM,NIONS),SL(MMAX,MMAX,0:LDIM), &
                ROTMAT(LMDIM,LMDIM,NIONS))

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        LMAX=MAXL1(P(IS))      ! maximum l for this species
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NR
          ! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LMAX,S(1,1,IROT),SL,A,B)
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! rotate the matrix and store result in ROTMAT
            CALL ROTATE_MATRIX(LMDIM,MAT(1,1,IA),ROTMAT(1,1,IA),MMAX,LMAX,SL,P(IS))
          ENDDO

          ! loop over all space group operations (translations+ rotations)
          DO IA=ISTART,NITYP(IS)+ISTART-1
            DO ITRANS=1,NP
               ! destination atom
               IAP=MAP(IA,IROT,ITRANS)
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT,ITRANS)
               TMP(:,:,IA)=TMP(:,:,IA)+ROTMAT(:,:,IAP)*SCALE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NP*NR)
      MAT=TMP*SCALE

      DEALLOCATE(TMP,ROTMAT,SL)

      END SUBROUTINE
!************** SUBROUTINE AUGSYM_NONCOL *******************************
!                                                                      *
!   Routine AUGSYM_NONCOL symmetrizes arrays like RHOLM or DLM         *
!   or other arrays with a similar structure                           *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!   see routine FSYM in symlib.F                                       *
!                                                                      *
!      MAP(NAX,1,48,NP) contains a table connecting the rotated and    *
!                   the unrotated positions (stored is the index of    *
!                   the rotated position)                              *
!      S(3,3,48) contains the INTEGER rotation matrices.               *
!      NR contains the number of given space group operations.         *
!      NP contains the number of "primitive" translations in the cell. *
!                                                                      *
!      NTYP is the number of atomic species.                           *
!      NITYP(NSP) contains the number of atoms per species.            *
!                                                                      *
!      A(:,I)  contains the direct (real space) lattice vectors.       *
!      B(:,I)  contains the reciprocal lattice vectors.                *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      Array MAT contains the symmetrized matrix on output (input is   *
!      overwritten!).                                                  *
!                                                                      *
!                                                                      *
!***********************************************************************

      SUBROUTINE AUGSYM_NONCOL(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
     &      NIOND,NR,NP,MAP,MAGROT,SAXIS,S,INVMAP,A,B)
      USE prec
      USE pseudo
      USE relativistic
      
      IMPLICIT NONE

! parameters required from the symmetry package
      INTEGER NIOND      ! dimension for number of ions
      INTEGER NIONS      ! number of ions
      INTEGER NP         ! number of primitive translations
      INTEGER MAP(NIOND,48,NIOND)
      INTEGER NR         ! number of rotations
      INTEGER NTYP       ! number of species
      INTEGER NITYP(NTYP)! number of atoms per species
      INTEGER ISP        ! spin component
      INTEGER S(3,3,48),INVMAP(48),I
      REAL(q) A(3,3),B(3,3),MAGROT(48,NP)
!
      INTEGER LMDIM      ! first dimension of MAT
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS,3)
      TYPE (potcar) P(NTYP)

! local varibales
      COMPLEX(q),ALLOCATABLE :: TMP(:,:,:,:),ROTMAT(:,:,:,:)
      COMPLEX(q),ALLOCATABLE :: ROTMAT_TEMP(:,:,:),TROTMAT(:,:,:)
      REAL(q),ALLOCATABLE :: SL(:,:,:)
      INTEGER NROT,LMAX,LDIM,ITRANS,IA,IAP,MMAX,IROT,IROTI,IS,ISTART,IDIR
      INTEGER,EXTERNAL :: MAXL,MAXL1
      REAL(q) SCALE,SAXIS(3),ALPHA,BETA
      INTEGER L,LP,NI
!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------

      LDIM=MAXL(NTYP,P  )         ! maximum l quantum number
      MMAX=(2*LDIM+1)             ! maximum m quantum number
      ALLOCATE (TMP(LMDIM,LMDIM,NIONS,3),SL(MMAX,MMAX,0:LDIM), &
                ROTMAT(LMDIM,LMDIM,NIONS,3),ROTMAT_TEMP(LMDIM,LMDIM,3), &
		TROTMAT(LMDIM,LMDIM,3))

      CALL EULER(SAXIS,ALPHA,BETA)

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        LMAX=MAXL1(P(IS))      ! maximum l for this species
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NR
          ! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LMAX,S(1,1,IROT),SL,A,B)
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! rotate the matrix and store result in ROTMAT
	    DO IDIR=1,3
            CALL ROTATE_MATRIX(LMDIM,MAT(1,1,IA,IDIR),ROTMAT(1,1,IA,IDIR),MMAX,LMAX,SL,P(IS))
	    ENDDO
          ENDDO
          ! loop over all space group operations (translations+ rotations)
          DO IA=ISTART,NITYP(IS)+ISTART-1
            DO ITRANS=1,NP
               ! destination atom
               IAP=MAP(IA,IROT,ITRANS)
               SCALE=1._q
! Transform from "SAXIS basis" to the system of cartesian axes
! in which the integer rotation matrices are defined
               ROTMAT_TEMP(:,:,1)=COS(BETA)*COS(ALPHA)*ROTMAT(:,:,IAP,1)- &
                    SIN(ALPHA)*ROTMAT(:,:,IAP,2)+ &
                    SIN(BETA)*COS(ALPHA)*ROTMAT(:,:,IAP,3)
               ROTMAT_TEMP(:,:,2)=COS(BETA)*SIN(ALPHA)*ROTMAT(:,:,IAP,1)+ &
                    COS(ALPHA)*ROTMAT(:,:,IAP,2)+ &
                    SIN(BETA)*SIN(ALPHA)*ROTMAT(:,:,IAP,3)
               ROTMAT_TEMP(:,:,3)=-SIN(BETA)*ROTMAT(:,:,IAP,1)+ &
                    COS(BETA)*ROTMAT(:,:,IAP,3)
! Bring to direct coordinates
               CALL MAT_KARDIR(LMDIM,ROTMAT_TEMP,B)
! Rotate in direct space
               TROTMAT(:,:,1)=ROTMAT_TEMP(:,:,1)*S(1,1,IROT)+ &
                    ROTMAT_TEMP(:,:,2)*S(2,1,IROT)+ &
                    ROTMAT_TEMP(:,:,3)*S(3,1,IROT)
               
               TROTMAT(:,:,2)=ROTMAT_TEMP(:,:,1)*S(1,2,IROT)+ &
                    ROTMAT_TEMP(:,:,2)*S(2,2,IROT)+ &
                    ROTMAT_TEMP(:,:,3)*S(3,2,IROT)
	       
               TROTMAT(:,:,3)=ROTMAT_TEMP(:,:,1)*S(1,3,IROT)+ &
                    ROTMAT_TEMP(:,:,2)*S(2,3,IROT)+ &
                    ROTMAT_TEMP(:,:,3)*S(3,3,IROT)
! bring TROTMAT  to cartesian coordinates
               CALL MAT_DIRKAR(LMDIM,TROTMAT,A)
! And back to SAXIS representation
               ROTMAT_TEMP(:,:,1)=COS(BETA)*COS(ALPHA)*TROTMAT(:,:,1)+ &
                    COS(BETA)*SIN(ALPHA)*TROTMAT(:,:,2)- &
                    SIN(BETA)*TROTMAT(:,:,3)
               ROTMAT_TEMP(:,:,2)=-SIN(ALPHA)*TROTMAT(:,:,1)+ &
                    COS(ALPHA)*TROTMAT(:,:,2)
               ROTMAT_TEMP(:,:,3)=SIN(BETA)*COS(ALPHA)*TROTMAT(:,:,1)+ &
                    SIN(BETA)*SIN(ALPHA)*TROTMAT(:,:,2)+ &
                    COS(BETA)*TROTMAT(:,:,3)
! And then we sum
               TMP(:,:,IA,:)=TMP(:,:,IA,:)+ROTMAT_TEMP(:,:,:)*SCALE

            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NP*NR)
      MAT=TMP*SCALE

      DEALLOCATE(TMP,ROTMAT,SL,ROTMAT_TEMP,TROTMAT)

      END SUBROUTINE


!**************** SUBROUTINE MAT_KARDIR ********************************
! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE MAT_KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(NMAX,NMAX,3),V1,V2,V3
      DIMENSION  BASIS(3,3)

      DO N=1,NMAX
      DO M=1,NMAX
        V1=V(N,M,1)*BASIS(1,1)+V(N,M,2)*BASIS(2,1)+V(N,M,3)*BASIS(3,1)
        V2=V(N,M,1)*BASIS(1,2)+V(N,M,2)*BASIS(2,2)+V(N,M,3)*BASIS(3,2)
        V3=V(N,M,1)*BASIS(1,3)+V(N,M,2)*BASIS(2,3)+V(N,M,3)*BASIS(3,3)
        V(N,M,1)=V1
        V(N,M,2)=V2
        V(N,M,3)=V3
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE MAT_DIRKAR ********************************
! transform a set of vectors from
! ) direct lattice      (BASIS must be equal to A direct lattice)
! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
! to cartesian coordinates
!***********************************************************************

      SUBROUTINE MAT_DIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(NMAX,NMAX,3),V1,V2,V3
      DIMENSION  BASIS(3,3)
      
      DO N=1,NMAX
      DO M=1,NMAX
        V1=V(N,M,1)*BASIS(1,1)+V(N,M,2)*BASIS(1,2)+V(N,M,3)*BASIS(1,3)
        V2=V(N,M,1)*BASIS(2,1)+V(N,M,2)*BASIS(2,2)+V(N,M,3)*BASIS(2,3)
        V3=V(N,M,1)*BASIS(3,1)+V(N,M,2)*BASIS(3,2)+V(N,M,3)*BASIS(3,3)
        V(N,M,1)=V1
        V(N,M,2)=V2
        V(N,M,3)=V3
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!-------------------------------------------------------------------
! subroutine to rotate (1._q,0._q) matrix MAT and store result in ROTMAT
!-------------------------------------------------------------------

      SUBROUTINE ROTATE_MATRIX(LMDIM,MAT,ROTMAT,MMAX,LMAX,SL,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER LMDIM,MMAX,LMAX

      COMPLEX(q) :: MAT(LMDIM,LMDIM)    ! initial matrix
      COMPLEX(q) :: ROTMAT(LMDIM,LMDIM) ! final matrix
      REAL(q) :: SL(MMAX,MMAX,0:LMAX)     ! rotation matrix (allways symmetric)
      TYPE (potcar) P
! local variables
      INTEGER CHANNEL,CHANNELS,IND,L,M,MP
      REAL(q) :: TMP(LMDIM,LMDIM)

      CHANNELS=P%LMAX
! left hand transformation
      IND=0
      TMP=0

      DO CHANNEL=1,CHANNELS
        ! l-qantum number of this channel
        L=P%LPS(CHANNEL)
        ! rotate this l-block
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          TMP(IND+M,:)=TMP(IND+M,:)+SL(M,MP,L)*MAT(IND+MP,:)
        ENDDO
        ENDDO

        IND=IND+(2*L+1)
      ENDDO
! right hand transformation
      IND=0
      ROTMAT=0

      DO CHANNEL=1,CHANNELS
        ! l-qantum number of this channel
        L=P%LPS(CHANNEL)
        ! rotate this l-block
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          ROTMAT(:,IND+M)=ROTMAT(:,IND+M)+SL(M,MP,L)*TMP(:,IND+MP)
        ENDDO
        ENDDO

        IND=IND+(2*L+1)
      ENDDO
      END SUBROUTINE

!*******************************************************************
!
! this subroutine builds up the transformation matrix for
! vectors transforming according to the quantum numbers
! L=0,...,LMAX from the integer transformation matrix S
!   Input parameters:
!   -----------------
!      S(3,3,48) the INTEGER rotation matrices.
!      A(:,I)    the direct (real space) lattice vectors.
!      B(:,I)    the reciprocal lattice vectors.
!
!   Output parameters:
!   ------------------
!      Array SL contains the transformation matrices
!
!*******************************************************************

      SUBROUTINE SETUP_SYM_LL(MMAX,LMAX,S,SL,A,B)
      USE prec
      USE asa
      IMPLICIT NONE

      INTEGER LMAX   !  maximum L quantum number
      INTEGER MMAX   !  first and second dimension of array SL
      INTEGER S(3,3)
      REAL(q) S_(3,3)
      REAL(q) SL(MMAX,MMAX,0:LMAX)
      REAL(q) A(3,3),B(3,3)
      INTEGER L,LP,LM,LNEW,LMINDX,ISTART,IEND,M,MP,M2,MP2,LMINDX2, &
              ISTART2,IEND2,LM2,IC,IC2,LP2,LSET,LMINDX0
      INTEGER I,J,K
      SL=0
!-----------------------------------------------------------------------
! transformation matrix for L=0 (scalar) is simply 1
!-----------------------------------------------------------------------
      SL(1,1,0) = 1
!-----------------------------------------------------------------------
! L=1 transforms like a vector (y,z,x) (see SETYLM in asa.F)
! build first the rotation matrix for a cartesian vector a
! the indexing should be
! a_trans(l) = \sum_i SL(l,i) a(i)
!-----------------------------------------------------------------------
      S_=0
! &*#%^!(#& FUCK THIS CRAY COMPILER AT HIGHER OPTIMIZATION LEVELS (*(@)*&!
! The compiler is really able to "optimize" these four loops in such a way
! that S_ contains either zeros or NaNs but no reasonable result at all ...
!DIR$ NOPATTERN
      DO L=1,3
      DO K=1,3
      DO J=1,3
      DO I=1,3
        S_(L,I)=S_(L,I)+A(L,K)*S(J,K)*B(I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ! shift components
      SL(1,1,1) = S_(2,2) ; SL(1,2,1) = S_(2,3) ; SL(1,3,1) = S_(2,1)
      SL(2,1,1) = S_(3,2) ; SL(2,2,1) = S_(3,3) ; SL(2,3,1) = S_(3,1)
      SL(3,1,1) = S_(1,2) ; SL(3,2,1) = S_(1,3) ; SL(3,3,1) = S_(1,1)
!      LNEW=1
!      WRITE(0,*) LNEW
!      DO M=1,2*LNEW+1
!         WRITE(0,'(10F10.6)') (SL(M,M2,LNEW),M2=1,2*LNEW+1)
!      ENDDO
!-----------------------------------------------------------------------
! for all remaining L the transformation matrix is build up from L=1
! using Clebsch-Gordan like coefficients
!-----------------------------------------------------------------------
      LSET=1

      LP=1
      DO L=LSET,LMAX-1
         CALL YLM3LOOKUP(L,LP,LMINDX0)
         LNEW=L+LP

         LMINDX=LMINDX0
         DO M = 1, 2*L +1
         DO MP= 1, 2*LP+1
            LMINDX=LMINDX+1

            ISTART=INDCG(LMINDX) ; IEND  =INDCG(LMINDX+1)

            DO IC=ISTART,IEND-1
               LM=JS(IC)
               IF (LM > LNEW*LNEW       .AND. &
                   LM <= (LNEW+1)*(LNEW+1)) THEN
                 LM=LM-LNEW*LNEW

                 LMINDX2=LMINDX0
                 DO M2 = 1, 2*L +1
                 DO MP2= 1, 2*LP+1
                    LMINDX2=LMINDX2+1
                    ISTART2=INDCG(LMINDX2) ; IEND2  =INDCG(LMINDX2+1)
                    DO IC2=ISTART2,IEND2-1
                       LM2=JS(IC2)
                       IF (LM2 > LNEW*LNEW       .AND. &
                           LM2 <= (LNEW+1)*(LNEW+1)) THEN
                           LM2=LM2-LNEW*LNEW
                           ! order of elements checked by recalculating L=1 term
                           SL(LM,LM2,LNEW)= SL(LM,LM2,LNEW)+ &
                           YLM3(IC)*SL(M,M2,L)*SL(MP,MP2,LP)*YLM3I(IC2)
                       ENDIF
                    ENDDO
                 ENDDO
                 ENDDO

               ENDIF
            ENDDO
         ENDDO
         ENDDO
!         WRITE(0,*) LNEW
!         DO M=1,2*LNEW+1
!         WRITE(0,'(10F10.6)') (SL(M,M2,LNEW),M2=1,2*LNEW+1)
!         ENDDO
      ENDDO

      END SUBROUTINE

!*******************************************************************
!
! this subroutine calculates the local electronic potential
! for atomic occupancies within the augmentation spheres;
! it can be (and is) called to reduce numerice errors introduced by
! different implementations of the exchange correlation funtional
! used in the pseudopotential generation code and VASP
!
!*******************************************************************

    SUBROUTINE SET_ATOM_POT( P , T_INFO, LOVERL, LMDIM, LEXCH, LEXCHG )

      USE pseudo
      USE asa
      USE poscar
      USE wave
      USE constant
      USE radial
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar),TARGET::      P(T_INFO%NTYP)
      INTEGER LMDIM,LEXCH, LEXCHG
      LOGICAL  LOVERL

    ! local variables
      TYPE (potcar),POINTER:: PP
      INTEGER NT,LYMAX,NI,NDIM,LMMAX,LOW,LHI,L,LP,LL,LM,M,MMAX,K
      INTEGER,PARAMETER :: ISPIN=1
      INTEGER, EXTERNAL :: MAXL_AUG,MAXL1
      REAL(q)  RHOLM(LMDIM*LMDIM)
      COMPLEX(q)  CRHODE(LMDIM,LMDIM)
      REAL(q),ALLOCATABLE :: POT(:,:,:), RHO(:,:,:),POTOLD(:)
      REAL(q) :: DOUBLEC_AE,DOUBLEC_PS,EXCG
      REAL(q) :: DOUBLEPS,DOUBLEAE,SCALE
      REAL(q) :: Z
    ! variables required to store core wavefunctions
      INTEGER MAXNL
      REAL(q), ALLOCATABLE :: W(:,:), EIG(:)
      INTEGER, ALLOCATABLE :: N(:), LC(:)
!=======================================================================
! quick return and allocation of work space
!=======================================================================
      DOUBLEC_AE=0
      DOUBLEC_PS=0
      IF (.NOT.LOVERL) RETURN

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2

      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO

      IF (NDIM == 0) RETURN
      ALLOCATE ( POT( NDIM, LMMAX, ISPIN ), RHO( NDIM, LMMAX, ISPIN ),POTOLD(NDIM))

!=======================================================================
! cycle all ions and add corrections to pseudopotential strength
!=======================================================================
      SCALE=2*SQRT(PI)

      DO NT=1,T_INFO%NTYP
         PP=> P(NT)
         IF ( .NOT. ASSOCIATED(P(NT)%QPAW )) CYCLE
         ! set the atomic occupancies
         CRHODE=0
         LYMAX =MAXL1(PP)*2

         LOW=1
         LM =1
         block: DO
           LL=PP%LPS(LOW)
         ! search block with same L
           DO LHI=LOW,PP%LMAX
              IF (LL/=PP%LPS(LHI)) EXIT
           ENDDO
           LHI=LHI-1
           MMAX=2*LL+1
  !-----------------------------------------------------------------------
  ! first set RHOLM (i.e. the on site occupancy matrix)
  !-----------------------------------------------------------------------
           DO L =LOW,LHI
           DO LP=LOW,LHI
           DO M =0,MMAX-1
              CRHODE(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M)=PP%QATO(L,LP)
           ENDDO; ENDDO; ENDDO

         ! set new LOW value and LM value and go on
            LM=LM+(LHI-LOW+1)*MMAX
            LOW=LHI+1
            IF (LOW > PP%LMAX) EXIT block
         ENDDO block
         ! transform CRHODE to llp,LM
         RHOLM=0
         CALL TRANS_RHOLM( CRHODE, RHOLM, PP )

  !-----------------------------------------------------------------------
  ! calculate the local radial potential
  !-----------------------------------------------------------------------
         POTOLD=PP%POTAE
         PP%POTAE=0
         RHO=0
         CALL RAD_CHARGE( RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WAE )
         CALL RAD_POT( PP%R, ISPIN, LEXCH, LEXCHG, 1, 1,  &
                 RHO, PP%RHOAE, PP%POTAE,  POT, DOUBLEAE, EXCG)

         CALL SIMPI(PP%R, PP%RHOAE, Z)
   ! calculate  total electronic core charge (well in future we should include this in POTCAR
   ! I think)
   ! we round up in all cases, since some electronic charge might be outside the
   ! augmentation sphere
         PP%ZCORE= AINT(Z*sqrt(4*PI)+0.9)

         DO K=1,P(NT)%R%NMAX
           PP%POTAE(K)=-POT(K,1,1)/SCALE
!           WRITE(77,'(I4,5F14.7)') K,PP%R%R(K)/AUTOA,RHO(K,1,1),PP%POTAE(K),PP%POTAE(K)-POTOLD(K)
         ENDDO
  ! calculate the local radial PS potential
  ! mind this is the valence only contribution
         POTOLD=PP%POTPS
         PP%POTPS=0
         RHO=0
         CALL RAD_CHARGE( RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WPS )
         CALL RAD_AUG_CHARGE(  RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
              LYMAX, PP%AUG, PP%QPAW )
         CALL RAD_POT( PP%R, ISPIN, LEXCH, LEXCHG, 1, 1,  &
                 RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCG)

         DO K=1,P(NT)%R%NMAX
           PP%POTPS(K)=-POT(K,1,1)/SCALE
!           WRITE(77,'(I4,5F14.7)') K,PP%R%R(K)/AUTOA,RHO(K,1,1)*AUTOA*SCALE,PP%POTPS(K),PP%POTPS(K)-POTOLD(K)
         ENDDO

         DOUBLEC_PS= DOUBLEC_PS-DOUBLEPS*T_INFO%NITYP(NT)
         DOUBLEC_AE= DOUBLEC_AE+(DOUBLEAE-PP%DEXCCORE)*T_INFO%NITYP(NT)

      ENDDO

      DOUBLEC_PS_ATOM=DOUBLEC_PS
      DOUBLEC_AE_ATOM=DOUBLEC_AE

      DEALLOCATE( POT, RHO , POTOLD )
    END SUBROUTINE SET_ATOM_POT

!*******************************************************************
!
! this is the central routine of the PAW method
! it calculates the on site terms and adds them to D(I,J)
! at the same time double counting corrections are calculated
!
!*******************************************************************

    SUBROUTINE SET_DD_PAW(WDES, P , T_INFO, LOVERL, &
         ISPIN, LMDIM, CDIJ, RHOLM_STORE, CRHODE, LEXCH, LEXCHG, &
         E, LMETA, LASPH, LCOREL )
      USE pseudo
      USE asa
      USE poscar
      USE wave
      USE constant
      USE radial
      USE base
      USE relativistic
      USE LDAPLUSU_MODULE
      USE cl
      IMPLICIT NONE

      TYPE (type_info) T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      TYPE (wavedes)  WDES
      TYPE (energy)   E
      INTEGER LMDIM, ISPIN
      INTEGER LEXCH, LEXCHG
      COMPLEX(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q)  RHOLM_STORE(:,:)
      LOGICAL  LOVERL
      LOGICAL  LMETA      ! calculate meta GGA contribution
      LOGICAL  LASPH      ! calculate aspherical corrections to potential
      LOGICAL  LCOREL     ! calculate accurate core level shifts
    ! local variables
      TYPE (potcar),POINTER:: PP
      INTEGER NT,LYMAX,NI,NDIM,LMMAX,NIP,ISP,IBASE,IADD,ISIZE,K, ITMP,NCDIJ,LMAX
      INTEGER, EXTERNAL :: MAXL_AUG,MAXL1
      REAL(q) DDLM(LMDIM*LMDIM),RHOLM(LMDIM*LMDIM),RHOLM_(LMDIM*LMDIM,WDES%NCDIJ)
!     COMPLEX(q) CTMP(LMDIM,LMDIM,WDES%NCDIJ),CSO(LMDIM,LMDIM,WDES%NCDIJ)
!     COMPLEX(q) COCC(LMDIM,LMDIM,WDES%NCDIJ)
      COMPLEX(q) CTMP(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),CSO(LMDIM,LMDIM,WDES%NCDIJ)
      COMPLEX(q) COCC(LMDIM,LMDIM,MAX(2,WDES%NCDIJ)),COCC_IM(LMDIM,LMDIM)
      
      REAL(q),ALLOCATABLE :: POT(:,:,:), RHO(:,:,:), POTAE(:,:,:), RHOAE(:,:,:)
      REAL(q),ALLOCATABLE :: RHOCOL(:,:,:)
      ! core level shifts
      REAL(q),ALLOCATABLE :: DRHOCORE(:)
      REAL(q) :: DOUBLEC_AE,DOUBLEC_PS
      REAL(q) :: DOUBLEPS,DOUBLEAE
      REAL(q) :: EXCM,EXCG,EXCGA
!OBengone modify start
      REAL(q) :: DOUBLEC_LDAU
!OBengone modify end
      ! euler angles of the global spin quantisation axis
      REAL(q) :: ALPHA,BETA
! #define robdeb
    ! kinetic energy density (true and Weizsaecker)
      REAL(q),ALLOCATABLE :: KINDENSAE(:,:),KINDENSPS(:,:)
      REAL(q),ALLOCATABLE :: WKDAE(:,:),WKDPS(:,:)
      REAL(q),ALLOCATABLE :: RHOUPD(:,:,:),RHOLMUPD(:,:)
      REAL(q),POINTER :: NULPOINTER(:)
      REAL(q) :: SPI2,TMP
      INTEGER II,RNMAX
! variables required to store core wavefunctions
      INTEGER MAXNL
      REAL(q), ALLOCATABLE :: W(:,:), EIG(:)
      INTEGER, ALLOCATABLE :: N(:), LC(:)
!=======================================================================
! quick return and allocation of work space
!=======================================================================

      DOUBLEC_AE=0
      DOUBLEC_PS=0
      E%PAWPSM=0; E%PAWAEM=0
      E%PAWPSG=0; E%PAWAEG=0
      E%PAWCORE=0
      E%PAWCOREM=0
      E%PAWPSAS=0; E%PAWAEAS=0

      CL_SHIFT= 0

      NULLIFY(NULPOINTER)

      IF (.NOT.LOVERL) RETURN
! mimic US-PP just set the double counting corrections correctly
      IF (MIMIC_US) THEN
         DOUBLEC_AE=DOUBLEC_AE_ATOM
         DOUBLEC_PS=DOUBLEC_PS_ATOM
         E%PAWAE=DOUBLEC_AE
         E%PAWPS=DOUBLEC_PS
         RETURN
      ENDIF

      SPI2= 2*SQRT(PI)

      LYMAX =MAXL_AUG(T_INFO%NTYP,P)

      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO

      IF (NDIM == 0) RETURN

      LMMAX=(LYMAX+1)**2
      NCDIJ = WDES%NCDIJ
      ALLOCATE ( POT( NDIM, LMMAX, NCDIJ ), RHO( NDIM, LMMAX, NCDIJ ), &
           POTAE( NDIM, LMMAX, NCDIJ ), RHOAE( NDIM, LMMAX, NCDIJ), DRHOCORE(NDIM))

      ALLOCATE (RHOCOL( NDIM, LMMAX, NCDIJ ))
! allocate kinetic energy density if metagga
      IF (LMETA) THEN
         ALLOCATE (KINDENSAE(NDIM,NCDIJ),KINDENSPS(NDIM,NCDIJ)) 
         ALLOCATE (WKDAE(NDIM,NCDIJ),WKDPS(NDIM,NCDIJ))
         ALLOCATE (RHOUPD(NDIM,1,NCDIJ),RHOLMUPD(LMDIM*LMDIM,NCDIJ))
      ENDIF
      
! for spin orbit coupling set the euler angles
      IF ( WDES%LSORBIT ) &
         CALL EULER(WDES%SAXIS, ALPHA, BETA)
!=======================================================================
! cycle all ions and add corrections to pseudopotential strength CDIJ
!=======================================================================
      IBASE=1

      ion: DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB) ! local storage index
         NT=T_INFO%ITYP(NI)
      ! if this element is not treated locally CYCLE
         IF (.NOT. DO_LOCAL(NI)) THEN
            ! for PAW, set CDIJ to (0._q,0._q) if it resides on local node
            ! and if the element is not treated locally
            IF (ASSOCIATED(P(NT)%QPAW)) THEN
               IF (NIP /= 0) THEN
                  DO ISP=1,NCDIJ
                     CDIJ(:,:,NIP,ISP)=0
                  ENDDO
               ENDIF
            ELSE
            ! US PP: initialize to (0._q,0._q) if we are not on first node in COMM_INTER
            ! (at the end, we use a global sum over COMM_INTER) 
            ENDIF
            CYCLE ion
         ENDIF

         PP=> P(NT)
         LYMAX =MAXL1(PP)*2
         RNMAX =PP%R%NMAX
  !-----------------------------------------------------------------------
  ! first set RHOLM (i.e. the on site occupancy matrix)
  ! and then the lm dependent charge densities RHO and RHOAE
  ! (excluding augmentation charges yet)
  !-----------------------------------------------------------------------
         RHOLM_=0

         IF ( LMAX_MIX < PP%LMAX_CALC) THEN
            DO ISP=1,NCDIJ
               ! transform CRHODE to llp,LM
               CALL TRANS_RHOLM( CRHODE(:,:,NIP,ISP), RHOLM_(:,ISP), PP )
            ENDDO
         ENDIF

  !      WRITE(*,'("RHOMIX",6F10.6)') RHOLM_STORE
         COCC=0
         ISIZE=UBOUND(RHOLM_STORE,1)

         DO ISP=1,NCDIJ
            RHOLM=RHOLM_(:,ISP)
    ! retrieve mixed elements from RHOLM_STORE and overwrite them in RHOLM
            CALL RETRIEVE_RHOLM( RHOLM, RHOLM_STORE(IBASE:,ISP), &
                           METRIC(IBASE:), IADD, PP, .FALSE.,  ITMP)
    !       WRITE(*,'(I4,(10F10.3))') IBASE,RHOLM(1:IADD)
    ! calculate the total radial charge
            IF (LMETA) RHOLMUPD(:,ISP)=RHOLM(:)
            LMMAX=(LYMAX+1)**2
            RHOAE(:,1:LMMAX,ISP)=0; RHO(:,1:LMMAX,ISP)=0
            CALL RAD_CHARGE( RHOAE(:,:,ISP), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WAE )
            CALL RAD_CHARGE( RHO(:,:,ISP), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WPS )

    ! for LDA+U we need the mixed occupancy matrix
            IF (USELDApU()) THEN            
               CALL TRANS_RHOLMI( COCC(:,:,ISP), RHOLM, PP )
               COCC_IM=0
               ! retrieve imaginary part of the mixed elements 
               CALL TRANS_RHOLM_IM( CRHODE(:,:,NIP,ISP), RHOLM, PP )
!              CALL STORE_RHOLM( RHOLM, RHOLM_STORE(IBASE+ISIZE/2:,ISP),  &
!                        METRIC(IBASE+ISIZE/2:), IADD, PP, .FALSE., ITMP )
!              CALL TRANS_RHOLM_IM(CRHODE(:,:,NIP,ISP), RHOLM,PP)
               CALL RETRIEVE_RHOLM( RHOLM, RHOLM_STORE(IBASE+ISIZE/2:,ISP), &
                                 METRIC(IBASE+ISIZE/2:), IADD, PP, .FALSE.,  ITMP)
               CALL TRANS_RHOLMI_IM( COCC_IM, RHOLM, PP )
               ! join imaginary and real parts
               COCC(:,:,ISP)=COCC(:,:,ISP)+(0._q,1._q)*COCC_IM(:,:)
!              WRITE(*,*) 'spin',ISP
!              CALL DUMP_DLLMM_IM(COCC(:,:,ISP),PP)
!              CALL DUMP_DLLMM_IM(CRHODE(:,:,NIP,ISP),PP)
!              WRITE(*,*) 'spin',ISP
!              CALL DUMP_DLLMM(COCC(:,:,ISP),PP)
!              CALL DUMP_DLLMM(CRHODE(:,:,NIP,ISP),PP)
            ENDIF
         ENDDO
  !-----------------------------------------------------------------------
  ! calculation and output of radial kinetic energy density
  ! (remember augmentation charges are still excluded)
  !-----------------------------------------------------------------------
         IF ( LMETA ) THEN
            KINDENSAE=0; WKDAE=0
            KINDENSPS=0; WKDPS=0
            IF (ISPIN==2) THEN
               CALL FLIP_RAD(RHOLMUPD,RHOLMUPD,LMDIM*LMDIM)
               CALL FLIP_RAD(RHOAE(:,1,1:2),RHOUPD(:,1,1:2),RNMAX)
            ELSE
               RHOUPD(:,1,1)=RHOAE(:,1,1)
            ENDIF

            DO ISP=1,ISPIN
               CALL RAD_KINETIC_EDENS(KINDENSAE(:,ISP),WKDAE(:,ISP),PP%R,RHOLMUPD(:,ISP), &
                    PP%LMAX,PP%LPS, PP%WAE, PP%RHOAE, PP%TAUAE, RHOUPD(:,:,ISP),ISPIN)
            ENDDO
            IF (ISPIN==2) THEN
               CALL FLIP_RAD(RHO(:,1,1:2),RHOUPD(:,1,1:2),RNMAX)
            ELSE
               RHOUPD(:,1,1)=RHO(:,1,1)
            ENDIF
            DO ISP=1,ISPIN
               CALL RAD_KINETIC_EDENS(KINDENSPS(:,ISP),WKDPS(:,ISP),PP%R,RHOLMUPD(:,ISP), &
                    PP%LMAX,PP%LPS, PP%WPS, NULPOINTER, NULPOINTER, RHOUPD(:,:,ISP),ISPIN)
            ENDDO
         ENDIF
  !-----------------------------------------------------------------------
  ! add augmentation charges now
  !-----------------------------------------------------------------------
         DO ISP=1,NCDIJ
            RHOLM=RHOLM_(:,ISP)
    ! retrieve mixed elements from RHOLM_STORE and overwrite them in RHOLM
            CALL RETRIEVE_RHOLM( RHOLM, RHOLM_STORE(IBASE:,ISP), &
                 METRIC(IBASE:), IADD, PP, .FALSE.,  ITMP)

            CALL RAD_AUG_CHARGE(  RHO(:,:,ISP), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
                  LYMAX, PP%AUG, PP%QPAW )
            CALL RAD_INT( PP%R,  LYMAX, RHO(:,:,ISP), RHOAE(:,:,ISP) )
         ENDDO
         IBASE=IBASE+IADD
  !-----------------------------------------------------------------------
  ! now finish the meta GGA stuff
  !-----------------------------------------------------------------------
            IF ( LMETA ) THEN
   ! output of kinetic energy density

      ! calculate E(xc) for meta-GGA
!               CALL RAD_META_GGA(PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC, &
!                    RHO, PP%RHOPS, PP%POTPS,POT,KINDENSPS,WKDPS,EXCM)
               CALL RAD_META_GGA_ASPH(PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC, &
                    RHO, PP%RHOPS, PP%POTPS,POT,KINDENSPS,WKDPS,EXCM)
               E%PAWPSM=E%PAWPSM-EXCM
!               CALL RAD_META_GGA(PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC, &
!                    RHOAE, PP%RHOAE, PP%POTAE,POTAE,KINDENSAE,WKDAE,EXCM)
               CALL RAD_META_GGA_ASPH(PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC, &
                    RHOAE, PP%RHOAE, PP%POTAE,POTAE,KINDENSAE,WKDAE,EXCM)
               E%PAWAEM=E%PAWAEM+EXCM-PP%DEXCCOREM
               E%PAWCOREM=E%PAWCOREM+PP%DEXCCOREM
            ENDIF
  !-----------------------------------------------------------------------
  ! calculate the local radial potential
  ! mind in the non-collinear case the potential V(r) = d E(r) / d rho (r)
  ! and the potential vec mu(r) = d E(r) / d vec m (r) are stored in
  ! POT and POTAE (potentials need to be real), whereas 
  ! in the collinear case the spin up and down potentials
  ! are stored in POT and POTAE
  ! probably (1._q,0._q) should rewrite this in the collinear case
  !-----------------------------------------------------------------------
    ! initialise the spin orbit contributions to D_ij to 0
         CSO=0

         IF ( WDES%LNONCOLLINEAR ) THEN
            CALL RAD_MAG_DENSITY( RHO, RHOCOL, LYMAX, PP%R)
            ! calculate aspherical energy first to make sure the potential
            ! is correct at the end
! do LDA+U instead of LSDA+U
            IF (L_NO_LSDA()) RHOCOL(:,:,2:WDES%NCDIJ)=0
! 
            IF (LASPH) THEN
               CALL RAD_GGA_ASPH( PP%R, 2, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
                    RHOCOL, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCGA)
               E%PAWPSAS=E%PAWPSAS-EXCGA
            ENDIF
            CALL RAD_POT( PP%R, 2, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
               RHOCOL, PP%RHOPS, PP%POTPS, POT, DOUBLEPS, EXCG)
            E%PAWPSG=E%PAWPSG-EXCG
            CALL RAD_MAG_DIRECTION( RHO, RHOCOL, POT, LYMAX, PP%R)

            CALL RAD_MAG_DENSITY( RHOAE, RHOCOL, LYMAX, PP%R)
! do LDA+U instead of LSDA+U
            IF (L_NO_LSDA()) RHOCOL(:,:,2:WDES%NCDIJ)=0
!             
            IF (LASPH) THEN 
               CALL RAD_GGA_ASPH( PP%R, 2, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
                    RHOCOL, PP%RHOAE, PP%POTAE,  POTAE, DOUBLEAE,EXCGA)
               E%PAWAEAS=E%PAWAEAS+EXCGA-PP%DEXCCORE
            ENDIF
            CALL RAD_POT( PP%R, 2, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
                 RHOCOL, PP%RHOAE, PP%POTAE,  POTAE, DOUBLEAE,EXCG)
            E%PAWAEG=E%PAWAEG+EXCG-PP%DEXCCORE
            E%PAWCORE=E%PAWCORE+PP%DEXCCORE

            IF (WDES%LSORBIT) &
              CALL SPINORB_STRENGTH(POTAE(:,1,:), PP%RHOAE, PP%POTAE, PP%R, CSO, &
                PP%LMAX, PP%LPS ,PP%WAE, PP%ZCORE+PP%ZVALF_ORIG, THETA=BETA, PHI=ALPHA)
            
            CALL RAD_MAG_DIRECTION( RHOAE, RHOCOL, POTAE, LYMAX, PP%R)
            
         ELSE
    ! collinear case
            DRHOCORE=0
! do LDA+U instead of LSDA+U
            RHOCOL=RHO
            IF (L_NO_LSDA()) RHOCOL(:,:,2)=0
!
            IF (LASPH) THEN
               CALL RAD_GGA_ASPH( PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
                    RHOCOL, PP%RHOPS-DRHOCORE(1:RNMAX), PP%POTPS, POT, DOUBLEPS,EXCGA)
               E%PAWPSAS=E%PAWPSAS-EXCGA
            ENDIF
            ! cl shifts DRHOCORE is (0._q,0._q) here, since for the pseudo terms
            ! we do not include the core electron in the exchange correlation term
            ! but only in the Hartree term
            CALL RAD_POT( PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
               RHOCOL, PP%RHOPS-DRHOCORE(1:RNMAX), PP%POTPS, POT, DOUBLEPS, EXCG)
            CALL SET_CL_DRHOCORE_PS(DRHOCORE, NT, PP%R, PP%AUG)
!RCORE_MAIN
!           this statement is usefull in conjunction with
!           #define RCORE_MAIN XXX
!            CALL SETAUG_CL( DRHOCORE, NT, PP)
!RCORE_MAIN
            CALL ADD_CL_HARTREE_POT(DRHOCORE, NT, ISPIN, POT, PP%R)

            E%PAWPSG=E%PAWPSG-EXCG

            DRHOCORE=0
! do LDA+U instead of LSDA+U
            RHOCOL=RHOAE
            IF (L_NO_LSDA()) RHOCOL(:,:,2)=0
!
            CALL SET_CL_DRHOCORE_AE(DRHOCORE, NT, PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG )
!RCORE_MAIN or soft augmentation charge
!           this statement allows to use a augmentation like charge
!           for the AE part
!            CALL SETAUG_CL( DRHOCORE, NT, PP)
!RCORE_MAIN or soft augmentation charge
            IF (LASPH) THEN
               CALL RAD_GGA_ASPH( PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
                    RHOCOL, PP%RHOAE-DRHOCORE(1:RNMAX), PP%POTAE,  POTAE, DOUBLEAE,EXCGA)
               E%PAWAEAS=E%PAWAEAS+EXCGA-PP%DEXCCORE
            ENDIF
            CALL RAD_POT( PP%R, ISPIN, LEXCH, LEXCHG, LYMAX, PP%LMAX_CALC,  &
               RHOCOL, PP%RHOAE-DRHOCORE(1:RNMAX), PP%POTAE,  POTAE, DOUBLEAE,EXCG)

            CALL ADD_CL_HARTREE_POT(DRHOCORE, NT, ISPIN, POTAE, PP%R)

            E%PAWAEG=E%PAWAEG+EXCG-PP%DEXCCORE
            E%PAWCORE=E%PAWCORE+PP%DEXCCORE
         ENDIF

         ! subtract from double counting corrections the atomic contributions
         ! (remember our reference are the atoms)
         DOUBLEC_PS= DOUBLEC_PS-DOUBLEPS
         DOUBLEC_AE= DOUBLEC_AE+DOUBLEAE-PP%DEXCCORE

  !-----------------------------------------------------------------------
  ! calculate core level shift for up potential (collinear case) or
  ! total potential (non collinear case, 
  ! SPIN IS STILL MISSING !!!!!!!!!!!!!
  !-----------------------------------------------------------------------

         IF (.NOT. LCOREL) THEN
           CALL RAD_CL_SHIFT( POT(:,:,1)/SPI2, POTAE(:,:,1)/SPI2, PP%R, CL_SHIFT(1,NI), PP%AUG(:,0))
           IF (NCDIJ==2) THEN
            ! collinear spinpolarised case
            ! average spin up and down potential
              CALL RAD_CL_SHIFT( POT(:,:,2)/SPI2, POTAE(:,:,2)/SPI2, PP%R, CL_SHIFT(1,NI), PP%AUG(:,0))
              CL_SHIFT(:,NI)=CL_SHIFT(:,NI)/2
           ENDIF
         ELSE
           CALL CALCULATE_MAX_N_L(PP%ZCORE, MAXNL)

           ALLOCATE( W(RNMAX,MAXNL), N(MAXNL), LC(MAXNL), EIG(MAXNL))

           ! first version for cl-shifts
           ! ---------------------------
           ! calculate core wavefunctions in PAW sphere for atomic reference potential

          ! CL_SHIFT(1:MAXNL,NI) =0
          ! CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
          !   W, N, LC, EIG)
          ! IF (MAXNL > SIZE(CL_SHIFT,DIM=1)) THEN
          !    WRITE(0,*) 'internal error: increase CL_MAXNL in main.F',MAXNL,SIZE(CL_SHIFT,DIM=1)
          !    STOP
          ! ENDIF

           ! now calculate the first order change caused by the current potential
           ! and subtract the pseudo contribution
          ! IF (NCDIJ==2) THEN
          !   CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ (POT(1:RNMAX,1,1)+POT(1:RNMAX,1,2))/2/SPI2 , &
          !         PP%R, CL_SHIFT(:,NI), &
          !         PP%AUG(:,0),  MAXNL, W, EIG, (POTAE(:,1,1)+POTAE(:,1,2))/SPI2/2)
          ! ELSE
          !   CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ POT(1:RNMAX,1,1)/SPI2 , &
          !         PP%R, CL_SHIFT(:,NI), &
          !         PP%AUG(:,0),  MAXNL, W, EIG, POTAE(:,1,1)/SPI2)
          ! ENDIF
          !  WRITE(0,*) CL_SHIFT(1:MAXNL,NI)
           CL_SHIFT(1:MAXNL,NI) =0

           ! version for cl-shifts that uses exact potential
           ! ------------------------------------------------
           ! solve radial Schroedinger equation for *current* potential
           IF (NCDIJ==2) THEN

              CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
              W, N, LC, EIG, (POTAE(1:RNMAX,1,1)+POTAE(1:RNMAX,1,2))/2 )

           ! subtract only the pseudo contribution
              CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ (POT(1:RNMAX,1,1)+POT(1:RNMAX,1,2))/2/SPI2 , &
                   PP%R, CL_SHIFT(:,NI), &
                   PP%AUG(:,0),  MAXNL, W, EIG)

           ELSE
              CALL SET_CORE_WF( PP%RHOAE, PP%POTAE , PP%R, PP%ZCORE, PP%ZVALF_ORIG , &
              W, N, LC, EIG, POTAE(1:RNMAX,1,1))


           ! subtract only the pseudo contribution
              CALL RAD_CL_SHIFT_AE( (PP%POTPSC-PP%POTPS)+ POT(1:RNMAX,1,1)/SPI2 , &
                   PP%R, CL_SHIFT(:,NI), &
                   PP%AUG(:,0),  MAXNL, W, EIG)
           ENDIF
           ! WRITE(0,*)  CL_SHIFT(1:MAXNL,NI)
           DEALLOCATE(W, N, LC, EIG)
         ENDIF

  !-----------------------------------------------------------------------
  ! calculate the correction terms to D
  ! I have defined the PAW contribution in a way that in the limit of
  ! atomic occupancies no contributions are added
  !-----------------------------------------------------------------------

         CALL RAD_POT_WEIGHT( PP%R, NCDIJ, LYMAX, POTAE)
         CALL RAD_POT_WEIGHT( PP%R, NCDIJ, LYMAX, POT)
         CTMP=0
         DO ISP=1,NCDIJ
            DDLM=0

            CALL RAD_PROJ(  POTAE(:,:,ISP), PP%R, 1._q, DDLM, PP%LMAX, PP%LPS, PP%WAE )
            CALL RAD_PROJ(  POT(:,:,ISP)  , PP%R,-1._q, DDLM, PP%LMAX, PP%LPS, PP%WPS )
            CALL RAD_AUG_PROJ( POT(:,:,ISP), PP%R, DDLM, PP%LMAX, PP%LPS, &
                  LYMAX, PP%AUG, PP%QPAW )
    ! transform them using Clebsch Gordan coefficients and add to CDIJ
            CALL TRANS_DLM( CTMP(:,:,ISP), DDLM , PP )
         ENDDO
    
         ! non-collinear case: strength parameters need to go to the spinor presentation now
         IF (WDES%LNONCOLLINEAR) CALL DIJ_FLIP(CTMP,LMDIM)

!OBengone modify start
         ! correction terms from LDA+U         
         IF (USELDApU()) THEN            
            IF (WDES%ISPIN==2) CALL DIJ_FLIP2(COCC,LMDIM)  ! calculate up and down occupancies

            IF (WDES%LNONCOLLINEAR) CALL DIJ_FLIP4(COCC,LMDIM) ! take occupancies to spinor representation

!           IF (WDES%LNONCOLLINEAR) THEN
!               COCC=REAL(COCC)
!               COCC=CRHODE(:,:,NIP,:)
!               DO ISP=1,WDES%NCDIJ
!               WRITE(*,*) 'spin',ISP
!               CALL DUMP_DLLMM_IM(COCC(:,:,ISP),PP)
!               CALL DUMP_DLLMM_IM(CRHODE(:,:,NIP,ISP),PP)
!               WRITE(*,*) 'spin',ISP
!               CALL DUMP_DLLMM(COCC(:,:,ISP),PP)
!               CALL DUMP_DLLMM(CRHODE(:,:,NIP,ISP),PP)
!               ENDDO 
!               STOP
!               CALL DIJ_FLIP4(COCC,LMDIM) ! occupancies need to go to spinor representation
!           ENDIF

            IF (WDES%ISPIN==1.AND.(.NOT.WDES%LNONCOLLINEAR)) THEN
               COCC(:,:,1)=COCC(:,:,1)/2
               COCC(:,:,2)=COCC(:,:,1)
            ENDIF
            ! handle the occupancies to the LDA+U code            
            CALL LDAPLUSU(LMDIM,NI,NT,COCC, CTMP, PP, DOUBLEC_LDAU)
            DOUBLEC_AE = DOUBLEC_AE + DOUBLEC_LDAU
         ENDIF
!OBengone modify end           

         IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) THEN
            COCC=CRHODE(:,:,NIP,:)
            CALL DIJ_FLIP4(COCC,LMDIM) ! go to spinor representation
            CALL CALC_ORBITAL_MOMENT(LMDIM,NI,NT,COCC,PP,0._q,0._q)
         ENDIF

         CDIJ(:,:,NIP,:)=CDIJ(:,:,NIP,:)+CTMP+CSO
      ENDDO ion
!=======================================================================
! now distribute the DIJ to all  which hold DIJ (using global sum)
!=======================================================================
      
      

      
      
      
      
      
      IF ( LMETA ) THEN
         
         
         
      ENDIF
      IF ( LASPH ) THEN
         
         
      ENDIF

! deallocate kinedens if metagga
      IF (LMETA) DEALLOCATE(KINDENSAE,KINDENSPS,RHOUPD,RHOLMUPD)

      DEALLOCATE(RHOCOL)
      DEALLOCATE( POTAE, RHOAE, POT, RHO, DRHOCORE )

      E%PAWAE=DOUBLEC_AE
      E%PAWPS=DOUBLEC_PS
    END SUBROUTINE SET_DD_PAW


!*******************************************************************
!
! a few small helper routines
!
!*******************************************************************

    SUBROUTINE DIJ_FLIP(CTMP,LMDIM)
      USE prec
      IMPLICIT NONE
      INTEGER LMDIM
      COMPLEX(q) CTMP(:,:,:)
      COMPLEX(q) C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER L,LP

      FAC=1.0
      DO L=1,LMDIM
         DO LP=1,LMDIM
            C00=CTMP(L,LP,1)
            CX =CTMP(L,LP,2)
            CY =CTMP(L,LP,3)
            CZ =CTMP(L,LP,4)
            
            CTMP(L,LP,1)= (C00+CZ)*FAC           
            CTMP(L,LP,2)= (CX-CY*(0._q,1._q))*FAC
            CTMP(L,LP,3)= (CX+CY*(0._q,1._q))*FAC
            CTMP(L,LP,4)= (C00-CZ)*FAC           
         ENDDO
      ENDDO
    END SUBROUTINE DIJ_FLIP
    
!
! this (1._q,0._q) calculate the up and down occupancies from
! total and magnetisation
!
    SUBROUTINE DIJ_FLIP2(CTMP,LMDIM)
      USE prec
      IMPLICIT NONE
      INTEGER LMDIM
      COMPLEX(q) CTMP(:,:,:)
      COMPLEX(q) CQU,CQD
      REAL(q) :: FAC
      INTEGER L,LP

      FAC=0.5
      DO L=1,LMDIM
         DO LP=1,LMDIM
            CQU=CTMP(L,LP,1)
            CQD=CTMP(L,LP,2)
            CTMP(L,LP,1)=FAC*(CQU+CQD)
            CTMP(L,LP,2)=FAC*(CQU-CQD)
         ENDDO
      ENDDO
    END SUBROUTINE DIJ_FLIP2

    
    SUBROUTINE DIJ_FLIP4(CTMP,LMDIM)
      USE prec
      IMPLICIT NONE
      INTEGER LMDIM
      COMPLEX(q) CTMP(:,:,:)
      COMPLEX(q) C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER L,LP

      FAC=0.5
      DO L=1,LMDIM
         DO LP=1,LMDIM
            C00=CTMP(L,LP,1)
            CX =CTMP(L,LP,2)
            CY =CTMP(L,LP,3)
            CZ =CTMP(L,LP,4)
            
            CTMP(L,LP,1)= (C00+CZ)*FAC           
            CTMP(L,LP,2)= (CX-CY*(0._q,1._q))*FAC
            CTMP(L,LP,3)= (CX+CY*(0._q,1._q))*FAC
            CTMP(L,LP,4)= (C00-CZ)*FAC           
         ENDDO
      ENDDO
    END SUBROUTINE DIJ_FLIP4
    
   
!*******************************************************************
!
!  calculate net moment of the soft
!  augmentation charges (compensation charges in PAW language)
!  RHO(LM) at (1._q,0._q) site (this routine is called from us.F)
!
!  RHO(lm,l'm') are the occupancies of each channel
!  RHO(LM) is given by
!  RHO(LM) =  sum C(LM,ll',mm')  RHO(lm,l'm') QPAW(llp)
!  also see TRANS_RHOLM
!*******************************************************************


    SUBROUTINE CALC_RHOLM( LYMAX, RHOLLMM, RHOLM, P)
      USE pseudo
      USE constant
      USE asa
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) RHOLLMM(:,:)   ! net augmentation charge
      REAL(q) RHOLM(:)       ! local charge for each L,M
      INTEGER LYMAX          ! maximum L
    ! local varible
      INTEGER LMYMAX,CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      REAL(q) FAKT

    ! initialize everything to 0

      LMYMAX=(LYMAX+1)**2
      RHOLM=0

   ! loop over all channels (l,epsilon)
      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)

   ! transform coefficients and multiply with QPAW (which gives the moment
   ! of the corresponding charge)
      FAKT=1
      IF (CH1 /= CH2) THEN
         FAKT=2   ! CH2 is >= CH1
      ENDIF

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)

         DO  IC=ISTART,IEND-1
            RHOLM(JS(IC))= RHOLM(JS(IC))+ YLM3(IC)*FAKT*P%QPAW(CH1,CH2,JL(IC))* &
     &         RHOLLMM(LM+M-1,LMP+MP-1)
         ENDDO
      ENDDO
      ENDDO

      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO


      END SUBROUTINE CALC_RHOLM



!*******************************************************************
!
!  transform the real part of the occupancies RHO(lm,l'm') 
!  to RHO(llp,L,M) using Clebsch Gordan coefficients
!
!' RHO(ll',LM) =  sum C(LM,ll',mm')  RHO(lm,l'm')
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of RHO(llp,LM) is somewhat akward
!  for each l lp pair, (2l+1) (2lp+1) elements must be stored
!  they are stored in the order
!    Lmin,M=0 ... Lmin,M=2*Lmin+1
!     ...
!    Lmax,M=0 ... Lmax,M=2*Lmax+1
!  where Lmin and Lmax are given by the triangular rule
!  Lmin = | l-l' |  and Lmax = | l+l'|
!  certain elements in this array will be always (0._q,0._q), because
!  the sum rule allows only L=Lmin,Lmin+2,...,Lmax
!*******************************************************************


    SUBROUTINE TRANS_RHOLM( RHOLLMM, RHOLM, P)
      USE pseudo
      USE asa
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) RHOLLMM(:,:)   ! net augmentation charge
      REAL(q) RHOLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX
      REAL(q) FAKT
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! transform coefficients
      FAKT=1
      IF (CH1 /= CH2) THEN
         FAKT=FAKT*2
      ENDIF
   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN

      RHOLM(IBASE+1:IBASE+(2*LL+1)*(2*LLP+1))=0

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            RHOLM(JS(IC)+JBASE)= RHOLM(JS(IC)+JBASE)+ YLM3(IC)*FAKT* &
     &         RHOLLMM(LM+M-1,LMP+MP-1)
         ENDDO
      ENDDO
      ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
    END SUBROUTINE TRANS_RHOLM


!*******************************************************************
!
!  transform the imaginary part of the occupancies RHO(lm,l'm') 
!  to RHO(llp,L,M) using Clebsch Gordan coefficients
!
!' RHO(ll',LM) =  sum C(LM,ll',mm')  RHO(lm,l'm')
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of RHO(llp,LM) is somewhat akward
!  for each l lp pair, (2l+1) (2lp+1) elements must be stored
!  they are stored in the order
!    Lmin,M=0 ... Lmin,M=2*Lmin+1
!     ...
!    Lmax,M=0 ... Lmax,M=2*Lmax+1
!  where Lmin and Lmax are given by the triangular rule
!  Lmin = | l-l' |  and Lmax = | l+l'|
!  certain elements in this array will be always (0._q,0._q), because
!  the sum rule allows only L=Lmin,Lmin+2,...,Lmax
!*******************************************************************

    SUBROUTINE TRANS_RHOLM_IM( RHOLLMM, RHOLM, P)
      USE pseudo
      USE asa
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) RHOLLMM(:,:)   ! net augmentation charge
      REAL(q) RHOLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX
      REAL(q) FAKT,FAKT2
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)
      
      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! transform coefficients
      FAKT=1
      IF (CH1 /= CH2) THEN
         FAKT=FAKT*2
      ENDIF
   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN

      RHOLM(IBASE+1:IBASE+(2*LL+1)*(2*LLP+1))=0

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1
         IF (M<MP) THEN
           FAKT2=1
         ELSE
           FAKT2=-1
         ENDIF
         
         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            RHOLM(JS(IC)+JBASE)= RHOLM(JS(IC)+JBASE)+ YLM3(IC)*FAKT*FAKT2* &
     &         AIMAG(RHOLLMM(LM+M-1,LMP+MP-1))
         ENDDO
      ENDDO
      ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
    END SUBROUTINE TRANS_RHOLM_IM
    
    SUBROUTINE TRANS_RHOLMI( DLLMM, DLM, P)
      USE pseudo
      USE asa
      USE constant
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)   ! net augmentation charge
      REAL(q) DLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX,INMIN,INMAX

   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)+ &
                            DLM(JS(IC)+JBASE)*YLM3I(IC)
         ENDDO
      ENDDO
      ENDDO
  ! fill symmetric components (CH2 is >= CH1)
      IF  (CH1 /= CH2) THEN
         DO M =1,2*LL+1
         DO MP=1,2*LLP+1
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)/2
            DLLMM(LMP+MP-1,LM+M-1)=DLLMM(LM+M-1,LMP+MP-1)
         ENDDO
         ENDDO
      ENDIF

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE TRANS_RHOLMI
    
    SUBROUTINE TRANS_RHOLMI_IM( DLLMM, DLM, P)
      USE pseudo
      USE asa
      USE constant
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)   ! net augmentation charge
      REAL(q) DLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX,INMIN,INMAX
      REAL(q) :: FAKT2
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1
         IF (M<MP) THEN
           FAKT2= 1
         ELSE
           FAKT2=-1
         ENDIF

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)+ &
                            DLM(JS(IC)+JBASE)*YLM3I(IC)*FAKT2
         ENDDO
      ENDDO
      ENDDO
  ! fill symmetric components (CH2 is >= CH1)
      IF  (CH1 /= CH2) THEN
         DO M =1,2*LL+1
         DO MP=1,2*LLP+1
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)/2
            DLLMM(LMP+MP-1,LM+M-1)=-DLLMM(LM+M-1,LMP+MP-1)
         ENDDO
         ENDDO
      ENDIF

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE TRANS_RHOLMI_IM
    
!*******************************************************************
!
!  transform D(llp,L,M) to the representation D(lm,l'm')
!  using Clebsch Gordan coefficients and add to another array
!
!  D(lm,l'm') =  sum C(LM,ll',mm') D(llp,LM)
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of D(llp,LM) is somewhat akward see above
!
!*******************************************************************


    SUBROUTINE TRANS_DLM( DLLMM, DLM, P)
      USE pseudo
      USE asa
      USE constant
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)   ! net augmentation charge
      REAL(q) DLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX,INMIN,INMAX

   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN

      INMIN=1000
      INMAX =-1000

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            INMIN=MIN(JS(IC)+JBASE,INMIN)
            INMAX=MAX(JS(IC)+JBASE,INMAX)
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)+ &
                            DLM(JS(IC)+JBASE)*YLM3(IC)
         ENDDO
      ENDDO
      ENDDO
!      WRITE(0,*) IBASE,IBASE+(2*LL+1)*(2*LLP+1),INMIN,INMAX
  ! fill symmetric components (CH2 is >= CH1)
      IF  (CH1 /= CH2) THEN
         DO M =1,2*LL+1
         DO MP=1,2*LLP+1
            DLLMM(LMP+MP-1,LM+M-1)=DLLMM(LM+M-1,LMP+MP-1)
         ENDDO
         ENDDO
      ENDIF


      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE TRANS_DLM

!*******************************************************************
!
!  transform D(LM) to D(lm,l'm')
!  where D(LM) is defined as
!  D(LM) = \int V Y(L,M) Q(L)(r) d^3 r
!  and \int Q(L)(r) r^2 dr is integrating to 1
!  (this routine is called from us.F)
!
!  the "pseudopotential strenght" is then given by
!  D(lm,l'm') =  sum C(LM,ll',mm')  D(LM) QPAW(llp)
!  also see TRANS_DLM
!*******************************************************************


    SUBROUTINE CALC_DLLMM( DLLMM, DLM, P)
      USE pseudo
      USE asa
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)  ! net augmentation charge
      REAL(q) DLM(:)      ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP

    ! initialize everything to 0
      DLLMM=0

   ! loop over all channels (l,epsilon)
      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

   ! quantum numbers l and lp of these two channels
      LL =P%LPS(CH1)
      LLP=P%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)

   ! transform coefficients and multiply with QPAW (which gives the moment
   ! of the corresponding charge)

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)

         DO  IC=ISTART,IEND-1
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)+ &
                      DLM(JS(IC))*YLM3(IC)*P%QPAW(CH1,CH2,JL(IC))
         ENDDO
      ENDDO
      ENDDO
  ! fill symmetric components (CH2 is >= CH1)
      IF  (CH1 /= CH2) THEN
         DO M =1,2*LL+1
         DO MP=1,2*LLP+1
            DLLMM(LMP+MP-1,LM+M-1)=DLLMM(LM+M-1,LMP+MP-1)
         ENDDO
         ENDDO
      ENDIF

      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
    END SUBROUTINE CALC_DLLMM

!*******************************************************************
!
! write out DIJ
!
!*******************************************************************


    SUBROUTINE DUMP_DLLMM( DLLMM, P)
      USE pseudo
      USE constant
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)   ! net augmentation charge
    ! local varible
      INTEGER LM,LMP

      WRITE(*,*) 'pseudopotential strength DIJ is'
      DO LM=1,P%LMMAX
         WRITE(*,'(18(F9.5,1X))') (REAL(DLLMM(LMP,LM),q),LMP=1,MIN(18,P%LMMAX))
      ENDDO

    END SUBROUTINE DUMP_DLLMM


    SUBROUTINE DUMP_DLLMM_IM( DLLMM, P)
      USE pseudo
      USE constant
      IMPLICIT NONE
      TYPE (potcar) P
      COMPLEX(q) DLLMM(:,:)   ! net augmentation charge
    ! local varible
      INTEGER LM,LMP
      WRITE(*,*) 'pseudopotential strength DIJ is'
      DO LM=1,P%LMMAX
         WRITE(*,'(18(F9.5,1X))') (AIMAG(DLLMM(LMP,LM)),LMP=1,MIN(18,P%LMMAX))
      ENDDO

    END SUBROUTINE DUMP_DLLMM_IM

!*******************************************************************
!
! retrieve the elements of the occupancy channels 
! that are mixed from RHOLM and store
! them in continous order in an array RHOLM_STORE
! the number of transfered elements in returned in IADD
! LCOMPACT = T return in a compact storage (i.e. no (0._q,0._q) elements
!              inbetween)
! LCOMPACT = F return in conventional mode (storage layout used in
!              radial.F and paw.F)
!
!*******************************************************************

    SUBROUTINE STORE_RHOLM( RHOLM, RHOLM_STORE, METRIC, IADD, P, &
              LCOMPACT, IBASE)
      USE pseudo
      IMPLICIT NONE
      TYPE (potcar) P
      REAL(q) RHOLM(:)       ! full RHO(llp,L,M) matrix
      REAL(q) RHOLM_STORE(:) ! storage point for elements with L=0
      REAL(q) METRIC(:)
      INTEGER IADD           ! number of elements stored in RHOLM_STORE
      INTEGER IBASE          ! maximum number of elements in RHOLM
      LOGICAL LCOMPACT       ! compact mode or not
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LMIN,LMAX,JBASE,LMAIN,MMAIN,LMMAIN

      IADD=0
      IBASE=0

      DO CH1=1,P%LMAX
      DO CH2=CH1,P%LMAX

       ! quantum numbers l and lp of these two channels
         LL =P%LPS(CH1)
         LLP=P%LPS(CH2)
       ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=MIN(LMAX_MIX,P%LMAX_CALC,ABS(LL+LLP))
         JBASE=IBASE-LMIN*LMIN

         LMMAIN=LMIN*LMIN
         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            IF (LCOMPACT) THEN
               LMMAIN=LMMAIN+1
            ELSE
               LMMAIN=LMAIN*LMAIN+MMAIN
            ENDIF

            IADD=IADD+1
            RHOLM_STORE(IADD)=RHOLM(JBASE+LMMAIN)*METRIC(IADD)
         ENDDO
         ENDDO

         IF (LCOMPACT) THEN
            DO LMAIN=LMIN,ABS(LL+LLP),2
               IBASE=IBASE+LMAIN*2+1
            ENDDO
         ELSE
            IBASE=IBASE+(2*LL+1)*(2*LLP+1)
         ENDIF
      ENDDO
      ENDDO

    END SUBROUTINE STORE_RHOLM

!*******************************************************************
!
! retrieve the elements of the occupancy channels that are mixed
! from RHOLM_STORE and store them in conventional order in RHOLM
! the number of transfered elements is returned in IADD
! the routine supports two mode
! LCOMPACT = T return in a compact storage (i.e. no (0._q,0._q) elements
!              inbetween)
! LCOMPACT = F return in conventional mode (storage layout used in
!              radial.F and paw.F)
!
!*******************************************************************


    SUBROUTINE RETRIEVE_RHOLM( RHOLM, RHOLM_STORE, METRIC, IADD, P, &
              LCOMPACT, IBASE)
      USE pseudo
      IMPLICIT NONE
      TYPE (potcar) P
      REAL(q) RHOLM(:)       ! full RHO(llp,L,M) matrix
      REAL(q) RHOLM_STORE(:) ! storage point for elements with L=0
      REAL(q) METRIC(:)      ! metric
      INTEGER IADD           ! number of elements transfered from RHOLM_STORE
      INTEGER IBASE          ! maximum number of elements in RHOLM
      LOGICAL LCOMPACT       ! compact mode or not
    ! local variable
      INTEGER CH1,CH2,LL,LLP,LMIN,LMAX,JBASE,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IADD=0
      IBASE=0

      DO CH1=1,P%LMAX
      DO CH2=CH1,P%LMAX

       ! quantum numbers l and lp of these two channels
         LL =P%LPS(CH1)
         LLP=P%LPS(CH2)
       ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=MIN(LMAX_MIX,P%LMAX_CALC,ABS(LL+LLP))
         JBASE=IBASE-LMIN*LMIN

         LMMAIN=LMIN*LMIN

         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            IF (LCOMPACT) THEN
               LMMAIN=LMMAIN+1
            ELSE
               LMMAIN=LMAIN*LMAIN+MMAIN
            ENDIF
               
            IADD=IADD+1
            RHOLM(JBASE+LMMAIN)=RHOLM_STORE(IADD)/METRIC(IADD)
         ENDDO
         ENDDO

         IF (LCOMPACT) THEN
            DO LMAIN=LMIN,ABS(LL+LLP),2
               IBASE=IBASE+2*LMAIN+1
            ENDDO
         ELSE
            IBASE=IBASE+(2*LL+1)*(2*LLP+1)
         ENDIF
      ENDDO
      ENDDO

    END SUBROUTINE RETRIEVE_RHOLM
  END MODULE paw
