#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: electron.F,v 1.12 2003/06/27 13:22:15 kresse Exp kresse $
!
! subroutine for performing electronic minimization in VASP
! I had to remove this subroutine from the main.F in order
! to get things compiled in the T3D in 16 Mwords
!
!**********************************************************************

      SUBROUTINE ELMIN( &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,WUP,WDW,LATT_CUR,LATT_INI,EXCTAB, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,DIP,E,E2, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,DESUM,XCSIF, &
          NSTEP,NELMLS,LMDIM,NIOND,IRDMAX,NBLK,NEDOS, &
          TOTEN,TOTENL,EFERMI,LDIMP,LMDIMP,LTRUNC)

      USE prec
      USE charge
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE pot
      USE force
      USE fileio
      USE nonl
      USE nonlr
      USE rmm_diis
      USE david
      USE ini
      USE ebs
!      USE rot
      USE dfast
      USE choleski
      USE mwavpre
      USE mwavpre_noio
      USE msphpro
      USE broyden
      USE msymmetry
      USE subrot
      USE melf
      USE base
      USE mpimy
      USE mgrid
      USE mkpoints
      USE constant
      USE setexm
      USE poscar
      USE wave
      USE hamil
      USE paw
      USE cl
      USE vaspxml
!      USE pawfock
!-MM- Added to accomodate constrained moments
      USE Constrained_M_modular
!-MM- end of additions
!-MM- Added for writing LDA+U occupancies
      USE LDAPLUSU_MODULE
!-MM- end of additions

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
!      IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavefun)     W_F        ! wavefunction for all bands simultaneous
      TYPE (wavefun)     W_G        ! same as above
      TYPE (wavefun)     WUP
      TYPE (wavefun)     WDW
      TYPE (latt)        LATT_CUR
      TYPE (exctable)    EXCTAB
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID       ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
      TYPE (grid_3d)     GRIDB      ! Broyden grid
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (dipol)       DIP
      TYPE (energy)      E,E2
      TYPE (latt)        LATT_INI
     
      INTEGER NSTEP,NELMLS,LMDIM,NIOND,IRDMAX,NBLK,NEDOS
      REAL(q) TOTEN,TOTENL,EFERMI

      COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
      COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
      RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
      COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!-----non-local pseudopotential parameters
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!---- needed temporary for aspherical GGA calculation 
      OVERLAP,ALLOCATABLE ::  CDIJ_TMP(:,:,:,:)
!-----paw sphere charge density
      INTEGER N_MIX_PAW
      REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!-----charge-density and potential on small grid
      COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      RGRID       SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
!-----description how to go from one grid to the second grid
!-----density of states
      REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
!-----local l-projected wavefunction characters (not really used here)
      REAL(q)    PAR(1,1,1,1,1),DOSPAR(1,1,1,1)
!  all-band-simultaneous-update arrays
      GDEF       CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS), &
                 CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS)
!----- energy at each step
      REAL(q)   DESUM(500)
      REAL(q)   XCSIF(3,3),DESUM1
      INTEGER  :: IONODE, NODE_ME

      REAL(q), EXTERNAL :: RHO0
      INTEGER N,ISP,ICONJU,IROT,ICEL,I,II,IRDMAA, &
              IERR,IDUM,IFLAG,ICOUEV,ICOUEV2,NN,NORDER,IERRBR,L,LP, &
              NCOUNT
      REAL(q) TV,TV0,TC,TC0,TVPUL,TVPUL0,TCPUL,TCPUL0, &
              BTRIAL,RDUM,RMS,ORT,TOTEN2,RMS2,RMST, &
              WEIGHT,BETATO,DESUM2,RMSC,RMSP
      REAL(q) RHOAUG(WDES%NCDIJ),RHOTOT(WDES%NCDIJ)
      COMPLEX(q) CDUM
      CHARACTER (LEN=1) CHARAC
      LOGICAL LDELAY
! METAGGA (Robin Hirschl)
      COMPLEX(q),ALLOCATABLE:: TAU(:,:)
      COMPLEX(q),ALLOCATABLE:: TAUW(:,:)      
!-----parameters for sphpro.f
      INTEGER :: LDIMP,LMDIMP,LTRUNC
      REAL(q),ALLOCATABLE:: PAR_DUMMY(:,:,:,:,:)

!R.S
      integer tiu6, tiu0, tiuvtot
        tiu6 = IO%IU6
        tiu0 = IO%IU0
      IONODE=0
      NODE_ME=0
#ifdef MPI
      IONODE  = WDES%COMM%IONODE
      NODE_ME = WDES%COMM%NODE_ME
#endif
      NELM=INFO%NELM
      ! to make timing more sensefull syncronize now
      CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror ))

      CALL VTIME(TVPUL0,TCPUL0)

     io_begin
      IF (INFO%LONESW) THEN
      IF (IO%IU0>=0) &
      WRITE(TIU0,141)
      WRITE(17,141)
  141 FORMAT('       N       E                     dE             ' &
            ,'d eps       ncg     rms          ort')

      ELSE
      IF (IO%IU0>=0) &
      WRITE(TIU0,142)
      WRITE(17,142)
  142 FORMAT('       N       E                     dE             ' &
            ,'d eps       ncg     rms          rms(c)')
      ENDIF
     io_end

      DESUM1=0
      INFO%LMIX=.FALSE.

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'----------------------------------------- Iteration ', &
     &I4,'(',I4,')  ---------------------------------------'//)
     DWRITE0 'electron entered'

      E%EENTROPY=0
      DOS=0
      DOSI=0
      CALL DIPOL_RESET()

!=======================================================================
      electron: DO N=1,NELM

      CALL XML_TAG("scstep")

!======================================================================
      io_begin
      WRITE(TIU6,140) NSTEP,N
      io_end
!=======================================================================
! if recalculation of total lokal potential is necessary (INFO%LPOTOK=.F.)
! call POTLOK: the subroutine calculates
! ) the hartree potential from the electronic  charge density
! ) the exchange correlation potential
! ) and the total lokal potential
!  in addition all double counting correction and forces are calculated
! &
! call SETDIJ
! calculates the Integral of the depletion charges * local potential
! and sets CDIJ
!=======================================================================

      CALL WVREAL(WDES,GRID,W) ! only for gamma some action
      CALL VTIME(TV0,TC0)
      IF (.NOT. INFO%LPOTOK) THEN
      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  EXCTAB,INFO,P,T_INFO,E,LATT_CUR,DIP, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
      CALL VTIME(TV,TC)
      IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'POTLOK',TV-TV0,TC-TC0
      DWRITE0 'potlok is ok'

      CALL VTIME(TV0,TC0)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX,.0_q,.0_q,.0_q)
#ifndef noPAW
      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), INFO%LEXCH, INFO%LEXCHG, &
         E,  LMETA=.FALSE., LASPH=.FALSE., LCOREL= .FALSE.  )

!      CALL TEST_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
!         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), INFO%LEXCH, INFO%LEXCHG, &
!         E,  .FALSE., .FALSE. , .FALSE. )

#endif
!-MM- write LDA+U occupancy matrices
      IF (USELDApU()) CALL LDAPLUSU_PRINTOCC(WDES,T_INFO%NIONS,T_INFO%ITYP,IO%IU6)
!-MM- end of addition

      CALL VTIME(TV,TC)
      IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'SETDIJ',TV-TV0,TC-TC0
      DWRITE0 'setdij is ok'

      INFO%LPOTOK=.TRUE.
      ENDIF
!======================== SUBROUTINE EDDSPX ============================
!
! these subroutines improve the electronic degrees of freedom
! using band by band schemes
! the Harris functional is used for the calculation
! of the total (free) energy so
! E  =  Tr[ H rho ] - d.c. (from input potential)
!
!=======================================================================

      DESUM1=0
      RMS   =0
      ICOUEV=0

      LDELAY=.FALSE.
      ! if Davidson and RMM are selected, use Davidsons algorithm during
      ! delay phase
      IF (INFO%LRMM  .AND. INFO%LDAVID .AND. (N <= ABS(INFO%NELMDL) .OR. N==1)) LDELAY=.TRUE.
      ! if LDELAY is set, subspace rotation and orthogonalisations can be bypassed
      ! since they are done by the Davidson algorithm

!
! sub space rotation before eigenvalue optimization
!
      IF (INFO%LPDIAG .AND. .NOT. LDELAY ) THEN

        IF (INFO%LDIAG) THEN
           IFLAG=3    ! exact diagonalization
        ELSE
           IFLAG=4    ! using Loewdin perturbation theory
        ENDIF
        IF (INFO%IALGO==4) THEN
           IFLAG=0
        ENDIF
        IF (N < ABS(INFO%NELMDL)) IFLAG=13
        CALL VTIME(TV0,TC0)
        CALL EDDIAG(GRID,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, IFLAG,INFO%LOVERL,INFO%LREAL,NBLK,SV, &
             IO%IU0,E%EXHF,.FALSE.)

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0)WRITE(TIU6,2300)'EDDIAG',TV-TV0,TC-TC0

        CALL XML_TIMING(TV-TV0, TC-TC0, "diag")

        DWRITE0 "eddiag is ok"

      ENDIF

      CALL VTIME(TV0,TC0)

      select_algo: IF (INFO%LRMM .AND. .NOT. LDELAY) THEN
!
! RMM-DIIS alogrithm
!

        ! CALL EDEXP(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,WUP,WDES, &
        !  LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV)
        ! GOTO 123

        CALL EDDRMM(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, IO%IU6,IO%IU0, &
             N < ABS(INFO%NELMDL)-ABS(INFO%NELMDL)/4)
        ! previous line selects  special algorithm during delay

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'RMM-DIIS',TV-TV0,TC-TC0
        CALL XML_TIMING(TV-TV0, TC-TC0, "diis")

        DWRITE0 "eddrmm is ok"

      ELSE IF (INFO%LDAVID) THEN
!
! blocked Davidson alogrithm,
!
        NSIM=WDES%NSIM*2
#ifdef MPI
        NSIM=((WDES%NSIM*2+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
        CALL EDDAV(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
              LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, NBLK, IO%IU6,IO%IU0, .FALSE., INFO%LDIAG)

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'EDDAV ',TV-TV0,TC-TC0
        CALL XML_TIMING(TV-TV0, TC-TC0, "dav")

        DWRITE0 "edddav is ok"

      ELSE IF (INFO%IALGO==5 .OR.INFO%IALGO==6 .OR. &
     &         INFO%IALGO==7 .OR.INFO%IALGO==8 .OR. INFO%IALGO==0) THEN select_algo

!
! CG (Teter, Alan, Payne) potential is fixed !!
!
        CALL EDSTEP(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV,IO%IU6, IO%IU0)

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'EDSTEP',TV-TV0,TC-TC0
        CALL XML_TIMING(TV-TV0, TC-TC0, "cg")

        DWRITE0 "edstep is ok"

      ENDIF select_algo
!
! orthogonalise all bands (necessary only for residuum-minimizer)
!
      IF (.NOT.INFO%LORTHO .AND. .NOT. LDELAY) THEN
        CALL VTIME(TV0,TC0)
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ,NBLK)

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'ORTHCH',TV-TV0,TC-TC0
        CALL XML_TIMING(TV-TV0, TC-TC0, "orth")

        DWRITE0 "ortch is ok"
      ENDIF
!
! sub space rotation after eigen value optimization
!
      IF (INFO%LCDIAG .AND. .NOT. LDELAY) THEN

        IF (INFO%LDIAG) THEN
           IFLAG=3
        ELSE
           IFLAG=4
        ENDIF

        CALL VTIME(TV0,TC0)
        CALL REDIS_PW_ALL(WDES, W)
        CALL EDDIAG(GRID,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, IFLAG,INFO%LOVERL,INFO%LREAL,NBLK,SV,IO%IU0, &
             E%EXHF,.FALSE.)

        CALL VTIME(TV,TC)
        IF (IO%NWRITE>=2.AND.IO%IU6>=0)WRITE(TIU6,2300)'EDDIAG',TV-TV0,TC-TC0
        CALL XML_TIMING(TV-TV0, TC-TC0, "diag")
      ENDIF

!=======================================================================
! recalculate the broadened density of states and fermi-weights
! recalculate depletion charge size
!=======================================================================
      CALL VTIME(TV0,TC0)
      CALL MRG_CEL(WDES,W)
      CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
               INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
!=======================================================================
! calculate free-energy and bandstructur-energy
! EBANDSTR = sum of the energy eigenvalues of the electronic states
!         weighted by the relative weight of the special k point
! TOTEN = total free energy of the system
!=======================================================================
      E%EBANDSTR=0.0_q

      DO ISP=1,WDES%ISPIN
      DO I=1,KPOINTS%NKPTS
      DO II=1,WDES%NB_TOT
        E%EBANDSTR=E%EBANDSTR+WDES%RSPIN* REAL( W%CELTOT(II,I,ISP) ,KIND=q) *KPOINTS%WTKPT(I)*W%FERTOT(II,I,ISP)
      ENDDO; ENDDO; ENDDO

      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT
!-MM- Added to accomodate constrained moment calculations
      IF (M_CONSTRAINED()) TOTEN=TOTEN+E_CONSTRAINT()
      io_begin
      CALL WRITE_CONSTRAINED_M(17,.FALSE.)
      io_end
!-MM- end of additions
!---- write total energy to OSZICAR file and stdout
      DESUM(N)=TOTEN-TOTENL
      NELMLS=N
      io_begin
  303 FORMAT('CG : ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
 1303 FORMAT('RMM: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
10303 FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)

      IF (INFO%LRMM .AND. .NOT. LDELAY) THEN
        WRITE(17,  1303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(TIU0, 1303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%LDAVID) THEN
        WRITE(17, 10303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(TIU0,10303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE
        WRITE(17,   303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(TIU0,  303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ENDIF
      IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'DOS   ',TV-TV0,TC-TC0
      io_end
!=======================================================================
!  Test for Break condition
!=======================================================================

      INFO%LABORT=.FALSE.
!-----conjugated gradient eigenvalue and energy must be converged
      IF(ABS(DESUM(N))<INFO%EDIFF.AND.ABS(DESUM1)<INFO%EDIFF) INFO%LABORT=.TRUE.
!-----charge-density not constant and in last cycle no change of charge
      IF (.NOT. INFO%LMIX .AND. .NOT. INFO%LCHCON .AND. MIX%IMIX/=0) INFO%LABORT=.FALSE.
!-----do not stop during the non-selfconsistent startup phase
      IF (N <= ABS(INFO%NELMDL)) INFO%LABORT=.FALSE.
!-----do not stop before minimum number of iterations is reached
      IF (N < ABS(INFO%NELMIN)) INFO%LABORT=.FALSE.
!-----but stop after INFO%NELM steps no matter where we are now
      IF (N>=INFO%NELM) INFO%LABORT=.TRUE.

      IF ((IO%LORBIT>=10).AND.(MOD(N,5)==0).AND.WDES%LNONCOLLINEAR) THEN
         ALLOCATE(PAR_DUMMY(WDES%NB_TOT,WDES%NKDIM,LDIMP,T_INFO%NIONP,WDES%NCDIJ))
         CALL SPHPRO_FAST( &
          GRID,LATT_CUR,LATT_INI, P,T_INFO,W, WDES, 71,IO%IU6,&
          INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,.FALSE., IO%LORBIT,PAR_DUMMY)
         DEALLOCATE(PAR_DUMMY)
      ENDIF
! ======================================================================
! If the end of the electronic loop is reached 
! calculate aspherical contributions to Exc and 
! METAGGA if applicable
!
! Robin Hirschl, 09.02.2001
! ======================================================================

      IF (INFO%LABORT .AND. (INFO%LASPH .OR. INFO%LMETAGGA)) THEN
! calling SET_DD_PAW again causes the array CDIJ to have wrong values
! (the Forces will be wrong in the end)
! This is avaoided by using a temporary array
! Robin Hirschl 11.04.2001 
         ALLOCATE(CDIJ_TMP(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ))
         CDIJ_TMP=CDIJ
           
         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
              WDES%NCDIJ, LMDIM, CDIJ_TMP(1,1,1,1), RHOLM,  &
              CRHODE, INFO%LEXCH, INFO%LEXCHG, &
              E, INFO%LMETAGGA, INFO%LASPH, LCOREL=.FALSE. )
         DEALLOCATE(CDIJ_TMP)
      ENDIF
! ======================================================================
! If the end of the electronic loop is reached
! calculate accurate initial state core level shifts
! if required
! ======================================================================
      IF (INFO%LABORT .AND. ACCURATE_CORE_LEVEL_SHIFTS()) THEN

         ALLOCATE(CDIJ_TMP(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ))
         CDIJ_TMP=CDIJ

         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
              WDES%NCDIJ, LMDIM, CDIJ_TMP(1,1,1,1), RHOLM,  &
              CRHODE, INFO%LEXCH, INFO%LEXCHG, &
              E, LMETA= .FALSE. , LASPH=.FALSE. , LCOREL=.TRUE. )
         DEALLOCATE(CDIJ_TMP)
      ENDIF
!-----------------------------------------------------------------------
! MetaGGA Exc energy on PW grid 
! Robin Hirschl 20001216 
!-----------------------------------------------------------------------
      IF (INFO%LABORT .AND. INFO%LMETAGGA) THEN

         ALLOCATE(TAU(GRIDC%MPLWV,WDES%NCDIJ),TAUW(GRIDC%MPLWV,WDES%NCDIJ))
         
         CALL VTIME(TV0,TC0)
         ! calculate kinetic energy density
         CALL TAU_PW_DIRECT(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
              NIOND,W,WDES,TAU,TAUW)      
         CALL METAEXC_PW(GRIDC,WDES,INFO,E,LATT_CUR,CHTOT,TAU,TAUW,DENCOR)
         CALL VTIME(TV,TC)

         IF (IO%NWRITE>=2.AND.IO%IU6>=0) WRITE(TIU6,2300)'MEXCPW',TV-TV0,TC-TC0
         DEALLOCATE(TAU,TAUW)
      ENDIF

! update total energies if aspherical calculations or METAGGA 
      IF (INFO%LABORT .AND. INFO%LASPH) THEN
         E%TOTENASPH=TOTEN-E%PAWPSG-E%PAWAEG+E%PAWPSAS+E%PAWAEAS
      ENDIF
      IF (INFO%LABORT .AND. INFO%LMETAGGA) &
              E%TOTENMGGA=TOTEN-E%EXCG-E%PAWPSG-E%PAWAEG+E%EXCM+E%PAWPSM+E%PAWAEM
     
!========================= subroutine CHSP  ============================
! if charge density is updated
!  ) first copy current charge to CHTOTL
!  ) set  INFO%LPOTOK to .F. this requires a recalculation of the local pot.
!  ) set INFO%LMIX to .T.
!  ) call subroutine CHSP+ DEPLE to generate the new charge density
!  ) then performe mixing
! MIND:
! ) if delay is selected  do not update
! ) if convergence corrections to forces are calculated do not update charge
!   in last iteration
!=======================================================================
      CALL DEPSUM(W,WDES, LMDIM,CRHODE, INFO%LOVERL)
      CALL US_FLIP(WDES, LMDIM, CRHODE, INFO%LOVERL, .FALSE.)

      INFO%LMIX=.FALSE.
      MIX%NEIG=0

      IF (.NOT. INFO%LCHCON .AND. .NOT. (INFO%LABORT .AND. INFO%LCORR ) &
     &    .AND. N >= ABS(INFO%NELMDL)  ) THEN
      CALL VTIME(TV0,TC0)

      DO ISP=1,WDES%NCDIJ
      CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
      ENDDO
      INFO%LPOTOK=.FALSE.

      RHOLM_LAST=RHOLM

      CALL SET_CHARGE(W, WUP, WDW, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      CALL VTIME(TV,TC)

      IF (IO%NWRITE>=2.AND.IO%IU6>=0)WRITE(TIU6,2300)'CHARGE',TV-TV0,TC-TC0

!-----------------------------------------------------------------------

      IF (MIX%IMIX/=0) THEN
      CALL VTIME(TV0,TC0)
      INFO%LMIX=.TRUE.

      IF (MIX%IMIX==4) THEN
!  broyden mixing ... :
        CALL BRMIX(GRIDB,GRIDC,IO,MIX,B_TO_C, &
           (2*GRIDC%MPLWV),CHTOT,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
           LATT_CUR%OMEGA, N_MIX_PAW, RHOLM, RHOLM_LAST, &
           RMST,RMSC,RMSP,WEIGHT,IERRBR)
        MIX%LRESET=.FALSE.
      ELSE
!  simple mixing ... :
        RMST=0
        CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT,CHTOTL, &
             N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
      ENDIF
      DWRITE0 "mixing is ok"

      CALL VTIME(TV,TC)
      IF (IO%NWRITE>=2.AND.IO%IU6>=0)WRITE(TIU6,2300)'MIXING',TV-TV0,TC-TC0
!---- ENDIF (MIX%IMIX/=0)     end of mixing
      ENDIF
!-----ENDIF (.NOT.INFO%LCHCON)   end of charge update
      ENDIF

      CALL VTIME(TV0,TC0)
      IF (W%OVER_BAND) THEN
         CALL REDIS_PW_ALL(WDES, W)
         CALL VTIME(TV,TC)
         IF (IO%NWRITE>=2.AND.IO%IU6>=0)WRITE(TIU6,2300)'REDIS',TV-TV0,TC-TC0
      ENDIF
!=======================================================================
! total time used for this step
!=======================================================================
     CALL VTIME(TVPUL,TCPUL)
     io_begin
      IF (IO%NWRITE>=2.AND.IO%IU6>=0) THEN
        WRITE(TIU6,2310)'LOOP',TVPUL-TVPUL0,TCPUL-TCPUL0
      ELSE
        WRITE(TIU6,2300)'LOOP',TVPUL-TVPUL0,TCPUL-TCPUL0 
      ENDIF
     io_end
      CALL XML_TIMING(TVPUL-TVPUL0,TCPUL-TCPUL0,name="total")  

      TVPUL0=TVPUL
      TCPUL0=TCPUL

 2300 FORMAT(2X,A8,':  VPU time',F8.2,': CPU time',F8.2)
 2310 FORMAT(2X,'  ------------------------------------------'/ &
     &       2X,A8,':  VPU time',F8.2,': CPU time',F8.2)

!=======================================================================
!  important write statements
!=======================================================================

 2440 FORMAT(/' eigenvalue-minimisations  :',I6,/ &
     &       ' total energy-change (2. order) :',E14.7,'  (',E14.7,')')
 2441 FORMAT(/ &
     &       ' Broyden mixing:'/ &
     &       '  rms(total) =',E12.5,'    rms(broyden)=',E12.5,/ &
     &       '  rms(prec ) =',E12.5/ &
     &       '  weight for this iteration ',F10.2)

 2442 FORMAT(/' eigenvalues of (default mixing * dielectric matrix)' / &
             '  average eigenvalue GAMMA= ',F8.4,/ (10F8.4))

 200  FORMAT(' number of electron ',F12.7,' magnetization ',3F12.7)
 201  FORMAT(' augmentation part  ',F12.7,' magnetization ',3F12.7)

      NORDER=0
      IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR

      DO I=1,WDES%NCDIJ
         RHOTOT(I)=RHO0(GRIDC, CHTOT(1,I))
         RHOAUG(I)=RHOTOT(I)-RHO0(GRID_SOFT, CHDEN(1,I))
      END DO

      io_begin

    ! iteration counts
      WRITE(TIU6,2440) ICOUEV,DESUM(N),DESUM1

    ! charge density
      WRITE(TIU6,200) RHOTOT
      IF (INFO%LOVERL) THEN
        WRITE(TIU6,201) RHOAUG
      ENDIF
    ! dipol moment
      IF (DIP%LCOR_DIP) CALL WRITE_DIP(IO%IU6, DIP)

    ! mixing
      IF ( INFO%LMIX .AND. MIX%IMIX==4 ) THEN
        IF (IERRBR/=0) THEN
          IF (IO%IU0>=0) &
          WRITE(TIU0,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                      'mixing'' now and reset mixing at next step!'
          WRITE(TIU6,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                       'mixing'' now and reset mixing at next step!'
        ENDIF

        IF (IO%NWRITE>=2 .OR. NSTEP==1) THEN
          WRITE(TIU6,2441) RMST,RMSC,RMSP,WEIGHT
          IF (ABS(RMST-RMSC)/RMST> 0.1_q) THEN
            WRITE(TIU6,*) ' WARNING: grid for Broyden might be to small'
          ENDIF
        ENDIF
        IF (IO%IU0>=0) &
        WRITE(TIU0,308) RMST
        WRITE(17,308) RMST
   308  FORMAT('   ',E10.3)
        IF (MIX%NEIG > 0) THEN
           WRITE(TIU6,2442) MIX%AMEAN,MIX%EIGENVAL(1:MIX%NEIG)
        ENDIF
      ELSE IF (INFO%LMIX) THEN
        IF (IO%IU0>=0) &
        WRITE(TIU0,308) RMST
        WRITE(17,308) RMST
      ELSE
        IF (IO%IU0>=0) &
        WRITE(TIU0,*)
        WRITE(17 ,*)
      ENDIF
 io1: IF (IO%NWRITE>=2 .OR. (NSTEP==1)) THEN
    ! energy
      WRITE(TIU6,7240) E%PSCENC,E%TEWEN,E%DENC,E%XCENC,E%PAWPS,E%PAWAE, &
                      E%EENTROPY,E%EBANDSTR,INFO%EALLAT,TOTEN, &
                      TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)
    ! Aspherical and METAGGA energies  
      IF (INFO%LABORT .AND. INFO%LASPH) THEN
         WRITE(TIU6,72611) E%PAWPSG,E%PAWAEG,E%PAWPSAS,E%PAWAEAS,E%PAWCORE, &
              E%TOTENASPH,E%TOTENASPH-E%EENTROPY, &
              E%TOTENASPH-E%EENTROPY/(2+NORDER)
      ENDIF
      IF (INFO%LABORT .AND. INFO%LMETAGGA) &
            WRITE(TIU6,72612) E%EXCG,E%PAWPSAS,E%PAWAEAS,E%PAWCORE,E%EXCM,E%PAWPSM,E%PAWAEM,E%PAWCOREM, &
            E%TOTENMGGA,E%TOTENMGGA-E%EENTROPY, &
            E%TOTENMGGA-E%EENTROPY/(2+NORDER)

7240  FORMAT(/ &
              ' Free energy of the ion-electron system (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  alpha Z        PSCENC = ',F18.8/ &
     &        '  Ewald energy   TEWEN  = ',F18.8/ &
     &        '  -1/2 Hartree   DENC   = ',F18.8/ &
     &        '  -V(xc)+E(xc)   XCENC  = ',F18.8/ &
     &        '  PAW double counting   = ',2F18.8/ &
     &        '  entropy T*S    EENTRO = ',F18.8/ &
     &        '  eigenvalues    EBANDS = ',F18.8/ &
     &        '  atomic energy  EATOM  = ',F18.8/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free energy    TOTEN  = ',F18.8,' eV'// &
     &        '  energy without entropy =',F18.8, &
     &        '  energy(sigma->0) =',F18.8)
72611 FORMAT(//&
     &        '  ASPHERICAL CONTRIBUTION TO EXCH AND CORRELATION IN SPHERES (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  standard PAW   PS : AE= ',2F18.6/ &
     &        '  Aspheric PAW   PS : AE= ',2F18.6/ &
     &        '  core xc             AE= ',1F18.6/ &
     &        '  ---------------------------------------------------'/ &
     &        '  Aspherical result:'/ &
     &        '  free  energy   TOTEN  = ',F18.6,' eV'// &
     &        '  energy  without entropy=',F18.6, &
     &        '  energy(sigma->0) =',F18.6)
72612 FORMAT(//&
     &        '  METAGGA EXCHANGE AND CORRELATION (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  LDA+GGA E(xc)  EXCG   = ',F18.6/ &
     &        '  LDA+GGA PAW    PS : AE= ',2F18.6/ &
     &        '  core xc             AE= ',1F18.6/ &
     &        '  metaGGA E(xc)  EXCM   = ',F18.6/ &
     &        '  metaGGA PAW    PS : AE= ',2F18.6/ &
     &        '  metaGGA core xc     AE= ',1F18.6/ &
     &        '  ---------------------------------------------------'/ &
     &        '  METAGGA result:'/ &
     &        '  free  energy   TOTEN  = ',F18.6,' eV'// &
     &        '  energy  without entropy=',F18.6, &
     &        '  energy(sigma->0) =',F16.6)
   ELSE io1
      WRITE(TIU6,7242) TOTEN,TOTEN-E%EENTROPY
 7242 FORMAT(/'  free energy = ',E20.12, &
     &        '  energy without entropy= ',E20.12)

   ENDIF io1

      IF (IO%LOPEN) CALL WFORCE(IO%IU6)
      IF (IO%LOPEN) CALL WFORCE(17)
      WRITE(TIU6,130)
      io_end
!=======================================================================
!  perform some additional write statments if required
!=======================================================================
!-----Eigenvalues and weights
      IF (((NSTEP==1 .OR.NSTEP==DYN%NSW).AND.INFO%LABORT).OR. &
     &     (IO%NWRITE>=1 .AND.INFO%LABORT).OR.IO%NWRITE>=3) THEN

        ! calculate the core level shifts
      IF (INFO%LOVERL) THEN        
        CALL CL_SHIFT_PW( GRIDC, LATT_CUR, IRDMAX,  &
           T_INFO, P, WDES%NCDIJ, CVTOT, INFO%ENAUG, IO%IU6)
      ELSE
        do_io WRITE(*,*) " **** core level shifts not calculated ****"
      ENDIF
 
      io_begin
        CALL RHOAT0(P,T_INFO, BETATO,LATT_CUR%OMEGA)

        WRITE(TIU6,2202) EFERMI,REAL( E%CVZERO ,KIND=q) ,E%PSCENC/INFO%NELECT+BETATO
 2202   FORMAT(' E-fermi : ', F8.4,'     XC(G=0): ',F8.4, &
     &         '     alpha+bet :',F8.4/ &
     &         '   add alpha+bet to get absolut eigen values')

        DO ISP=1,WDES%ISPIN
        IF (WDES%ISPIN==2) WRITE(TIU6,'(/A,I1)') ' spin component ',ISP
        DO NN=1,KPOINTS%NKPTS
        WRITE(TIU6,2201)NN,WDES%VKPT(1,NN),WDES%VKPT(2,NN),WDES%VKPT(3,NN), &
     &      (I,REAL( W%CELTOT(I,NN,ISP) ,KIND=q) ,W%FERTOT(I,NN,ISP)*WDES%RSPIN,I=1,WDES%NB_TOT)
        ENDDO
        ENDDO

 2201   FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
     &         '  band No.  band energies     occupation '/ &
     &           (3X,I4,3X,F10.4,3X,F10.5))

!-----Charge-density along one line
        WRITE(TIU6,130)
        DO I=1,WDES%NCDIJ
           WRITE(TIU6,*)'soft charge-density along one line, spin component',I
           WRITE(TIU6,'(10(6X,I4))') (II,II=0,9)
           CALL WRT_RC_LINE(IO%IU6,GRID, CHDEN(1,I))
           IF (INFO%LOVERL) THEN
              WRITE(TIU6,*)'total charge-density along one line'
              CALL WRT_RC_LINE(IO%IU6,GRIDC, CHTOT(1,I))
           ENDIF
           WRITE(TIU6,*)
        ENDDO
!-----pseudopotential strength and augmentation charge
        DO I=1,WDES%NCDIJ
           WRITE(TIU6,*) 'pseudopotential strength for first ion, spin component:',I
           DO LP=1,P(1)%LMMAX
              WRITE(TIU6,'(16(F7.3,1X))') &
     &             (CDIJ(L,LP,1,I),L=1,MIN(8,P(1)%LMMAX))
!     &             (REAL(CDIJ(L,LP,1,I),q),L=1,MIN(16,P(1)%LMMAX))
           ENDDO
        ENDDO

        IF (INFO%LOVERL) THEN
        DO I=1,WDES%NCDIJ
           WRITE(TIU6,*) 'total augmentation occupancy for first ion, spin component:',I
           DO LP=1,P(1)%LMMAX
              WRITE(TIU6,'(16(F7.3,1X))') &
     &             (REAL(CRHODE(L,LP,1,I),q),L=1,MIN(16,P(1)%LMMAX))
           ENDDO
        ENDDO
        ENDIF
      io_end

      ENDIF
!=======================================================================
!  xml related output
!=======================================================================
      CALL XML_TAG("energy")
      IF (INFO%LABORT .OR. N==1) THEN
         CALL XML_TAG_REAL("alphaZ",E%PSCENC)
         CALL XML_TAG_REAL("ewald", E%TEWEN)
         CALL XML_TAG_REAL("hartreedc",E%DENC)
         CALL XML_TAG_REAL("XCdc",E%XCENC)
         CALL XML_TAG_REAL("pawpsdc",E%PAWPS)
         CALL XML_TAG_REAL("pawaedc",E%PAWAE)
         CALL XML_TAG_REAL("eentropy",E%EENTROPY)
         CALL XML_TAG_REAL("bandstr",E%EBANDSTR)
         CALL XML_TAG_REAL("atom",INFO%EALLAT)
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ELSE
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ENDIF
      CALL XML_CLOSE_TAG

      IF (INFO%LABORT .AND. INFO%LASPH) THEN
         CALL XML_TAG("aspherical")
         CALL XML_ENERGY(E%TOTENASPH,E%TOTENASPH-E%EENTROPY, E%TOTENASPH-E%EENTROPY/(2+NORDER))
         CALL XML_CLOSE_TAG
      ENDIF
      IF (INFO%LABORT .AND. INFO%LMETAGGA) THEN
         CALL XML_TAG("metagga")
         CALL XML_TAG_REAL("e_fr_energy",E%TOTENMGGA)
         CALL XML_TAG_REAL("e_wo_entrp", E%TOTENMGGA-E%EENTROPY)
         CALL XML_TAG_REAL("e_0_energy", E%TOTENMGGA-E%EENTROPY/(2+NORDER))
         CALL XML_CLOSE_TAG
      ENDIF

      CALL XML_CLOSE_TAG("scstep")
!======================== end of loop ENDLSC ===========================
! This is the end of the selfconsistent calculation loop
!=======================================================================
      IF (INFO%LABORT) THEN
         io_begin
        WRITE(TIU6,131)
 131    FORMAT (5X, //, &
     &  '------------------------ aborting loop because EDIFF', &
     &  ' is reached ----------------------------------------'//)
         io_end
        EXIT electron
      ENDIF
      INFO%LSOFT=.FALSE.
#if ! defined(noSTOPCAR) && ! defined(F90_T3D)
      CALL RDATAB(IO%LOPEN,'STOPCAR',99,'LABORT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,INFO%LSOFT,CHARAC,NCOUNT,1,IERR)
#endif
      IF (INFO%LSOFT) THEN
        io_begin
        IF (IO%IU0>=0) &
        WRITE(TIU0,*) 'hard stop encountered!  aborting job ...'
        WRITE(TIU6,13131)
13131   FORMAT (5X, //, &
     &  '------------------------ aborting loop because hard', &
     &  ' stop was set ---------------------------------------'//)
        io_end
        EXIT electron
      ENDIF
      TOTENL=TOTEN

      ENDDO electron

!
! calculate dipol corrections now
!
      IF ( DIP%IDIPCO >0 ) THEN
         IF (.NOT. DIP%LCOR_DIP) THEN
            CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, DIP, &
            CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT )
            
            CALL WRITE_DIP(IO%IU6, DIP)
            IF (IO%IU6>0) THEN
               WRITE(TIU6,*)
               WRITE(TIU6,*) &
               " *************** adding dipol energy to TOTEN NOW **************** "
            ENDIF
         TOTEN=TOTEN+ DIP%EDIPOL
         ENDIF
      ENDIF

      DWRITE0 'electron left'

      RETURN
      END SUBROUTINE
