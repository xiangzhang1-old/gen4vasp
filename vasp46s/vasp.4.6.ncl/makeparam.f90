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





!****************** PROGRAM VASP  Version 4.6.20 (f90)*****************
! RCS:  $Id: main.F,v 1.18 2003/06/27 13:22:18 kresse Exp kresse $
! Vienna Ab initio total energy and Molecular-dynamics Program
!            written  by   Kresse Georg
!                     and  Juergen Furthmueller
! Georg Kresse                       email: Georg.Kresse@univie.ac.at
! Juergen Furthmueller               email: furth@ifto.physik.uni-jena.de
! Institut fuer Materialphysik         voice: +43-1-4277-51402
! Uni Wien, Sensengasse 8/12           fax:   +43-1-4277-9514 (or 9513)
! A-1090 Wien, AUSTRIA                 http://cms.mpi.univie.ac.at/kresse
!
! This program comes without any waranty.
! No part of this program must be distributed, modified, or supplied
! to any other person for any reason whatsoever
! without prior written permission of the Institut of Theoretical Physics
! TU Vienna, Austria.
!
! This program performs total energy calculations using
! a selfconsistency cylce (i.e. mixing + iterative matrix diagonal.)
! most of the algorithms implemented are described in
! G. Kresse and J. Furthmueller
!  Efficiency of ab--initio total energy calculations for
!   metals and semiconductors using a plane--wave basis set
!  Comput. Mat. Sci. 6,  15-50 (1996)
! G. Kresse and J. Furthmueller
!  Efficient iterative schemes for ab--initio total energy
!   calculations using a plane--wave basis set
!   Phys. Rev. B 54, 11169 (1996)
!
! The iterative matrix diagonalization is based
! a) on the conjugated gradient eigenvalue minimisation proposed by
!  D.M. Bylander, L. Kleinmann, S. Lee, Phys Rev. B 42, 1394 (1990)
! and is a variant of an algorithm proposed by
!  M.P. Teter, M.C. Payne and D.C. Allan, Phys. Rev. B 40,12255 (1989)
!  T.A. Arias, M.C. Payne, J.D. Joannopoulos, Phys Rev. B 45,1538(1992)
! b) or the residual vector minimization method (RMM-DIIS) proposed by
!  P. Pulay,  Chem. Phys. Lett. 73, 393 (1980).
!  D. M. Wood and A. Zunger, J. Phys. A, 1343 (1985)
! For the mixing a Broyden/Pulay like method is used (see for instance):
!  D. Johnson, Phys. Rev. B 38, 12807 (1988)
!
! The program works with normconserving PP, 
! generalised ultrasoft-PP (Vanderbilt-PP Vanderbilt Phys Rev B 40,  
! 12255 (1989)) and PAW (P.E. Bloechl, Phys. Rev. B{\bf 50}, 17953 (1994))
! datasets. Partial core corrections can be handled
! Spin and GGA are implemented
!
! The units used in the programs are electron-volts and angstroms.
! The unit cell is arbitrary, and arbitrary species of ions are handled.
! A full featured symmetry-code is included, and calculation of
! Monkhorst-Pack special-points is possible (tetrahedron method can be
! used as well).
!
! The integretion of the ionic movements is performed using a predictor-
! corrector (Nose) dynamics, a  conjugate gradient techniques,
! or a Broyden/Pulay like technique (RMM-DIIS)
!
! The original version was written by  M.C. Payne
! at Professor J. Joannopoulos research  group at the MIT
! (3000 lines, excluding FFT, July 1989)
! The program was completely rewritten and vasply extended by
! Kresse Georg (gK) and Juergen Furthmueller. Currently the
! code has about 60000 source lines
! Some of the additions made by gK:
!  nose-dynamic,  predictor-corrector scheme for
!  ionic movements, conjugate-gradient scheme,
!  sub-space alignment, non-local pseudopotentials in real and
!  reciprocal space, generalized Vanderbild pseudopotentials
!  general cellshapes, arbitrary species of ions, stresstensor
!  residual minimization,
!  all bands simultaneous update, US-PP, PAW method
! Juergen Furthmueller(jF) wrote the symmetry code, the special-kpoint
! generation code, Broyden/Pulay mixing (with support of gK),
! and implemented the first GGA version
!
!** The following parts have been taken from other programs
!   please contact the authors of these programs
!   before using them
! - Tetrahedron method (original author unknown, parts have been written
!                     by Peter Bloechl probably)
!
! please refer to the README file to learn about new features
! notes on singe-precision:
! USAGE NOT RECOMMENDED DUE TO FINITE DIFFERENCES IN SOME
! FORCES-ROUTINES
! (except for native 64-bit-REAL machines like the CRAYs ...)
!**********************************************************************

      PROGRAM VAMP
      USE prec

      USE charge
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE paw
      USE pot
      USE force
      USE fileio
      USE nonl
      USE nonlr
      USE rmm_diis
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
      USE main_mpi
      USE chain
      USE pardens
      USE finite_differences
      USE LDAPLUSU_MODULE
      USE cl
!-MM- Added to accomodate constrained moments etc
      USE Constrained_M_modular
      USE writer
!-MM- end of additions
      USE sym_prec
      USE elpol
      USE mdipol
      USE compat_gga
      USE vaspxml

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

!=======================================================================
!  a small set of parameters might be set here
!  but this is really rarely necessary :->
!=======================================================================
!-----hard limits for k-point generation package
!     NTETD  is the maximum number of tetrahedra which can be
!            stored when using the tetrahedron integration and
!     IKPTD  is the maximum number of k-mesh points in each
!            direction that can be stored in the 'connection'
!            tables for the k-points used for the symmetry
!            reduction of the tetrahedron tiling. Finally
      INTEGER, PARAMETER :: NKDIMD=10000,NTETD=90000,IKPTD=45
!----I/O-related things (adapt on installation or for special purposes)
!     IU6    overall output ('console protocol'/'OUTCAR' I/O-unit)
!     IU0    very important output ('standard [error] output I/O-unit')
!     IU5    input-parameters ('standard input'/INCAR I/O-unit)
!     ICMPLX size of complex items (in bytes/complex item)
!     MRECL  maximum record length for direct access files
!            (if no restictions set 0 or very large number)
      INTEGER,PARAMETER :: ICMPLX=16,MRECL=10000000

!=======================================================================
!  structures
!=======================================================================
      TYPE (potcar),ALLOCATABLE :: P(:)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavefun)     W_F        ! wavefunction for all bands simultaneous
      TYPE (wavefun)     W_G        ! same as above
      TYPE (wavefun)     WUP
      TYPE (wavefun)     WDW
      TYPE (wavefun)     WTMP       ! temporary
      TYPE (wavespin)    WTMP_SPIN  ! temporary
      TYPE (latt)        LATT_CUR
      TYPE (latt)        LATT_INI
      TYPE (exctable)    EXCTAB
      TYPE (type_info)   T_INFO
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID       ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! very find grid temporarily used in us.F
      TYPE (grid_3d)     GRIDHF     ! grid used to present the potential
      TYPE (grid_3d)     GRIDB      ! Broyden grid
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRID_SOFT and GRIDC
      TYPE( prediction)  PRED
      TYPE (dipol)       DIP
      TYPE (smear_struct) SMEAR_LOOP
      TYPE (paco_struct) PACO
      TYPE (energy)      E,E2

      INTEGER :: NGX,NGY,NGZ,NGXC,NGYC,NGZC
      INTEGER :: NRPLWV,NBLK,LDIM,LMDIM,LDIM2,LMYDIM
      INTEGER :: IRMAX,IRDMAX,ISPIND
      INTEGER :: NPLWV,MPLWV,NPLWVC,MPLWVC,NTYPD,NIOND,NIONPD,NTYPPD
      INTEGER :: NBANDS
      INTEGER :: NEDOS
!R.S
      LOGICAL junk
      integer tiu6, tiu0, tiuvtot

!=======================================================================
!  begin array dimensions ...
!=======================================================================
!-----charge-density in real reciprocal space, partial core charge
      COMPLEX(q),ALLOCATABLE:: CHTOT(:,:)    ! charge-density in real / reciprocal space
      COMPLEX(q),ALLOCATABLE:: CHTOTL(:,:)   ! old charge-density
      COMPLEX(q)     ,ALLOCATABLE:: DENCOR(:)     ! partial core
      COMPLEX(q),ALLOCATABLE:: CVTOT(:,:)    ! local potential
      COMPLEX(q),ALLOCATABLE:: CSTRF(:,:)    ! structure-factor
!-----METAGGA logical flag - obsolete: now in INFO
!      LOGICAL LMETAGGA
!-----non-local pseudopotential parameters
      COMPLEX(q),ALLOCATABLE:: CDIJ(:,:,:,:) ! strength of PP
      COMPLEX(q),ALLOCATABLE:: CQIJ(:,:,:,:) ! overlap of PP
      COMPLEX(q),ALLOCATABLE:: CRHODE(:,:,:,:) ! augmentation occupancies
!-----elements required for mixing in PAW method
      REAL(q)   ,ALLOCATABLE::   RHOLM(:,:),RHOLM_LAST(:,:)
!-----charge-density and potential on small grid
      COMPLEX(q),ALLOCATABLE:: CHDEN(:,:)    ! soft part of charge density
      COMPLEX(q)  ,ALLOCATABLE:: SV(:,:)       ! soft part of local potential
!-----description how to go from (1._q,0._q) grid to the second grid
!-----density of states
      REAL(q)   ,ALLOCATABLE::  DOS(:,:),DOSI(:,:)
      REAL(q)   ,ALLOCATABLE::  DDOS(:,:),DDOSI(:,:)
!-----local l-projected wavefunction characters
      REAL(q)   ,ALLOCATABLE::   PAR(:,:,:,:,:),DOSPAR(:,:,:,:)
!  all-band-simultaneous-update arrays
      COMPLEX(q)   ,ALLOCATABLE::   CHF(:,:,:),CHAM(:,:,:)
!  optics stuff
      COMPLEX(q)   ,ALLOCATABLE::   NABIJ(:,:)

!-----remaining mainly work arrays
      COMPLEX(q), ALLOCATABLE,TARGET :: CWORK1(:),CWORK2(:),CWORK(:,:)
      TYPE (wavefun1)    W1            ! current wavefunction
      TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point

      COMPLEX(q), ALLOCATABLE  ::  CPROTM(:),CMAT(:,:)
!=======================================================================
!  a few fixed size (or small) arrays
!=======================================================================
!----- energy at each step
      REAL(q)      DESUM(500)
!-----Forces and stresses
      REAL(q)      TFORNL(3),TEIFOR(3),TEWIFO(3),THARFO(3),VTMP(3), &
     &          SIKEF(3,3),EISIF(3,3),DSIF(3,3),XCSIF(3,3), &
     &          PSCSIF(3,3),EWSIF(3,3),FNLSIF(3,3),AUGSIF(3,3), &
     &          TSIF(3,3),D2SIF(3,3),PARSIF(3,3)
!-----forces on ions
      REAL(q)   ,ALLOCATABLE::  EIFOR(:,:),EINL(:,:),EWIFOR(:,:), &
     &         HARFOR(:,:),TIFOR(:,:),PARFOR(:,:)
      REAL(q)  STM(5)
!-----Temporary data for tutorial messages ...
      INTEGER,PARAMETER :: NTUTOR=1000
      REAL(q)     RTUT(NTUTOR),RDUM
      INTEGER  ITUT(NTUTOR),IDUM
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!=======================================================================
!  end array dimensions ...
!=======================================================================
      INTEGER NTYP_PP      ! number of types on POTCAR file

      INTEGER I,J,N,NT,K
!---- used for creation of param.inc
      REAL(q)    WFACT,PSRMX,PSDMX
      REAL(q)    XCUTOF,YCUTOF,ZCUTOF

!---- timing information
      INTEGER MINPGF,MAJPGF,ISWPS,IOOPS,IVCSW,IERR
      REAL(q)    UTIME,STIME,ETIME,RSIZM,AVSIZ

      INTEGER NORDER   !   order of smearing
!---- timing of individual calls
      REAL(q)    TC,TV,TC0,TV0

!---- a few logical and string variables
      LOGICAL    LTMP,LSTOP2
      LOGICAL    LPAW           ! paw is used 
      LOGICAL    LPARD          ! partial band decomposed charge density
      LOGICAL    LREALLOCATE    ! reallocation of proj operators required
      LOGICAL    L_NO_US        ! no ultrasoft PP
      LOGICAL    LADDGRID       ! additional support grid


      LOGICAL    LBERRY         ! calculate electronic polarisation
      CHARACTER (LEN=40)  SZ
      CHARACTER (LEN=1)   CHARAC
      CHARACTER (LEN=5)   IDENTIFY
!-----parameters for sphpro.f
      INTEGER :: LDIMP,LMDIMP,LTRUNC=3
!=======================================================================
! All COMMON blocks
!=======================================================================
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
      COMMON /WAVCUT/ IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX

      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      REAL(q)  RHOTOT(4)
      INTEGER(8) IL,I1,I2_0,I3,I4
      CHARACTER (LEN=*),PARAMETER :: VASP='vasp.4.6.28 25Jul05 complex '

      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
!=======================================================================
!  initialise / set constants and parameters ...
!=======================================================================
      IO%LOPEN =.TRUE.  ! open all files with file names
      IO%IU0   = 6
      IO%IU6   = 8
      IO%IU5   = 5
!R.S
      tiu6 = IO%IU6
      tiu0 = IO%IU0

      IO%ICMPLX=ICMPLX
      IO%MRECL =MRECL
      PRED%ICMPLX=ICMPLX

      CALL TIMING(0,UTIME,STIME,ETIME,MINPGF,MAJPGF, &
     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) ETIME=0._q
! switch off kill
!     CALL sigtrp()

      NPAR=1
      IUXML_SET=20
      CALL START_XML( IUXML_SET, "vasprun.xml" )
!-----------------------------------------------------------------------
!  open Files
!-----------------------------------------------------------------------
      IF (IO%IU0>=0) WRITE(TIU0,*) VASP
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTPAR',STATUS='UNKNOWN')
      OPEN(UNIT=18,FILE=DIR_APP(1:DIR_LEN)//'CHGCAR',STATUS='UNKNOWN')

!R.S
      junk=.TRUE.
      INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',EXIST=junk)
      IO%LFOUND=junk

! first reopen with assumed (wrong) record length ICMPLX
      OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
                   FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=ICMPLX)
! the first record contains the record length, get it ...
      RDUM=0._q
      READ(12,REC=1,ERR=17421) RDUM,RISPIN ; IDUM=NINT(RDUM)
! in the worst case IDUM could contain completely useless data and useless is
! all <=0 or all >10000000 (since on output we use RECL=ICMPLX*MAX(NRPLWV,6)
! or RECL=(NB_TOT+2)*ICMPLX more than ten millions sounds not very realistic)
      IF ((IDUM<=0).OR.(IDUM>10000000)) IDUM=ICMPLX  ! -> error reading WAVECAR
      GOTO 17422
17421 CONTINUE
      IDUM=ICMPLX  ! -> error reading WAVECAR
17422 CONTINUE
      CLOSE(12)
! reopen with correct record length (clumsy all that, I know ...)
      OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
                   FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IDUM)





      IF (IO%IU6>=0) WRITE(IO%IU6,*) VASP
      CALL XML_GENERATOR



      CALL PARSE_GENERATOR_XML(VASP//" serial")

      CALL MY_DATE_AND_TIME(IO%IU6)
      CALL XML_CLOSE_TAG

      CALL WRT_DISTR(IO%IU6)


! unit for extrapolation of wavefunction
      PRED%IUDIR =21
! unit for broyden mixing
      MIX%IUBROY=23
! unit for total potential
      IO%IUVTOT=62

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'----------------------------------------- Iteration ', &
     &I4,'(',I4,')  ---------------------------------------'//)
!-----------------------------------------------------------------------
! read header of POSCAR file to get NTYPD, NTYPDD, NIOND and NIONPD
!-----------------------------------------------------------------------
      CALL RD_POSCAR_HEAD(LATT_CUR, T_INFO, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

      ALLOCATE(T_INFO%ATOMOM(3*NIOND),T_INFO%RWIGS(NTYPPD),T_INFO%ROPT(NTYPD),T_INFO%POMASS(NTYPD))

      IF (IO%IU6>=0) THEN
         WRITE(TIU6,130)
         WRITE(TIU6,*)'INCAR:'
      ENDIF
!  first scan of POSCAR to get LDIM, LMDIM, LDIM2 ...
      LDIM =11
      LDIM2=(LDIM*(LDIM+1))/2
      LMDIM=32

      ALLOCATE(P(NTYPD))
      T_INFO%POMASS=0
      T_INFO%RWIGS=0

!-----------------------------------------------------------------------
! read pseudopotentials
!-----------------------------------------------------------------------
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS, &
     &           IO%IU0,IO%IU6,-1,LPAW)

!-----------------------------------------------------------------------
! read INCAR
!-----------------------------------------------------------------------
      CALL XML_TAG("incar")

      CALL READER( &
          IO%IU5,IO%IU0,INFO%SZNAM1,INFO%ISTART,INFO%IALGO,MIX%IMIX,MIX%MAXMIX,MIX%MREMOVE, &
          MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
          MIX%WC,MIX%INIMIX,MIX%MIXPRE,IO%LFOUND,INFO%LDIAG,INFO%LREAL,IO%LREALD,IO%LPDENS, &
          DYN%IBRION,INFO%ICHARG,INFO%INIWAV,INFO%NELM,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,DYN%EDIFFG, &
          DYN%NSW,DYN%ISIF,PRED%IWAVPR,SYMM%ISYM,DYN%NBLOCK,DYN%KBLOCK,INFO%ENMAX,DYN%POTIM,DYN%TEBEG, &
          DYN%TEEND,DYN%NFREE, &
          PACO%NPACO,PACO%APACO,T_INFO%NTYP,NTYPD,DYN%SMASS,T_INFO%POMASS, &
          T_INFO%RWIGS,INFO%NELECT,INFO%NUP_DOWN,INFO%TIME,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%ISMEAR,DYN%PSTRESS,INFO%NDAV, &
          KPOINTS%SIGMA,KPOINTS%LTET,INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,IO%NWRITE,INFO%LCORR, &
          IO%IDIOT,T_INFO%NIONS,T_INFO%NTYPP,IO%LMUSIC,IO%LOPTICS,STM, &
          INFO%ISPIN,T_INFO%ATOMOM,NIOND,IO%LWAVE,IO%LCHARG,IO%LVTOT,INFO%SZPREC,INFO%SZGGA, &
          DIP%LCOR_DIP,DIP%IDIPCO,DIP%POSCEN,INFO%ENAUG,IO%LORBIT,IO%LELF,T_INFO%ROPT,INFO%ENINI, &
          NGX,NGY,NGZ,NGXC,NGYC,NGZC,NBANDS,NEDOS,NBLK,LATT_CUR, &
          LPLANE_WISE,LCOMPAT,LMAX_CALC,SET_LMAX_MIX_TO,WDES%NSIM,LFCI,LPARD,LPAW,LADDGRID,WDES%LCRITICAL_MEM, &
          WDES%LNONCOLLINEAR,WDES%LSORBIT,WDES%SAXIS,INFO%LMETAGGA, &
          WDES%LSPIRAL,WDES%LZEROZ,WDES%QSPIRAL, &
          INFO%LASPH,INFO%LSECVAR,INFO%IGGA2,INFO%SZGGA2)


       CALL LREAL_COMPAT_MODE(IO%IU5, IO%IU0, LCOMPAT)
       CALL GGA_COMPAT_MODE(IO%IU5, IO%IU0, LCOMPAT)

       IF (WDES%LNONCOLLINEAR) THEN
          INFO%ISPIN = 1
       ENDIF
! METAGGA not implemented for non collinear magnetism
       IF (WDES%LNONCOLLINEAR .AND. INFO%LMETAGGA) THEN
          WRITE(*,*) 'METAGGA for non collinear magnetism not supported.' 
          WRITE(*,*) 'exiting VASP; sorry for the inconveniences.'
          STOP
       ENDIF
!-MM- Spin spirals require LNONCOLLINEAR=.TRUE.
       IF (.NOT.WDES%LNONCOLLINEAR .AND. WDES%LSPIRAL) THEN
          WRITE(*,*) 'Spin spirals require LNONCOLLINEAR=.TRUE. '
          WRITE(*,*) 'exiting VASP; sorry dude!'
          STOP
       ENDIF
!-MM- end of addition

       IF (LCOMPAT) THEN
               CALL VTUTOR('W','VASP.4.4',RTUT,1, &
     &                  ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
               CALL VTUTOR('W','VASP.4.4',RTUT,1, &
     &                  ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
       ENDIF
! WRITE out an advice if some force dependent ionic algorithm and METAGGA
! or ASPH
       IF ((INFO%LMETAGGA .OR. INFO%LASPH) .AND. &
     &       (DYN%IBRION>0 .OR. (DYN%IBRION==0 .AND. DYN%SMASS/=-2))) THEN
          CALL VTUTOR('A','METAGGA and forces',RTUT,1, &
     &                  ITUT,1,CDUM,1,(/INFO%LASPH, INFO%LMETAGGA /),2, &
     &                  IO%IU0,IO%IDIOT)
       ENDIF
!-----------------------------------------------------------------------
! core level shift related items (parses INCAR)
!-----------------------------------------------------------------------
      CALL INIT_CL_SHIFT(IO%IU5,IO%IU0, T_INFO%NIONS, T_INFO%NTYP )

      CALL READER_ADD_ON(IO%IU5,IO%IU0,LBERRY,IGPAR,NPPSTR, &
            INFO%ICHARG,KPOINTS%ISMEAR,KPOINTS%SIGMA)

      ISPIND=INFO%ISPIN


      DYN%TEMP =DYN%TEBEG
      INFO%RSPIN=3-INFO%ISPIN

!-----------------------------------------------------------------------
! loop over different smearing parameters
!-----------------------------------------------------------------------

      SMEAR_LOOP%ISMCNT=0
      IF (KPOINTS%ISMEAR==-3) THEN
        IF(IO%IU6>=0)   WRITE(TIU6,7219)
 7219   FORMAT('   Loop over smearing-parameters in INCAR')
        CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SMEARINGS','=','#',';','F', &
     &            IDUM,SMEAR_LOOP%SMEARS(1),CDUM,LDUM,CHARAC,N,200,IERR)
        IF ((IERR/=0).OR.((IERR==0).AND. &
     &          ((N<2).OR.(N>200).OR.(MOD(N,2)/=0)))) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*)'Error reading item ''SMEARINGS'' from file INCAR.'
           STOP
        ENDIF
        SMEAR_LOOP%ISMCNT=N/2
        DYN%NSW   =SMEAR_LOOP%ISMCNT+1
        DYN%KBLOCK=DYN%NSW
        KPOINTS%LTET  =.TRUE.
        DYN%IBRION=-1
        KPOINTS%ISMEAR=-5
      ENDIF
!=======================================================================
!  now read in Pseudopotential
!=======================================================================
      LMDIM=0
      LDIM=0
      DO NT=1,NTYP_PP
        LMDIM=MAX(LMDIM,P(NT)%LMMAX)
        LDIM =MAX(LDIM ,P(NT)%LMAX)
      END DO
      CALL DEALLOC_PP(P,NTYP_PP)

      LDIM2=(LDIM*(LDIM+1))/2
      LMYDIM=9
! second scan with correct setting
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS, &
     &           IO%IU0,IO%IU6,IO%NWRITE,LPAW)
      CALL CL_MODIFY_PP( NTYP_PP, P, ENAUG )
! now check everything
      CALL POST_PSEUDO(NTYPD,NTYP_PP,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,P,INFO, &
     &        IO%LREALD,T_INFO%ROPT, IO%IDIOT,IO%IU6,IO%IU0,LMAX_CALC,L_NO_US)
      CALL LDIM_PSEUDO(IO%LORBIT, NTYPD, P, LDIMP, LMDIMP)
! setup PAW
      IF (.NOT.LPAW) IO%LOPTICS=.FALSE.

!-----------------------------------------------------------------------
! LDA+U initialisation (parses INCAR)
!-----------------------------------------------------------------------
      CALL LDAU_READER(T_INFO%NTYP,IO%IU5,IO%IU0)
      IF (USELDApU().OR.LCALC_ORBITAL_MOMENT()) &
     &   CALL INITIALIZE_LDAU(T_INFO%NIONS,T_INFO%NTYP,P,WDES%LNONCOLLINEAR,IO%IU0,IO%IDIOT)

      CALL SET_AUG(T_INFO%NTYP, P, IO%IU6, INFO%LEXCH, INFO%LEXCHG, LMAX_CALC, INFO%LMETAGGA, LCOMPAT)
!-----------------------------------------------------------------------
! optics initialisation (parses INCAR)
!-----------------------------------------------------------------------
      IF (IO%LOPTICS) CALL SET_NABIJ_AUG(P,T_INFO%NTYP)

!-----------------------------------------------------------------------
! exchange correlation table
!-----------------------------------------------------------------------
      IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
         CALL RD_EX(EXCTAB,2,INFO%LEXCH,IO%LEXCHF,IO%IU6,IO%IU0,IO%IDIOT)
      ELSE
         CALL RD_EX(EXCTAB,1,INFO%LEXCH,IO%LEXCHF,IO%IU6,IO%IU0,IO%IDIOT)
      ENDIF


!-----------------------------------------------------------------------
!  Read UNIT=15: POSCAR Startjob and Continuation-job
!-----------------------------------------------------------------------
      CALL RD_POSCAR(LATT_CUR, T_INFO, DYN, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, &
     &           IO%IU0,IO%IU6)

!-----------------------------------------------------------------------
! constrained moment reader (INCAR reader)
!-----------------------------------------------------------------------
      CALL CONSTRAINED_M_READER(T_INFO%NIONS,IO%IU0,IO%IU5)
      CALL WRITER_READER(IO%IU0,IO%IU5)
!      CALL WANNIER_READER(IO%IU0,IO%IU5,IO%IU6,IO%IDIOT)
      CALL FIELD_READER(DIP,IO%IU0,IO%IU5)

!-----------------------------------------------------------------------
! init all chains (INCAR reader)
!-----------------------------------------------------------------------
      CALL chain_init( T_INFO, IO)

!-----------------------------------------------------------------------
!xml finish copying parameters from INCAR to xml file
! no INCAR reading from here 
      CALL XML_CLOSE_TAG("incar")
!-----------------------------------------------------------------------

      CALL COUNT_DEGREES_OF_FREEDOM( T_INFO, NDEGREES_OF_FREEDOM, &
          IO%IU6, IO%IU0, DYN%IBRION)

!-----for static calculations or relaxation jobs DYN%VEL is meaningless
      IF (DYN%INIT == -1) THEN
        CALL INITIO(T_INFO%NIONS,T_INFO%LSDYN,NDEGREES_OF_FREEDOM, &
               T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP, &
               T_INFO%POMASS,DYN%POTIM, &
               DYN%POSION,DYN%VEL,T_INFO%LSFOR,LATT_CUR%A,LATT_CUR%B,DYN%INIT,IO%IU6)
        DYN%INIT=0
      ENDIF
      IF (DYN%IBRION/=0) THEN
          DYN%VEL=0._q
      ENDIF
      IF (IO%IU6>=0) THEN
         WRITE(TIU6,*)
         WRITE(TIU6,130)
      ENDIF

      IF ( T_INFO%LSDYN ) THEN
         CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
      ELSE
         CALL SYMVEL_WARNING( T_INFO%NIONS, T_INFO%NTYP, T_INFO%ITYP, &
             T_INFO%POMASS, DYN%VEL, IO%IU6, IO%IU0 )
      ENDIF
      CALL NEAREST_NEIGHBOAR(IO%IU6, IO%IU0, T_INFO, LATT_CUR, P%RWIGS)
!-----------------------------------------------------------------------
!  initialize the symmetry stuff
!-----------------------------------------------------------------------
      ALLOCATE(SYMM%ROTMAP(NIOND,48,NIOND),SYMM%TAU(NIOND,3), &
     &         SYMM%TAUROT(NIOND,3),SYMM%WRKROT(3*(NIOND+2)), &
     &         SYMM%PTRANS(NIOND+2,3),SYMM%INDROT(NIOND+2))
      IF (INFO%ISPIN==2) THEN
         ALLOCATE(SYMM%MAGROT(48,NIOND))
      ELSE
         ALLOCATE(SYMM%MAGROT(1,1))
      ENDIF
      ! break symmetry parallel to IGPAR
      IF (LBERRY) THEN
         LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)*(1+TINY*10)
         CALL LATTIC(LATT_CUR)
      ENDIF
! Rotate the initial magnetic moments to counter the clockwise
! rotation brought on by the spin spiral
      IF (WDES%LSPIRAL) CALL COUNTER_SPIRAL(WDES%QSPIRAL,T_INFO%NIONS,T_INFO%POSION,T_INFO%ATOMOM)

      IF (SYMM%ISYM>0) THEN
! Finite temperature allows no symmetry by definition ...
         NCDIJ=INFO%ISPIN
         IF (WDES%LNONCOLLINEAR) NCDIJ=4
         CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,T_INFO%LSFOR, &
                     T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND, &
                     SYMM%PTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                     SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
      ELSE
! ... so take nosymm!
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,NIOND,SYMM%PTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
      END IF

!=======================================================================
!  Read UNIT=14: KPOINTS
!  number of k-points and k-points in reciprocal lattice
!=======================================================================
      IF(IO%IU6>=0)  WRITE(TIU6,*)
      NKDIM=NKDIMD

      IF (LBERRY) THEN
         CALL RD_KPOINTS_BERRY(KPOINTS,NPPSTR,IGPAR, &
        &   LATT_CUR, NKDIM,IKPTD,NTETD, &
        &   SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
        &   IO%IU6,IO%IU0)
          IF (LBERRY) THEN
            LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)/(1+TINY*10)
            CALL LATTIC(LATT_CUR)
         ENDIF
      ELSE
         CALL RD_KPOINTS(KPOINTS,LATT_CUR, NKDIM,IKPTD,NTETD, &
           SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
           IO%IU6,IO%IU0)
      ENDIF





      NKDIM=KPOINTS%NKPTS
!=======================================================================
!  at this point we have enough information to
!  create a param.inc file
!=======================================================================
      XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
!
!  setup NGX, NGY, NGZ if required
!
! high precission do not allow for wrap around
      IF (INFO%SZPREC(1:1)=='h' .OR. INFO%SZPREC(1:1)=='a') THEN
        WFACT=4
      ELSE
! medium-low precission allow for small wrap around
        WFACT=3
      ENDIF
      GRID%NGPTAR(1)=XCUTOF*WFACT+0.5_q
      GRID%NGPTAR(2)=YCUTOF*WFACT+0.5_q
      GRID%NGPTAR(3)=ZCUTOF*WFACT+0.5_q
      IF (NGX /= -1)   GRID%NGPTAR(1)=  NGX
      IF (NGY /= -1)   GRID%NGPTAR(2)=  NGY
      IF (NGZ /= -1)   GRID%NGPTAR(3)=  NGZ
      CALL FFTCHK(GRID%NGPTAR)
!
!  setup NGXC, NGYC, NGZC if required
!
      IF (INFO%LOVERL) THEN
        IF (INFO%ENAUG==0) INFO%ENAUG=INFO%ENMAX*1.5_q
        IF (INFO%SZPREC(1:1)=='h') THEN
          WFACT=16._q/3._q
        ELSE IF (INFO%SZPREC(1:1)=='l') THEN
          WFACT=3
        ELSE
          WFACT=4
        ENDIF
        XCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
        YCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
        ZCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
        GRIDC%NGPTAR(1)=XCUTOF*WFACT
        GRIDC%NGPTAR(2)=YCUTOF*WFACT
        GRIDC%NGPTAR(3)=ZCUTOF*WFACT
        ! prec Accurate and Medium double grids
        IF (INFO%SZPREC(1:1)=='a' .OR. INFO%SZPREC(1:1)=='n') THEN
           GRIDC%NGPTAR(1)=GRID%NGPTAR(1)*2
           GRIDC%NGPTAR(2)=GRID%NGPTAR(2)*2
           GRIDC%NGPTAR(3)=GRID%NGPTAR(3)*2
        ENDIF
        IF (NGXC /= -1)  GRIDC%NGPTAR(1)=NGXC
        IF (NGYC /= -1)  GRIDC%NGPTAR(2)=NGYC
        IF (NGZC /= -1)  GRIDC%NGPTAR(3)=NGZC
        CALL FFTCHK(GRIDC%NGPTAR)
      ELSE
        GRIDC%NGPTAR(1)= 1
        GRIDC%NGPTAR(2)= 1
        GRIDC%NGPTAR(3)= 1
      ENDIF

      GRIDC%NGPTAR(1)=MAX(GRIDC%NGPTAR(1),GRID%NGPTAR(1))
      GRIDC%NGPTAR(2)=MAX(GRIDC%NGPTAR(2),GRID%NGPTAR(2))
      GRIDC%NGPTAR(3)=MAX(GRIDC%NGPTAR(3),GRID%NGPTAR(3))
      GRIDUS%NGPTAR=GRIDC%NGPTAR
      IF (LADDGRID) GRIDUS%NGPTAR=GRIDC%NGPTAR*2

      NGX = GRID %NGPTAR(1); NGY = GRID %NGPTAR(2); NGZ = GRID %NGPTAR(3)
      NGXC= GRIDC%NGPTAR(1); NGYC= GRIDC%NGPTAR(2); NGZC= GRIDC%NGPTAR(3)
!
      IF (NBANDS == -1) THEN
         IF (WDES%LNONCOLLINEAR)  THEN
             NMAG=MAX(SUM(T_INFO%ATOMOM(1:T_INFO%NIONS*3-2:3)), &
                      SUM(T_INFO%ATOMOM(2:T_INFO%NIONS*3-1:3)), &
                      SUM(T_INFO%ATOMOM(3:T_INFO%NIONS*3:3)))
         ELSE IF (INFO%ISPIN > 1) THEN
             NMAG=SUM(T_INFO%ATOMOM(1:T_INFO%NIONS))
         ELSE
             NMAG=0
         ENDIF
         NMAG = (NMAG+1)/2
         NBANDS=MAX(NINT(INFO%NELECT+2)/2+MAX(T_INFO%NIONS/2,3),INT(0.6*INFO%NELECT))+NMAG
         IF (WDES%LNONCOLLINEAR) NBANDS = NBANDS*2
      ENDIF
!rS    IF (NBANDS == -1) NBANDS=0.6*INFO%NELECT + 4

      IF (INFO%EBREAK == -1) INFO%EBREAK=0.25_q*MIN(INFO%EDIFF,ABS(DYN%EDIFFG)/10)/NBANDS

      INFO%NBANDTOT=((NBANDS+NPAR-1)/NPAR)*NPAR


      IF ((.NOT. WDES%LNONCOLLINEAR) .and. INFO%NELECT>REAL(INFO%NBANDTOT*2,KIND=q)) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'ERROR: Number of bands NBANDS too small to hold', &
                           ' electrons',INFO%NELECT,INFO%NBANDTOT*2
         STOP
      ELSEIF((WDES%LNONCOLLINEAR) .and. ((INFO%NELECT*2)>REAL(INFO%NBANDTOT*2,KIND=q))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'ERROR: Number of bands NBANDS too small to hold', &
                           ' electrons',INFO%NELECT,INFO%NBANDTOT
         STOP
      ENDIF

      NRPLWV=4*PI*SQRT(INFO%ENMAX /RYTOEV)**3/3* &
     &     LATT_CUR%OMEGA/AUTOA**3/(2*PI)**3*1.1_q+50







      IF (NBLK==-1) NBLK=MIN(256,MAX(32,(NRPLWV/320)*32))


      PSRMX=0
      PSDMX=0
      DO NT=1,T_INFO%NTYP
        PSRMX=MAX(PSRMX,P(NT)%PSRMAX)
        PSDMX=MAX(PSDMX,P(NT)%PSDMAX)
      ENDDO
      IF (INFO%LREAL) THEN
       IRMAX=4*PI*PSRMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRID%NGPTAR(1)*GRID%NGPTAR(2)*GRID%NGPTAR(3)))+50
      ELSE
       IRMAX=1
      ENDIF
      IRDMAX=1
      IF (INFO%LOVERL) THEN
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDC%NGPTAR(1)*GRIDC%NGPTAR(2)*GRIDC%NGPTAR(3)))+200
      ENDIF

       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDUS%NGPTAR(1)*GRIDUS%NGPTAR(2)*GRIDUS%NGPTAR(3)))+200


      NPLWV =NGX *NGY *NGZ;
      MPLWV =NGX *NGY *NGZ
      NPLWVC=NGXC*NGYC*NGZC;
      MPLWVC=NGXC*NGYC*NGZC


      ISP =INFO%ISPIN
      I0=6    ! CHTOT, CHTOTL, DENCOR, CVTOT + 2 work arrays
      IF (INFO%LEXCHG> 0)   I0=I0+4 ! 4 more for GGA
      IF (PRED%IWAVPR>=12 ) I0=I0+2 ! 2 more for prediction without file
      I1=1_8*(I0*ISP+T_INFO%NTYPD)*MPLWVC*IO%ICMPLX



      I2_0=1_8*(NRPLWV+T_INFO%NIOND*LMDIM)*NBANDS*KPOINTS%NKPTS*IO%ICMPLX*ISP

      I2=1
      IF (INFO%IALGO>=50)   I2=I2+2
      IF (PRED%IWAVPR>=12 ) I2=I2+2
      I2_0=I2_0*I2
      IF (INFO%LREAL) THEN
        I3= 1_8*IRMAX*(LMDIM+2)*T_INFO%NIOND*IO%ICMPLX/2
        I4= 0
      ELSE
        I3= 0
        I4= 1_8*NRPLWV*T_INFO%NIONS*IO%ICMPLX+ &
            NRPLWV*LMDIM*T_INFO%NTYPD*KPOINTS%NKPTS*IO%ICMPLX/2
      ENDIF
      IF (IO%IU0>=0) &
      WRITE(TIU0,20) INFO%ENMAX,INFO%ENAUG,KPOINTS%NKPTS,NBANDS,MPLWVC,NRPLWV,IRMAX,IRDMAX,I0, &
          I1,I2,I2_0,I3,I4,(I1+I2_0+I3+I4)/1000000._q,(I1+I2_0+I3+I4)

   20 FORMAT(' memory requirements are:',// &
             ' energy cutoff:                   ENCUT     ',F10.2,/ &
             ' energy cutoff (augmentation):    ENAUG     ',F10.2,/ &
             ' memory requirements are:',// &
             ' number of kpoints                NKPTS     ',I12,/  &
             ' number of bands                  NBANDS    ',I12,/  &
             ' number of grid points            MPLWV     ',I12,/  &
             ' number of plane waves            NRPLWV    ',I12,/  &
             ' real space IRMAX   ',I12,     '  IRDMAX    ',I12,/ &
             ' arrays on large grid (',I2,' + NTYP)*MPLWVM*16 ',I12,/  &
             ' ',I1,' sets of wavefunctions                    ',I12,/  &
             ' projectors in real space                   ',I12,/  &
             ' projectors in reciprocal space             ',I12,/  &
             '---------------------------------------------------------',/ &
             ' gives a total  of ',F8.1,' Mbytes       or ',I12,/ &
             ' add 5-10 Mbytes for executable (depends on OS)',//)


      IF (IO%IU0>=0) &
      WRITE(TIU0,10) GRID%NGPTAR(1),GRID%NGPTAR(2),GRID%NGPTAR(3), &
     &            GRIDC%NGPTAR(1),GRIDC%NGPTAR(2),GRIDC%NGPTAR(3), &
     &            NBANDS,NBLK

   10 FORMAT('  NGX =',I3,'; NGY =',I3,'; NGZ =',I3,/ &
     & '  NGXF=',I3,'; NGYF=',I3,'; NGZF=',I3,/ &
     & '  NBANDS=',I4,'; NBLK=',I4/)
      
      CALL STOP_XML

      END
