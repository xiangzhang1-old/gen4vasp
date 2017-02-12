#include "symbol.inc"
      SUBROUTINE READER &
     &       (IU5,IU0,SZNAM1,ISTART,IALGO,IMIX,MAXMIX,MREMOVE, &
     &        AMIX,BMIX,AMIX_MAG,BMIX_MAG,AMIN, &
     &        WC,INIMIX,MIXPRE,LFOUND,LDIAG,LREAL,LREALD, &
     &        LPDENS,IBRION,ICHARG,INIWAV,NELM,NELMIN,NELMDL,EDIFF, &
     &        EDIFFG,NSW,ISIF,IWAVPR,ISYM,NBLOCK,KBLOCK,ENMAX,POTIM, &
     &        TEBEG,TEEND,NFREE, &
     &        NPACO,APACO,NTYPIN,NTYPD,SMASS,POMASS, &
     &        RWIGS,NELECT,NUP_DOWN,TIME,EMIN,EMAX,ISMEAR,PSTRESS,NDAV, &
     &        SIGMA,LTET,WEIMIN,EBREAK,DEPER,NWRITE,LCORR, &
     &        IDIOT,NIONS,NTYPP,lmusic,LOPTICS,STM, &
     &        ISPIN,ATOMOM,NIOND,LWAVE,LCHARG,LVTOT,SZPREC,SZGGA, &
     &        LCOR_DIP,IDIPCO,POSCEN,ENAUG,LORBIT,LELF,ROPT,ENINI, &
     &        NGX,NGY,NGZ,NGXF,NGYF,NGZF,NBANDS,NEDOS,NBLK,LATT_CUR, &
     &        LPLANE_WISE,LCOMPAT,LMAX_CALC,LMAX_MIX,NSIM,LFCI,LPARD,LPAW,LADDGRID, &
     &        LCRITICAL_MEM,LNONCOLLINEAR,LSORBIT,SAXIS,LMETAGGA, &
     &        LSPIRAL,LZEROZ,QSPIRAL, &
     &        LASPH,LSECVAR,IGGA2,SZGGA2)

      USE prec
      USE sym_prec
      USE ini
      USE lattice
      USE scala
      USE wave_mpi
      USE constant
      USE pseudo   ! for subroutine EXTYP
      USE vaspxml
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (latt)        LATT_CUR
#ifdef LONGCHAR
      CHARACTER*32767  INPLIN
#else
      CHARACTER*255  INPLIN
#endif
      CHARACTER*1    CHARAC
      CHARACTER*40   SZNAM1
      CHARACTER*80   SZNAM
      CHARACTER*6    SZPREC
      CHARACTER*2    SZGGA,SZGGA2


      LOGICAL   LDUM,LFOUND,LDIAG,LREAL,LREALD,LPDENS,LTET,LOPTICS, &
     &          LCORR,LOPEN,lmusic,LWAVE,LCHARG,LVTOT,LCOR_DIP, &
     &          LORBIT_,LELF,LCOMPAT,LPARD,LPAW,LADDGRID,LCRITICAL_MEM, &
     &          LNONCOLLINEAR,LSORBIT,LMETAGGA,LPLANE_WISE,LASPH,LSECVAR
      INTEGER   IGGA2
      DIMENSION POMASS(NTYPD),RWIGS(NTYPP), &
     &          ROPT(NTYPD)
      DIMENSION ATOMOM(*)
      DIMENSION POSCEN(3)
      REAL(q)   SAXIS(3)
      REAL(q)   NELECT,NUP_DOWN
      REAL(q)   STM(5)
!-MM- Spin spiral stuff      
      LOGICAL   LSPIRAL,LZEROZ
      REAL(q)   QSPIRAL(3)
!-MM- end of addition

! 'title'-string (defaults to 'unknown system'), keyword 'SYSTEM'
      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')

      SZNAM='unknown system'
      CALL RDATAB(LOPEN,'INCAR',IU5,'SYSTEM','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,SZNAM,N,40,IERR)
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SYSTEM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('SYSTEM','S',IDUM,RDUM,CDUM,LDUM,SZNAM,N)

      CALL STRIP(SZNAM,N,'L')
      SZNAM1=SZNAM
! start flag ISTART: a default value ISTART=1 should do the best job!
      ISTART=1
! ... of course if 'WAVECAR' doesnt exist --> take ISTART=0 ...
      IF (.NOT.LFOUND) ISTART=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISTART','=','#',';','I', &
     &            ISTART,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISTART'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ISTART','I',ISTART,RDUM,CDUM,LDUM,CHARAC,N)
! the 'idiot flag' (for VTUTOR ...), defaults to 3 ('complete idiot')
      IDIOT=3
      CALL RDATAB(LOPEN,'INCAR',IU5,'IDIOT','=','#',';','I', &
     &            IDIOT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IDIOT'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (IDIOT<0) IDIOT=0
      IF (IDIOT>3) IDIOT=3
      CALL XML_INCAR('IDIOT','I',IDIOT,RDUM,CDUM,LDUM,CHARAC,N)
! ... read in the required precission (low - medium - high)
      SZNAM='MEDIUM'
      CALL RDATAB(LOPEN,'INCAR',IU5,'PREC','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,SZNAM,N,40,IERR)
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''PREC'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL STRIP(SZNAM,N,'L')
      CALL LOWER(SZNAM)
      SZPREC=SZNAM
      CALL XML_INCAR('PREC','S',IDUM,RDUM,CDUM,LDUM,SZNAM,N)
! ... read in the GGA tag (enforces GGA, even if GGA was not used for
!     PP creation)
      SZNAM='--'
      CALL RDATAB(LOPEN,'INCAR',IU5,'GGA','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,SZNAM,N,40,IERR)
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''GGA'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('GGA','S',IDUM,RDUM,CDUM,LDUM,SZNAM,N)

      CALL STRIP(SZNAM,N,'L')
      CALL UPPER(SZNAM)
      SZGGA=SZNAM
! algorithm: default is 8 (prec. CG)
      IALGO=38
      CALL RDATAB(LOPEN,'INCAR',IU5,'IALGO','=','#',';','I', &
     &            IALGO,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IALGO'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('IALGO','I',IALGO,RDUM,CDUM,LDUM,CHARAC,N)
! algorithm: tag ALGO overwrites IALGO 
      CALL RDATAB(LOPEN,'INCAR',IU5,'ALGO','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,INPLIN,N,40,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ALGO'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ALGO','S',IDUM,RDUM,CDUM,LDUM,INPLIN,N)
      CALL STRIP(INPLIN,N,'L')
      IF ( INPLIN(1:1)=='F' .OR. INPLIN(1:1)=='f') THEN
         IALGO=68
      ELSE IF ( INPLIN(1:1)=='V' .OR. INPLIN(1:1)=='v') THEN
         IALGO=48
      ELSE IF ( INPLIN(1:1)=='N' .OR. INPLIN(1:1)=='N') THEN
         IALGO=38
      ENDIF
! LDIAG -- use subspace diagonalization or not (default is TRUE):
      LDIAG=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LDIAG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LDIAG,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LDIAG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LDIAG','L',IDUM,RDUM,CDUM,LDIAG,CHARAC,N)
! LADDGRID -- use an additional grid for the calculation of the US-PP
      LADDGRID=.FALSE.
!      IF (SZPREC(1:1)=='a') LADDGRID=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'ADDGRID','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LADDGRID,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LADDGRID'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ADDGRID','L',IDUM,RDUM,CDUM,LADDGRID,CHARAC,N)
! read in flag LSORBIT
      LSORBIT=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LSORBIT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LSORBIT,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LSORBIT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LSORBIT','L',IDUM,RDUM,CDUM,LSORBIT,CHARAC,N)
! read in flag LNONCOLLINEAR
      LNONCOLLINEAR=LSORBIT
      CALL RDATAB(LOPEN,'INCAR',IU5,'LNONCOLLINEAR','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LNONCOLLINEAR,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LNONCOLLINEAR'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (LSORBIT) LNONCOLLINEAR=LSORBIT
      CALL XML_INCAR('LNONCOLLINEAR','L',IDUM,RDUM,CDUM,LNONCOLLINEAR,CHARAC,N)

! ... read spin quantisation axis
      SAXIS(1)=0
      SAXIS(2)=0
      SAXIS(3)=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'SAXIS','=','#',';','F', &
     &            IDUM,SAXIS,CDUM,LDUM,CHARAC,N,3,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N/=3))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SAXIS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('SAXIS','F',IDUM,SAXIS,CDUM,LDUM,CHARAC,N)

#if defined(NGXhalf) ||  defined(NGZhalf)
      IF (LNONCOLLINEAR) THEN
         WRITE(*,*) 'ERROR: non collinear calculations require that VASP is compiled'
         WRITE(*,*) ' without the flag -DNGXhalf and -DNGZhalf'
         STOP
      ENDIF
#endif
! spin polarized calculation? (1 is default)
      ISPIN=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISPIN','=','#',';','I', &
     &            ISPIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISPIN'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (ISPIN>=2) ISPIN=2
      IF (ISPIN<=1) ISPIN=1

      CALL XML_INCAR('ISPIN','I',ISPIN,RDUM,CDUM,LDUM,CHARAC,N)

! Mixing parameters: by default use IMIX=4 (Broyden) with AMIX=0.8,
! BMIX=1.0 (should work almost always ...), WC=100, INIMIX=1, MIXPRE=1
      IMIX=4
      CALL RDATAB(LOPEN,'INCAR',IU5,'IMIX','=','#',';','I', &
     &            IMIX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('IMIX','I',IMIX,RDUM,CDUM,LDUM,CHARAC,N)

      MAXMIX=-45
      CALL RDATAB(LOPEN,'INCAR',IU5,'MAXMIX','=','#',';','I', &
     &            MAXMIX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''MAXMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('MAXMIX','I',MAXMIX,RDUM,CDUM,LDUM,CHARAC,N)

      MREMOVE=5
      CALL RDATAB(LOPEN,'INCAR',IU5,'MREMOVE','=','#',';','I', &
     &            MREMOVE,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''MREMOVE'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('MREMOVE','I',MREMOVE,RDUM,CDUM,LDUM,CHARAC,N)


      AMIX=0.8_q; IF (ISPIN == 2) AMIX = 0.4_q 
      IF (LPAW) AMIX=0.4_q

      CALL RDATAB(LOPEN,'INCAR',IU5,'AMIX','=','#',';','F', &
     &            IDUM,AMIX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''AMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('AMIX','F',IDUM,AMIX,CDUM,LDUM,CHARAC,N)

      BMIX=1.0_q;
      IF (LPAW) BMIX=1.0_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'BMIX','=','#',';','F', &
     &            IDUM,BMIX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''BMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('BMIX','F',IDUM,BMIX,CDUM,LDUM,CHARAC,N)

      AMIN=0.1_q;
      CALL RDATAB(LOPEN,'INCAR',IU5,'AMIN','=','#',';','F', &
     &            IDUM,AMIN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''AMIN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('AMIN','F',IDUM,AMIN,CDUM,LDUM,CHARAC,N)

      AMIX_MAG=AMIX*4
      CALL RDATAB(LOPEN,'INCAR',IU5,'AMIX_MAG','=','#',';','F', &
     &            IDUM,AMIX_MAG,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''AMIX_MAG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('AMIX_MAG','F',IDUM,AMIX_MAG,CDUM,LDUM,CHARAC,N)

      BMIX_MAG=BMIX
      CALL RDATAB(LOPEN,'INCAR',IU5,'BMIX_MAG','=','#',';','F', &
     &            IDUM,BMIX_MAG,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''BMIX_MAG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('BMIX_MAG','F',IDUM,BMIX_MAG,CDUM,LDUM,CHARAC,N)

      WC=100._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'WC','=','#',';','F', &
     &            IDUM,WC,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''WC'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('WC','F',IDUM,WC,CDUM,LDUM,CHARAC,N)

      INIMIX=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'INIMIX','=','#',';','I', &
     &            INIMIX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''INIMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('INIMIX','I',INIMIX,RDUM,CDUM,LDUM,CHARAC,N)

      MIXPRE=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'MIXPRE','=','#',';','I', &
     &            MIXPRE,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''MIXPRE'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('MIXPRE','I',MIXPRE,RDUM,CDUM,LDUM,CHARAC,N)

! initial charge density ICHARG (default 0, if startjob: default 2)
      ICHARG=0
      IF (ISTART==0) ICHARG=2
      CALL RDATAB(LOPEN,'INCAR',IU5,'ICHARG','=','#',';','I', &
     &            ICHARG,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ICHARG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ICHARG','I',ICHARG,RDUM,CDUM,LDUM,CHARAC,N)

      LPDENS=.FALSE.
      IF (ICHARG<0) THEN
         ICHARG=0
         LPDENS=.TRUE.
      ENDIF
! initial wavefunctions (defaults is 1, warning: keyword is 'INIWAV')
      INIWAV=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'INIWAV','=','#',';','I', &
     &            INIWAV,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''INIWAV'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('INIWAV','I',INIWAV,RDUM,CDUM,LDUM,CHARAC,N)

! max/min. number of electronic minimization steps, delay ... (default
! shall be NELM=60, NELMIN=2, NELMDL=-5 if ISTART=0 and
! NELM=60, NELMIN=2, NELMDL=0 if ISTART/=0 ...):
      NELM=60
      CALL RDATAB(LOPEN,'INCAR',IU5,'NELM','=','#',';','I', &
     &            NELM,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NELM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NELM','I',NELM,RDUM,CDUM,LDUM,CHARAC,N)

      NELMDL=0
      IF (ISTART==0 .AND. INIWAV==1) THEN
         NELMDL=-5
         IF (IALGO>=40 .AND. IALGO<=50) NELMDL=-12
      ENDIF
      CALL RDATAB(LOPEN,'INCAR',IU5,'NELMDL','=','#',';','I', &
     &            NELMDL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NELMDL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NELMDL','I',NELMDL,RDUM,CDUM,LDUM,CHARAC,N)

      NELMIN=2
      CALL RDATAB(LOPEN,'INCAR',IU5,'NELMIN','=','#',';','I', &
     &            NELMIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NELMIN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NELMIN','I',NELMIN,RDUM,CDUM,LDUM,CHARAC,N)

! energy tolerances (defaults: EDIFF=1E-4, EDIFFG=1E-3)
      EDIFF=1.E-4_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'EDIFF','=','#',';','F', &
     &            IDUM,EDIFF,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EDIFF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('EDIFF','F',IDUM,EDIFF,CDUM,LDUM,CHARAC,N)

! for reasons of safety (crazy user are present all over the world):
      EDIFF=MAX(ABS(EDIFF),1.E-9_q)
      EDIFFG=EDIFF*10
      CALL RDATAB(LOPEN,'INCAR',IU5,'EDIFFG','=','#',';','F', &
     &            IDUM,EDIFFG,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EDIFFG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('EDIFFG','F',IDUM,EDIFFG,CDUM,LDUM,CHARAC,N)

! conjugate gradient or quasi-Newton method? (default IBRION=0)
      IBRION=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'IBRION','=','#',';','I', &
     &            IBRION,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IBRION'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('IBRION','I',IBRION,RDUM,CDUM,LDUM,CHARAC,N)

! number of degrees of freedom
      IF (IBRION==2) THEN
         NFREE=1
      ELSE
         NFREE=0
      ENDIF
      CALL RDATAB(LOPEN,'INCAR',IU5,'NFREE','=','#',';','I', &
     &            NFREE,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NFREE'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NFREE','I',NFREE,RDUM,CDUM,LDUM,CHARAC,N)

! number of ionic steps, calculate stresses? (default NSW=0, ISIF=2):
      NSW=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'NSW','=','#',';','I', &
     &            NSW,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NSW'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NSW','I',NSW,RDUM,CDUM,LDUM,CHARAC,N)

! IBRION is 'useless' if NSW=0, set this flag to -1 in this case ...
      IF (NSW==0) IBRION=-1
      ISIF=2
! if MD is selected dont calculate stress
      IF (IBRION==0) ISIF=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISIF','=','#',';','I', &
     &            ISIF,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISIF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ISIF','I',ISIF,RDUM,CDUM,LDUM,CHARAC,N)
! prediction of wavefunction:
      IWAVPR=0
! MDs
      IF (IBRION==0) IWAVPR=2
#ifdef MPI
      IF (IWAVPR > 0) THEN
         IWAVPR=IWAVPR+10
      ENDIF
#endif
! relaxation: IWAVPR=1
      IF (IBRION>0) IWAVPR=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'IWAVPR','=','#',';','I', &
     &            IWAVPR,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IWAVPR'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (IWAVPR==10) THEN
        ! MD: IWAVPR=12
        IF (IBRION==0) IWAVPR=12
        ! relaxation: IWAVPR=11
        IF (IBRION>0) IWAVPR=11
      ENDIF
      IF (IWAVPR==1) IWAVPR=11 ! makes the same but requires less memory :->

#ifdef MPI
      IWAVPR=MOD(IWAVPR,10)+10
#endif
      CALL XML_INCAR('IWAVPR','I',IWAVPR,RDUM,CDUM,LDUM,CHARAC,N)

! switch on symmetry (default ISYM=1):
      ISYM=1 ; IF (LPAW) ISYM=2
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISYM','=','#',';','I', &
     &            ISYM,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISYM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ISYM','I',ISYM,RDUM,CDUM,LDUM,CHARAC,N)

! for reasons of safety (crazy user are present all over the world):
      TINY=1E-5
      CALL RDATAB(LOPEN,'INCAR',IU5,'SYMPREC','=','#',';','F', &
     &            IDUM,TINY,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SYMPREC'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('SYMPREC','F',IDUM,TINY,CDUM,LDUM,CHARAC,N)

! how often to write some data; defaults to KBLOCK=NSW, NBLOCK=1:
      NBLOCK=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NBLOCK','=','#',';','I', &
     &            NBLOCK,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NBLOCK'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NBLOCK','I',NBLOCK,RDUM,CDUM,LDUM,CHARAC,N)

      KBLOCK=NSW
      CALL RDATAB(LOPEN,'INCAR',IU5,'KBLOCK','=','#',';','I', &
     &            KBLOCK,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''KBLOCK'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('KBLOCK','I',KBLOCK,RDUM,CDUM,LDUM,CHARAC,N)

! plane wave cutoff energy for wavefunctions ..., no default!!!!
      ENMAX=-1._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'ENMAX','=','#',';','F', &
     &            IDUM,ENMAX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ENMAX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ENMAX','F',IDUM,ENMAX,CDUM,LDUM,CHARAC,N)

      CALL RDATAB(LOPEN,'INCAR',IU5,'ENCUT','=','#',';','F', &
     &            IDUM,ENMAX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ENCUT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ENCUT','F',IDUM,ENMAX,CDUM,LDUM,CHARAC,N)

! plane wave cutoff energy for wavefunctions ..., no default!!!!
      ENINI=-1._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'ENINI','=','#',';','F', &
     &            IDUM,ENINI,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ENINI'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ENINI','F',IDUM,ENINI,CDUM,LDUM,CHARAC,N)

! cutoff for augmentation charge
      ENAUG=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'ENAUG','=','#',';','F', &
     &            IDUM,ENAUG,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ENAUG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ENAUG','F',IDUM,ENAUG,CDUM,LDUM,CHARAC,N)

! read in NGX, NGY, NGZ, NBANDS
      NGX=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGX','=','#',';','I', &
     &            NGX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGX','I',NGX,RDUM,CDUM,LDUM,CHARAC,N)

      NGY=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGY','=','#',';','I', &
     &            NGY,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGY'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGY','I',NGY,RDUM,CDUM,LDUM,CHARAC,N)

      NGZ=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGZ','=','#',';','I', &
     &            NGZ,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGZ'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGZ','I',NGZ,RDUM,CDUM,LDUM,CHARAC,N)

      NGXF=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGXF','=','#',';','I', &
     &            NGXF,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGXF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGXF','I',NGXF,RDUM,CDUM,LDUM,CHARAC,N)

      NGYF=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGYF','=','#',';','I', &
     &            NGYF,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGYF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGYF','I',NGYF,RDUM,CDUM,LDUM,CHARAC,N)

      NGZF=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NGZF','=','#',';','I', &
     &            NGZF,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NGZF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NGZF','I',NGZF,RDUM,CDUM,LDUM,CHARAC,N)

      NBANDS=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NBANDS','=','#',';','I', &
     &            NBANDS,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NBANDS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NBANDS','I',NBANDS,RDUM,CDUM,LDUM,CHARAC,N)

! ionic time step, default is POTIM=0.5. for IBRION/=0, else no default!
      POTIM=0.5_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'POTIM','=','#',';','F', &
     &            IDUM,POTIM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''POTIM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('POTIM','F',IDUM,POTIM,CDUM,LDUM,CHARAC,N)

! if IBRION=0 (MD) then POTIM must be given, otherwise error ... !
      IF (((IERR==3).OR.(N<1)).AND.(IBRION==0)) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Fatal error! IBRION=0, but no entry for POTIM'// &
     &               ' on file INCAR. MUST be specified!!'
         IF (IU0>=0) &
         WRITE(IU0,*)'                                          '// &
     &               '                ----'
         STOP
      ENDIF
! start temperature and end temperature (default is 1E-4 for both),
      TEBEG=1.E-4_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'TEBEG','=','#',';','F', &
     &            IDUM,TEBEG,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''TEBEG'' from file INCAR.'
         GOTO 150
      ENDIF
      TEEND=TEBEG
      CALL XML_INCAR('TEBEG','F',IDUM,TEBEG,CDUM,LDUM,CHARAC,N)

      CALL RDATAB(LOPEN,'INCAR',IU5,'TEEND','=','#',';','F', &
     &            IDUM,TEEND,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''TEEND'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('TEEND','F',IDUM,TEEND,CDUM,LDUM,CHARAC,N)

! pair-correlation functions ..., defaults are NPACO=256, APACO=10
      NPACO=256
      CALL RDATAB(LOPEN,'INCAR',IU5,'NPACO','=','#',';','I', &
     &            NPACO,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NPACO'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NPACO','I',NPACO,RDUM,CDUM,LDUM,CHARAC,N)

      APACO=16._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'APACO','=','#',';','F', &
     &            IDUM,APACO,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''APACO'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('APACO','F',IDUM,APACO,CDUM,LDUM,CHARAC,N)

! NEDOS subdivisions for DOS
      NEDOS=301
      CALL RDATAB(LOPEN,'INCAR',IU5,'NEDOS','=','#',';','I', &
     &            NEDOS,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NEDOS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NEDOS','I',NEDOS,RDUM,CDUM,LDUM,CHARAC,N)

! NBLK blocking for some DGEMM commands
      NBLK=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'NBLK','=','#',';','I', &
     &            NBLK,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NBLK'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NBLK','I',NBLK,RDUM,CDUM,LDUM,CHARAC,N)

! default for SMASS is -3 (micro canonical MD)
      SMASS=-3
      CALL RDATAB(LOPEN,'INCAR',IU5,'SMASS','=','#',';','F', &
     &            IDUM,SMASS,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SMASS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('SMASS','F',IDUM,SMASS,CDUM,LDUM,CHARAC,N)

! Well, we supply the atomic masses on file POTCAR, but in some cases
! one might wish to change them artificially (--> for example trying
! some kind of 'pre-conditioning' by hand for relaxation runs ...):
! by default we set all masses to negative numbers (this shall be the
! 'signal' to take the values from file POTCAR ...).
      POMASS=-1._q
      RWIGS=-1._q
      ROPT=0

      CALL RDATAB(LOPEN,'INCAR',IU5,'POMASS','=','#',';','F', &
     &            IDUM,POMASS,CDUM,LDUM,CHARAC,N,NTYPIN,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<NTYPIN))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''POMASS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('POMASS','F',IDUM,POMASS,CDUM,LDUM,CHARAC,N)

! "Cutoff radii" (Wigner-Seitz-radii) for l-projections (default is -1.)
      CALL RDATAB(LOPEN,'INCAR',IU5,'RWIGS','=','#',';','F', &
     &            IDUM,RWIGS,CDUM,LDUM,CHARAC,N,NTYPP,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<NTYPP))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''RWIGS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('RWIGS','F',IDUM,RWIGS,CDUM,LDUM,CHARAC,N)

! number of up down electrons
      NUP_DOWN=-1._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'NUPDOWN','=','#',';','F', &
     &            IDUM,NUP_DOWN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NUPDOWN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NUPDOWN','F',IDUM,NUP_DOWN,CDUM,LDUM,CHARAC,N)

! Initial magnetic moments for each atom (default is 1. for all ions)
      AINI=1
      IF (NUP_DOWN >=0) THEN
        AINI=NUP_DOWN/ NIONS
      ENDIF

      NMAGMOM=NIONS
      IF (LNONCOLLINEAR) NMAGMOM=3*NIONS
      DO NI=1,NMAGMOM
         ATOMOM(NI)=AINI
      ENDDO
      CALL RDATAB(LOPEN,'INCAR',IU5,'MAGMOM','=','#',';','F', &
     &            IDUM,ATOMOM,CDUM,LDUM,CHARAC,N,NMAGMOM,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<NMAGMOM))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''MAGMOM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('MAGMOM','F',IDUM,ATOMOM,CDUM,LDUM,CHARAC,N)

! number of electrons ..., default is NELECT=0 (= neutral cell)
      NELECT=0._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'NELECT','=','#',';','F', &
     &            IDUM,NELECT,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NELECT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NELECT','F',IDUM,NELECT,CDUM,LDUM,CHARAC,N)

! Real-space projection: default should be POTCAR-dependent ... (if
! one finds 'optimization flag' then set LREAL=.TRUE., else .FALSE.)
      LREAL=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LREAL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LREAL,CHARAC,N,1,IERR)
      CALL XML_INCAR('LREAL','L',IDUM,RDUM,CDUM,LREAL,CHARAC,N)

! no input --> remind it and choose later the appropriate value ...
      LREALD=(IERR==3)
      IF (IERR==5) THEN
        CALL RDATAB(LOPEN,'INCAR',IU5,'LREAL','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,INPLIN,N,40,IERR)
        LREAL=.TRUE.
        CALL XML_INCAR('LREAL','S',IDUM,RDUM,CDUM,LDUM,INPLIN,N)

        CALL STRIP(INPLIN,N,'L')
        IF (INPLIN(1:1)=='O' .OR. INPLIN(1:1)=='o' .OR. &
            INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
          IF ( INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
            ROPTV=-2E-3
            IF  (SZPREC(1:1)=='l') ROPTV=-1E-2
            IF  (SZPREC(1:1)=='n') ROPTV=-5E-4
            IF  (SZPREC(1:1)=='h') ROPTV=-4E-4
            IF  (SZPREC(1:1)=='a') ROPTV=-2.5E-4
          ELSE
             ROPTV=1.0_q
             IF  (SZPREC(1:1)=='l') ROPTV=1/1.5
             IF  (SZPREC(1:1)=='h') ROPTV=1.5
          ENDIF
          DO NTYP=1,NTYPIN
            ROPT(NTYP)=ROPTV
          ENDDO
          CALL RDATAB(LOPEN,'INCAR',IU5,'ROPT','=','#',';','F', &
     &            IDUM,ROPT,CDUM,LDUM,CHARAC,N,NTYPIN,IERR)
          IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<NTYPIN))) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''ROPT'' from file INCAR.'
            GOTO 150
          ENDIF
          IF ( INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
            ROPT=-ABS(ROPT)
          ELSE
            ROPT=ABS(ROPT)
          ENDIF
          CALL XML_INCAR_V('ROPT','F',IDUM,ROPT,CDUM,LDUM,CHARAC,N)
        ELSE
          IERR=5
        ENDIF
      ENDIF
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LREAL'' from file INCAR.'
         GOTO 150
      ENDIF
! plane by plane distribution of data
!     LPLANE_WISE=.NOT. LREAL
      LPLANE_WISE=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LPLANE','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LPLANE_WISE,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LPLANE_WISE'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LPLANE','L',IDUM,RDUM,CDUM,LPLANE_WISE,CHARAC,N)

! LCOMPAT .TRUE. means full compatibility
      LCOMPAT = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LCOMPAT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LCOMPAT,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LCOMPAT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LCOMPAT','L',IDUM,RDUM,CDUM,LCOMPAT,CHARAC,N)
! electronic timestep
      TIME=0.1_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'TIME','=','#',';','F', &
     &            IDUM,TIME,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''TIME'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('TIME','F',IDUM,TIME,CDUM,LDUM,CHARAC,N)

! energy range for DOS (default is EMIN=10.,EMAX=-10. = automatic):
      EMIN=10._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'EMIN','=','#',';','F', &
     &            IDUM,EMIN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EMIN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('EMIN','F',IDUM,EMIN,CDUM,LDUM,CHARAC,N)

      EMAX=-10._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'EMAX','=','#',';','F', &
     &            IDUM,EMAX,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EMAX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('EMAX','F',IDUM,EMAX,CDUM,LDUM,CHARAC,N)

! z range for STM data
      STM=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'STM','=','#',';','F', &
     &            IDUM,STM,CDUM,LDUM,CHARAC,N,5,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<5))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''STM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('STM','F',IDUM,STM,CDUM,LDUM,CHARAC,N)

! BZ-integration type, default is ISMEAR=1 and SIGMA=0.2 ...
      ISMEAR=1
      CALL RDATAB(LOPEN,'INCAR',IU5,'ISMEAR','=','#',';','I', &
     &            ISMEAR,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''ISMEAR'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('ISMEAR','I',ISMEAR,RDUM,CDUM,LDUM,CHARAC,N)

      SIGMA=0.2_q
! If we provide fermi-weights on file INCAR the main intention is mostly
! to do calculations at given fixed occupancies -> this requires SIGMA=0
      IF (ISMEAR==-2) SIGMA=0._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'SIGMA','=','#',';','F', &
     &            IDUM,SIGMA,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''SIGMA'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('SIGMA','F',IDUM,SIGMA,CDUM,LDUM,CHARAC,N)

! ISMEAR<=-4 and ISMEAR>=30 means tetrahedron method for DOS ...,
! ISMEAR==-4,-5 and <=-7: also tetrahedron method for occ. numbers
      LTET=((ISMEAR<=-4).OR.(ISMEAR>=30))
      IF (ISMEAR==-6) ISMEAR=-1
      IF (ISMEAR>=0) ISMEAR=MOD(ISMEAR,30)
! min. occupation number for 'high quality update' (default: WEIMIN=0)
      WEIMIN=0._q
! MD and relaxation: take WEIMIN=0.001
      IF (IBRION>=0) WEIMIN=0.001_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'WEIMIN','=','#',';','F', &
     &            IDUM,WEIMIN,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''WEIMIN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('WEIMIN','F',IDUM,WEIMIN,CDUM,LDUM,CHARAC,N)

! break condition for intra-band min. (default: 0.25*EDIFF/NBANDS)
! because we allow to be EDIFFG smaller than EDIFF also consider
! EDIFFG
      EBREAK=-1
      CALL RDATAB(LOPEN,'INCAR',IU5,'EBREAK','=','#',';','F', &
     &            IDUM,EBREAK,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''EBREAK'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('EBREAK','F',IDUM,EBREAK,CDUM,LDUM,CHARAC,N)

! relative break condition for intra-band minimization (default is 0.3)
      DEPER=0.3_q
      CALL RDATAB(LOPEN,'INCAR',IU5,'DEPER','=','#',';','F', &
     &            IDUM,DEPER,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''DEPER'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('DEPER','F',IDUM,DEPER,CDUM,LDUM,CHARAC,N)

! 'verbosity' (default: 2):
      NWRITE=2
      CALL RDATAB(LOPEN,'INCAR',IU5,'NWRITE','=','#',';','I', &
     &            NWRITE,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NWRITE'' from file INCAR.'
         GOTO 150
      ENDIF
   ! allowed range is 0...4, if <0 assume 0, if >4 assume 4 ...
      IF (NWRITE<0) NWRITE=0
      IF (NWRITE>4) NWRITE=4
      CALL XML_INCAR('NWRITE','I',NWRITE,RDUM,CDUM,LDUM,CHARAC,N)

! Harris corrections for Hellman-Feynman forces ... (default yes):
      LCORR=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LCORR','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LCORR,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LCORR'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LCORR','L',IDUM,RDUM,CDUM,LCORR,CHARAC,N) 
! Pullay pressure (one could also say external pressure), default 0.
      PSTRESS=0._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'PSTRESS','=','#',';','F', &
     &            IDUM,PSTRESS,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''PSTRESS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('PSTRESS','F',IDUM,PSTRESS,CDUM,LDUM,CHARAC,N)

! max. number of iterations NRMM in RMM-DIIS (NDAV), default usually 4
      NDAV=4
      CALL RDATAB(LOPEN,'INCAR',IU5,'NRMM','=','#',';','I', &
     &            NDAV,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NRMM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NRMM','I',NDAV,RDUM,CDUM,LDUM,CHARAC,N)

! band blocking in RMM-DIIS and Davidson (and some other subroutines)
      NSIM=4
      ! for the Davidson it is advisable to increase blocking
      CALL RDATAB(LOPEN,'INCAR',IU5,'NSIM','=','#',';','I', &
     &            NSIM,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''NSIM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('NSIM','I',NSIM,RDUM,CDUM,LDUM,CHARAC,N)

! max. L for onsite charge expansion in PAW method
      LMAX_CALC=-100
      CALL RDATAB(LOPEN,'INCAR',IU5,'LMAXPAW','=','#',';','I', &
     &            LMAX_CALC,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LMAXPAW'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LMAXPAW','I',LMAX_CALC,RDUM,CDUM,LDUM,CHARAC,N)
! max. L for the mixing and CHGCAR for the onsite charge expansion in PAW method
      LMAX_MIX=2
      CALL RDATAB(LOPEN,'INCAR',IU5,'LMAXMIX','=','#',';','I', &
     &            LMAX_MIX,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LMAXMIX'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LMAXMIX','I',LMAX_MIX,RDUM,CDUM,LDUM,CHARAC,N)

! some "music" ? (--> default is no ...)
      lmusic=.false.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LMUSIC','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LMUSIC,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LMUSIC'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LMUSIC','L',IDUM,RDUM,CDUM,LMUSIC,CHARAC,N)

! Sometimes we not interested in any WAVECAR file at all ...
      LWAVE=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LWAVE','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LWAVE,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LWAVE'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LWAVE','L',IDUM,RDUM,CDUM,LWAVE,CHARAC,N)

! ... and maybe not even in any CHGCAR / CHG file ...
      LCHARG=.TRUE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LCHARG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LCHARG,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LCHARG'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LCHARG','L',IDUM,RDUM,CDUM,LCHARG,CHARAC,N)

! ... interested in partial charge density ?
      LPARD = .FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LPARD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LPARD,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LPARD'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LPARD','L',IDUM,RDUM,CDUM,LPARD,CHARAC,N)

! ... a WAVECAR must exist
      IF (.NOT.LFOUND) LPARD = .FALSE.
! ... but maybe in the total potential?
      LVTOT=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LVTOT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LVTOT,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LVTOT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LVTOT','L',IDUM,RDUM,CDUM,LVTOT,CHARAC,N)

! ... Dipol required
!
      IDIPCO=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'IDIPOL','=','#',';','I', &
     &            IDIPCO,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N/=1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''IDIPOL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('IDIPOL','I',IDIPCO,RDUM,CDUM,LDUM,CHARAC,N)

      LCOR_DIP=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LDIPOL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LCOR_DIP,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LDIPOL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LDIPOL','L',IDUM,RDUM,CDUM,LCOR_DIP,CHARAC,N)

! ... Dipol required
!
      POSCEN(1)=-100
      POSCEN(2)=-100
      POSCEN(3)=-100
      CALL RDATAB(LOPEN,'INCAR',IU5,'DIPOL','=','#',';','F', &
     &            IDUM,POSCEN,CDUM,LDUM,CHARAC,N,3,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N/=3))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''DIPOL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR_V('DIPOL','F',IDUM,POSCEN,CDUM,LDUM,CHARAC,N)

! read in flag LORBIT
      LORBIT_=.FALSE.
      LORBIT=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'LORBIT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LORBIT_,CHARAC,N,1,IERR)
      IF (IERR==5) THEN
      CALL RDATAB(LOPEN,'INCAR',IU5,'LORBIT','=','#',';','I', &
     &            LORBIT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      ELSE
         IF (LORBIT_) LORBIT=5
      ENDIF

      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LORBIT'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LORBIT','L',IDUM,RDUM,CDUM,LORBIT_,CHARAC,N)

! spin interpolation according to Vosko Wil and Nusair
      LFCI=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'VOSKOWN','=','#',';','I', &
     &            LFCI,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''VOSKOWN'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('VOSKOWN','I',LFCI,RDUM,CDUM,LDUM,CHARAC,N)

! read in flag LELF
      LELF=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LELF','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LELF,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LELF'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LELF','L',IDUM,RDUM,CDUM,LELF,CHARAC,N)

! read in flag LCRITICAL_MEM
      LCRITICAL_MEM=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LCRITICAL_MEM','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LCRITICAL_MEM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LCRITICAL_MEM'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LCRITICAL_MEM','L',IDUM,RDUM,CDUM,LCRITICAL_MEM,CHARAC,N)

! read in flag LOPTICS
      LOPTICS=.FALSE.
#ifndef noPAW
      CALL RDATAB(LOPEN,'INCAR',IU5,'LOPTICS','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LOPTICS,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LOPTICS'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LOPTICS','L',IDUM,RDUM,CDUM,LOPTICS,CHARAC,N)
#endif
! if scaLAPACK is used it can be switched of in the INCAR file
      IF (LscaLAPACK) THEN
      CALL RDATAB(LOPEN,'INCAR',IU5,'LSCALAPACK','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LscaLAPACK,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LSCALPACK'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (IU0>0) THEN
         IF (LscaLAPACK) THEN
            WRITE(IU0,*) 'scaLAPACK will be used'
         ELSE
            WRITE(IU0,*) 'scaLAPACK is switched off'
         ENDIF
      ENDIF
      ENDIF
      CALL XML_INCAR('LSCALAPACK','L',IDUM,RDUM,CDUM,LscaLAPACK,CHARAC,N)

      LSCALU= LSCALAPACK
! the parallel LU decomposition might be slower than the serial 
! one, hence we can switch it off
      CALL RDATAB(LOPEN,'INCAR',IU5,'LSCALU','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LSCALU,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LSCALU'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LSCALU','L',IDUM,RDUM,CDUM,LSCALU,CHARAC,N)

! try of overlap communication with calculations ?
      CALL RDATAB(LOPEN,'INCAR',IU5,'LASYNC','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LASYNC,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LASYNC'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LASYNC','L',IDUM,RDUM,CDUM,LASYNC,CHARAC,N)

! read in flag LASPH
      LASPH=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LASPH','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LASPH,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LASPH'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LASPH','L',IDUM,RDUM,CDUM,LASPH,CHARAC,N)

! read in flag LMETAGGA 
      LMETAGGA=.FALSE.
      CALL RDATAB(LOPEN,'INCAR',IU5,'LMETAGGA','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LMETAGGA,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LMETAGGA'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LMETAGGA','L',IDUM,RDUM,CDUM,LMETAGGA,CHARAC,N)

! set LASPH if LMETAGGA is chosen (metagga only calculated aspherically)
      IF (LMETAGGA) LASPH=.TRUE.
! ... read in the GGA tag for second variation
      LSECVAR=.FALSE.
      SZNAM='--'
      CALL RDATAB(LOPEN,'INCAR',IU5,'GGA2','=','#',';','S', &
     &            IDUM,RDUM,CDUM,LDUM,SZNAM,N,40,IERR)
      SZGGA2=SZNAM
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''GGA2'' from file INCAR.'
         GOTO 150
      ENDIF
      IF (IERR/=3) THEN
         CALL STRIP(SZNAM,N,'L')
         CALL UPPER(SZNAM)
         SZGGA2=SZNAM
         CALL EXTYP(SZGGA2,IGGA2)
         IF (IGGA2<0) THEN
            IF (IU0>=0) &
                 WRITE(IU0,*)'Invalid second variational GGA in INCAR.'
            GOTO 152
         ENDIF
         LSECVAR=.TRUE.
      ENDIF
      CALL XML_INCAR('GGA2','S',IDUM,RDUM,CDUM,LDUM,SZNAM,N)

!-MM- spin spiral stuff
! if LSPIRAL
      LSPIRAL=.FALSE.
      LZEROZ =.FALSE.
      QSPIRAL=0._q
      CALL RDATAB(LOPEN,'INCAR',IU5,'LSPIRAL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LSPIRAL,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Error reading item ''LSPIRAL'' from file INCAR.'
         GOTO 150
      ENDIF
      CALL XML_INCAR('LSPIRAL','L',IDUM,RDUM,CDUM,LSPIRAL,CHARAC,N)

! if LSPIRAL=.TRUE. we also need QSPIRAL, and possibly LZEROZ
      IF (LSPIRAL) THEN
! ... read propagation vector of spin spiral
         CALL RDATAB(LOPEN,'INCAR',IU5,'QSPIRAL','=','#',';','F', &
     &               IDUM,QSPIRAL,CDUM,LDUM,CHARAC,N,3,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                       ((IERR==0).AND.(N/=3))) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''QSPIRAL'' from file INCAR.'
            GOTO 150
         ENDIF
         CALL XML_INCAR_V('QSPIRAL','F',IDUM,QSPIRAL,CDUM,LDUM,CHARAC,N)
! ... look for LZEROZ
         LZEROZ=.TRUE.
         CALL RDATAB(LOPEN,'INCAR',IU5,'LZEROZ','=','#',';','L', &
     &               IDUM,RDUM,CDUM,LZEROZ,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                       ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''LZEROZ'' from file INCAR.'
            GOTO 150
         ENDIF
         CALL XML_INCAR('LZEROZ','L',IDUM,RDUM,CDUM,LZEROZ,CHARAC,N)
      ENDIF           
! Thats all from INCAR (for the first ...):
      CLOSE(IU5)
      RETURN

  150 CONTINUE
      IF (IU0>=0) &
      WRITE(IU0,151) IERR,N
  151 FORMAT(' Error code was IERR=',I1,' ... . Found N=',I5,' data.')
  152 CONTINUE
      STOP

      END
