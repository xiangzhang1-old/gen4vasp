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






      MODULE PSEUDO
      USE prec
      USE radial
      INCLUDE "pseudo.inc"
      CONTAINS

!**************** SUBROUTINE RD_PSEUDO *********************************
! RCS:  $Id: pseudo.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
!  reads in all pseudopotential from the POTCAR file
!
!  check if LDIM and LMDIM is sufficient
!  if not tell the user to increase those numbers
!
!***********************************************************************

      SUBROUTINE RD_PSEUDO(INFO,P, &
     &           NTYP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           POMASS,RWIGS, &
     &           IU0,IU6,NWRITE,LPAW)
      USE prec
      USE base
      USE ini
      IMPLICIT NONE

      TYPE(INFO_STRUCT) INFO
      INTEGER  NTYP,NTYPD      ! number of types (acutal/dimension)
      INTEGER LDIM,LDIM2,LMDIM ! see pseudo.inc
      TYPE (potcar) P(NTYPD)   ! PP information
      LOGICAL LPAW
      INTEGER IU6,IU0          ! where I/O goes
! wigner seitz radius and mass found on INCAR
      REAL(q)    RWIGS(NTYPD),POMASS(NTYPD)
! dynamical work space
      INTEGER,PARAMETER :: ISDIM=100
      CHARACTER*80 STRING(ISDIM)
      CHARACTER*40 FORM
! temporary varibales
      INTEGER IDUMM1,IDUMM2,IDUMM3,L2,NREAD,NL,IC,L,LP,I,J,K,CHANNELS,LMAX
      LOGICAL LDUM
      CHARACTER*1  CSEL
      INTEGER IREAD,IWRITE,NWRITE,NMAX

      OPEN(UNIT=10,FILE='POTCAR',STATUS='OLD')

      LPAW = .FALSE.
      INFO%LOVERL=.FALSE.
      INFO%LCORE =.FALSE.
!-----------------------------------------------------------------------
! loop over all pseudpotentials on POTCAR file
!-----------------------------------------------------------------------
      NTYP=1

      forever: DO NTYP=1,NTYPD
      P(NTYP)%LDIM =LDIM
      P(NTYP)%LMDIM=LMDIM
      P(NTYP)%LDIM2=LDIM2
      NULLIFY(P(NTYP)%QPAW)
      NULLIFY(P(NTYP)%QTOT)
      NULLIFY(P(NTYP)%QATO)
      NULLIFY(P(NTYP)%NABLA)
      
      READ(10,'(A40)',END=100,ERR=100) P(NTYP)%SZNAMP

      IF (IU6>=0) &
      WRITE(IU6,*)'POTCAR:  ',P(NTYP)%SZNAMP

      READ(10,*) P(NTYP)%ZVALF
      P(NTYP)%ZVALF_ORIG=P(NTYP)%ZVALF
!-----------------------------------------------------------------------
!  if there is PSCTR header read in and parse information from this
!  header
!-----------------------------------------------------------------------
      READ(10,'(1X,A1)') CSEL
      IF (CSEL /= 'p') THEN
        IF (IU6>=0) &
        WRITE(IU6,*)'this version requires full pseudpotential ', &
     &              ' generation information'
        IF (IU0>=0) &
        WRITE(IU0,*)'this version requires full pseudpotential ', &
     &              ' generation information'
        STOP
      ENDIF
! read contents to STRING
        CALL RDPARA(10,ISDIM,STRING,IREAD,IWRITE)

        IF (NWRITE>=2.AND. IU6>=0 ) &
           WRITE(IU6,'(A80)') (STRING(I),I=1,IWRITE)

! parse from STRING necessary information
        CALL RDPARS(ISDIM,STRING,IREAD,P(NTYP),IU6)
! set and check POMASS
        IF (POMASS(NTYP)<0) POMASS(NTYP)=P(NTYP)%POMASS
        IF (ABS(POMASS(NTYP)-P(NTYP)%POMASS)>1E-3_q) THEN

        IF (NWRITE>=0.AND. IU0>=0) THEN
          WRITE(IU0,*)'WARNING: mass on POTCAR and INCAR are incompatible'
          WRITE(IU0,*)' typ',NTYP,' Mass',POMASS(NTYP),P(NTYP)%POMASS
        ENDIF
        ENDIF

! set RWIGS
        IF (RWIGS(NTYP)==0) RWIGS(NTYP)=P(NTYP)%RWIGS

!-----------------------------------------------------------------------
! local potential, gradient corrections and type of gradient corrections
! partial core and atomic charge density
!-----------------------------------------------------------------------
      READ(10,'(1X,A1)') CSEL
      READ(10,*) P(NTYP)%PSGMAX
      ALLOCATE(P(NTYP)%PSP(NPSPTS,5))

      READ(10,*) (P(NTYP)%PSP (I,2),I=1,NPSPTS)
      DO I=1,NPSPTS
          P(NTYP)%PSP(I,1)=(P(NTYP)%PSGMAX/NPSPTS)*(I-1)
      ENDDO
! PSCORE is equal to V(q) + 4 Pi / q^2
      P(NTYP)%PSCORE=P(NTYP)%PSP(1,2)

      CALL SPLCOF(P(NTYP)%PSP(1,1) ,NPSPTS,NPSPTS,0._q)
      IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' local pseudopotential read in'

!  gradient corrections this is now overwritten by PSCTR header
      READ(10,'(1X,A1)') CSEL
      IF (CSEL=='g') THEN
         READ(10,*)
         READ(10,'(1X,A1)') CSEL
      ENDIF
!
!  partial core and atomic chargedensity
!
      IF (CSEL=='c') THEN
        INFO%LCORE=.TRUE.
        ALLOCATE(P(NTYP)%PSPCOR(NPSPTS))
        IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' partial core-charges read in'

        READ(10,*) (P(NTYP)%PSPCOR (I),I=1,NPSPTS)
        READ(10,*)
      ELSE
         NULLIFY(P(NTYP)%PSPCOR)
      ENDIF

      ALLOCATE(P(NTYP)%PSPRHO(NPSPTS))
      READ(10,*) (P(NTYP)%PSPRHO (I),I=1,NPSPTS)
      IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' atomic valenz-charges read in'

      CHANNELS=0
      P(NTYP)%LMMAX =0
      P(NTYP)%PSMAXN=0
!-----------------------------------------------------------------------
! depletion charges to 0
!-----------------------------------------------------------------------
      ALLOCATE(P(NTYP)%DION(LDIM,LDIM),P(NTYP)%QION(LDIM,LDIM), &
     &         P(NTYP)%LPS(LDIM),P(NTYP)%NLPRO(LDIM), &
     &         P(NTYP)%PSPNL(0:NPSNL,LDIM),P(NTYP)%PSPRNL(NPSRNL,5,LDIM))

! Intel efc compiler workaround (at least version 6.X)
!      P(NTYP)%DION=0
!      P(NTYP)%QION=0
      DO I=1,LDIM
         DO J=1,LDIM
            P(NTYP)%DION(I,J)=0
            P(NTYP)%QION(I,J)=0
         ENDDO
      ENDDO

      P(NTYP)%PSDMAX=0
      NULLIFY(P(NTYP)%QDEP)
      NULLIFY(P(NTYP)%NDEP)
!-----------------------------------------------------------------------
!  read reciprocal projection operators
!-----------------------------------------------------------------------
      READ(10,*,ERR=620,END=600)  P(NTYP)%PSMAXN,LDUM
  620 IF (P(NTYP)%PSMAXN==0) GOTO 600

  610 READ(10,'(1X,A1)',ERR=600,END=600) CSEL
      IF (CSEL== 'D'  .OR. CSEL== 'A' .OR. CSEL=='P' ) GOTO 650
      IF (CSEL== 'E' ) GOTO 600

      READ(10,*)  P(NTYP)%LPS(CHANNELS+1),P(NTYP)%NLPRO(CHANNELS+1), &
     &            P(NTYP)%PSRMAX

      IC= P(NTYP)%NLPRO(CHANNELS+1)
      IF (CHANNELS+IC>LDIM) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'RD_PSEUDO: internal ERROR: increase LDIM'
         STOP
      ENDIF

      DO 615 L=CHANNELS+1,CHANNELS+IC
        P(NTYP)%LPS(L)   =P(NTYP)%LPS   (CHANNELS+1)
        P(NTYP)%NLPRO(L) =P(NTYP)%NLPRO (CHANNELS+1)
  615 CONTINUE

!-----Multipliers DION
      READ(10,*) &
     & ((P(NTYP)%DION(L,LP),L=CHANNELS+1,CHANNELS+IC),LP=CHANNELS+1,CHANNELS+IC)

      P(NTYP)%LMMAX=P(NTYP)%LMMAX+(2*P(NTYP)%LPS(CHANNELS+1)+1)*IC
      IF (P(NTYP)%LMMAX>LMDIM) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'MAIN: ERROR: increase LMDIM to ',P(NTYP)%LMMAX
         STOP
      ENDIF

      DO 630 NL=1,IC
!-----reciprocal projection operator
         READ(10,*)
         READ(10,*)(P(NTYP)%PSPNL  (I,CHANNELS+NL),I=1,NPSNL)

        IF (NWRITE>=0.AND.IU6>=0) &
            WRITE(IU6,*)' non local Contribution for L=', &
                        P(NTYP)%LPS(CHANNELS+1),' read in'
          IF (MOD(P(NTYP)%LPS(CHANNELS+1),2)==0) THEN
            P(NTYP)%PSPNL(0,CHANNELS+NL) =  P(NTYP)%PSPNL(2,CHANNELS+NL)
          ELSE
            P(NTYP)%PSPNL(0,CHANNELS+NL) = -P(NTYP)%PSPNL(2,CHANNELS+NL)
          ENDIF
!-----real space projection operator
        READ(10,*)
        READ(10,*)(P(NTYP)%PSPRNL (I,2,CHANNELS+NL),I=1,NPSRNL)

        IF (NWRITE>=0.AND.IU6>0) &
               WRITE(IU6,*)'   real space projection operators read in'

        DO I=1,NPSRNL
           P(NTYP)%PSPRNL (I,1,CHANNELS+NL)=(P(NTYP)%PSRMAX/NPSRNL)*(I-1)
        ENDDO
        CALL SPLCOF(P(NTYP)%PSPRNL (1,1,CHANNELS+NL) ,NPSRNL,NPSRNL,0._q)

  630 CONTINUE

      CHANNELS=CHANNELS+IC
      P(NTYP)%LMAX=CHANNELS
      GOTO 610

!-----------------------------------------------------------------------
!  read depletion charges
!  old dataset start with 'D' new (1._q,0._q) with tag 'A'
!-----------------------------------------------------------------------
  650 CONTINUE
   aug: IF ( CSEL /= 'P') THEN

      INFO%LOVERL=.TRUE.
      ALLOCATE(P(NTYP)%QDEP(NPSRNL,5,LDIM2),P(NTYP)%NDEP(LDIM,LDIM))

      IF (NWRITE>=0.AND.IU6>=0) &
           WRITE(IU6,*)'   augmentation charges read in'

      L2=1
      DO L =1 ,CHANNELS
      DO LP=L ,CHANNELS
      NREAD=0

  661   CONTINUE

        IF (L2> LDIM2) THEN
           IF (IU0>=0) &
           WRITE(IU0,*)'RD_PSEUDO: internal ERROR increase LDIM2'
           STOP
        ENDIF

        IF (CSEL=='D') THEN
         READ(10,*) IDUMM1,IDUMM2,P(NTYP)%QION(L,LP),P(NTYP)%PSDMAX
         P(NTYP)%NDEP(L,LP)=-1
        ELSE
         READ(10,*) IDUMM1,IDUMM2,IDUMM3,P(NTYP)%NDEP(L,LP), &
     &              P(NTYP)%QION(L,LP),P(NTYP)%PSDMAX
         IF (P(NTYP)%LPS(L)/=IDUMM1 .OR. P(NTYP)%LPS(LP)/=IDUMM2) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'ERROR: on reading POTCAR augmentation charges', &
            '       are wrong'
            STOP
         ENDIF
        ENDIF
        P(NTYP)%QION(LP,L) = P(NTYP)%QION(L,LP)

        IF (P(NTYP)%NDEP(L,LP)==0) GOTO 667

        READ(10,*) (P(NTYP)%QDEP(I,2,L2),I=1,NPSRNL)

        DO 665 I=1,NPSRNL
         P(NTYP)%QDEP(I,1,L2)=(I-1)*P(NTYP)%PSDMAX/NPSRNL
  665   CONTINUE
        ! changed  12Dec96 (gK)
        P(NTYP)%QDEP(NPSRNL,2,L2)=0 ! fix last value to (0._q,0._q)

        CALL SPLCOF(P(NTYP)%QDEP(1,1,L2),NPSRNL,NPSRNL,0._q)

        L2   =L2   +1
        NREAD=NREAD+1

  667   READ(10,*,ERR=600,END=600)
        IF (NREAD<ABS(P(NTYP)%NDEP(L,LP))) GOTO 661

      ENDDO
      ENDDO
      ELSE aug
!-----------------------------------------------------------------------
!  read paw data sets
!-----------------------------------------------------------------------
         INFO%LOVERL=.TRUE.
         LPAW=.TRUE.
         READ(10,*) NMAX,P(NTYP)%PSDMAX

         LMAX=0
         DO I=1,CHANNELS
            LMAX=MAX( P(NTYP)%LPS(I),LMAX )
         ENDDO
         LMAX=LMAX*2
         P(NTYP)%LMAX_CALC=LMAX

         READ(10,'(A)') FORM            ! format of remaining entities
         ALLOCATE( P(NTYP)%QPAW(CHANNELS,CHANNELS,0:LMAX),P(NTYP)%QATO(CHANNELS,CHANNELS), &
                   P(NTYP)%R%R(NMAX), P(NTYP)%POTAE(NMAX), P(NTYP)%POTPS(NMAX),  &
                   P(NTYP)%POTPSC(NMAX),P(NTYP)%TAUAE(NMAX), &
                   P(NTYP)%RHOAE(NMAX), P(NTYP)%RHOPS(NMAX),P(NTYP)%QTOT(CHANNELS,CHANNELS), &
                   P(NTYP)%WAE(NMAX,CHANNELS), P(NTYP)%WPS(NMAX,CHANNELS), &
                   P(NTYP)%NABLA(3,LMDIM,LMDIM))

! Intel efc compiler workaround (at least version 6.X)
!         P(NTYP)%NABLA=0
         DO I=1,LMDIM
            DO J=1,LMDIM
               P(NTYP)%NABLA(1,I,J)=0
               P(NTYP)%NABLA(2,I,J)=0
               P(NTYP)%NABLA(3,I,J)=0
            ENDDO
         ENDDO

         READ(10,*)
! Intel efc compiler workaround (at least version 6.X)
!         P(NTYP)%QPAW=0
         DO I=1,CHANNELS
            DO J=1,CHANNELS
               DO K=0,LMAX
                  P(NTYP)%QPAW(I,J,K)=0
               ENDDO
            ENDDO
         ENDDO

         READ(10,FORM) P(NTYP)%QPAW(:,:,0)

         READ(10,'(1X,A1)') CSEL
         IF (CSEL=='t') THEN
            READ(10,*) P(NTYP)%QTOT
            READ(10,*)
         ELSE
            DO I=1,CHANNELS
               DO J=1,CHANNELS
                  P(NTYP)%QTOT(I,J)=0
               ENDDO
               P(NTYP)%QTOT(I,I)=1
            ENDDO
         ENDIF
         READ(10,FORM) P(NTYP)%QATO

         CALL READGRD(10,FORM,P(NTYP)%R%R,'g',IU0)
         CALL READGRD(10,FORM,P(NTYP)%POTAE,'a',IU0)
         CALL READGRD(10,FORM,P(NTYP)%RHOAE,'c',IU0)
         READ(10,'(1X,A1)') CSEL
         IF (CSEL=='k') THEN
            CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%TAUAE,'k',CSEL,IU0)
            READ(10,'(1X,A1)') CSEL
         ELSE
            P(NTYP)%TAUAE=0
         ENDIF
         CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%POTPS,'p',CSEL,IU0)
         CALL READGRD(10,FORM,P(NTYP)%RHOPS,'c',IU0)

         P(NTYP)%R%NMAX  =NMAX
         P(NTYP)%R%RSTART=P(NTYP)%R%R(1)
         P(NTYP)%R%REND  =P(NTYP)%R%R(NMAX)
         P(NTYP)%R%D     =(P(NTYP)%R%REND/P(NTYP)%R%RSTART)**(1._q/(NMAX-1))
         P(NTYP)%R%H     =LOG(P(NTYP)%R%D)
         P(NTYP)%R%RMAX=P(NTYP)%PSDMAX
 
         CALL POTTORHO( P(NTYP)%ZVALF, NPSPTS, P(NTYP)%PSP(:,2), P(NTYP)%PSGMAX/NPSPTS, &
                   .TRUE. , NMAX, P(NTYP)%R%R ,  P(NTYP)%POTPSC )                        

         DO I=1,CHANNELS
            CALL READGRD(10,FORM,P(NTYP)%WPS(:,I),'p',IU0)
            CALL READGRD(10,FORM,P(NTYP)%WAE(:,I),'a',IU0)
         ENDDO
         READ(10,*,ERR=600,END=600) CSEL
         IF (NWRITE>=0.AND.IU6>=0) &
              WRITE(IU6,*)'   PAW grid and wavefunctions read in'
      ENDIF aug

  600 CONTINUE

      IF (NWRITE>=0.AND.IU6>=0) THEN
      WRITE(IU6,*)
      WRITE(IU6,*)'  number of l-projection  operators is LMAX  =', &
     &            P(NTYP)%LMAX
      WRITE(IU6,*)'  number of lm-projection operators is LMMAX =', &
     &            P(NTYP)%LMMAX
      WRITE(IU6,*)
      ENDIF

      ENDDO forever

  100 CONTINUE
      NTYP=NTYP-1

      CLOSE(10)
      RETURN
      END SUBROUTINE

!******************* SUBROUTINE READGRD *******************************
!
! small helper routine to read a grid based entry from the POTCAR
! file
!
!**********************************************************************

      SUBROUTINE READGRD(IUNIT,FORM,A,STRING,IU0)
        IMPLICIT NONE
        INTEGER  :: IUNIT,IU0
        CHARACTER*(*) :: FORM
        REAL(q)  :: A(:)
        CHARACTER*(1) :: STRING,CSEL

        READ(IUNIT,'(1X,A1)') CSEL
        IF (CSEL /= STRING) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'READGRD: POTCAR file has an error: '
              WRITE(IU0,*)' expected: ',STRING,' found: ',CSEL
           ENDIF
           STOP
        ENDIF

        READ(IUNIT,FORM) A

      END SUBROUTINE READGRD

      ! selector already read, just compare CSEL and STRING

      SUBROUTINE READGRD_NO_SELECTOR(IUNIT,FORM,A,STRING,CSEL,IU0)
        IMPLICIT NONE
        INTEGER  :: IUNIT,IU0
        CHARACTER*(*) :: FORM
        REAL(q)  :: A(:)
        CHARACTER*(1) :: STRING,CSEL

        IF (CSEL /= STRING) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'READGRD: POTCAR file has an error: '
              WRITE(IU0,*)' expected: ',STRING,' found: ',CSEL
           ENDIF
           STOP
        ENDIF

        READ(IUNIT,FORM) A

      END SUBROUTINE READGRD_NO_SELECTOR

!******************* SUBROUTINE DEALLOC_PP ****************************
!
!  deallocate PP arrays
!  (I always try to do this in reverse order,
!   might help the memory subsystem to remerge free junks)
!**********************************************************************
      SUBROUTINE DEALLOC_PP(P,NTYPD)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (potcar) P(NTYPD)

      DO NTYP=NTYPD,1,-1
        CALL DEALLOC_PP1(P(NTYP))
      ENDDO
      RETURN

      END SUBROUTINE



      SUBROUTINE DEALLOC_PP1(P)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (potcar) P

      IF (ASSOCIATED(P%QDEP)) THEN
        DEALLOCATE(P%QDEP,P%NDEP)
      ENDIF
      DEALLOCATE(P%DION,P%QION, &
     &         P%LPS,P%NLPRO, &
     &         P%PSPNL,P%PSPRNL)

      DEALLOCATE(P%PSPRHO)

      IF (ASSOCIATED(P%PSPCOR)) THEN
        DEALLOCATE(P%PSPCOR)
      ENDIF
      DEALLOCATE(P%PSP)

      RETURN
      END SUBROUTINE

!******************* SUBROUTINE RDPARA ********************************
!
!  this subroutine  reads the implicit PSCTR description from the
!  POTCAR file (UNIT=IU)
!  information is stored in the string array
!   STRING(ISDIM)
!  the number of lines is returned in IREAD
!
!**********************************************************************

      SUBROUTINE RDPARA(IU,ISDIM,STRING,IREAD,IWRITE)
      USE prec
      USE mpimy
      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER*80 STRING(ISDIM)
      CHARACTER*6  TAG

      I=0
  100 I=I+1
      IF (I>ISDIM) THEN
        WRITE(*,*)'ERROR: Description of pseudopotential is too long'
        STOP
      ENDIF
      READ(IU,'(A)') STRING(I)
      IF (STRING(I)(1:3)=='END') GOTO 120
      GOTO 100
  120 CONTINUE
      IREAD=I

      DO 200 I=1,IREAD
        TAG=STRING(I)(4:9)
        IF (TAG(1:5)=='Error') GOTO 210
  200 CONTINUE
  210 IWRITE=I-1

      RETURN
      END SUBROUTINE

!******************* SUBROUTINE RDPARS ********************************
!
!  this subroutine interprets the STRING array read in the
!  subroutine RDPARA
!  only a few number of items are retreeved from the STRING array
!   LREAL
!
!**********************************************************************


      SUBROUTINE RDPARS(ISDIM,STRING,IREAD,P,IU6)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P
      CHARACTER*80 STRING(ISDIM)

      CHARACTER*6 TAG
      CHARACTER*80 VALUE
      LOGICAL LDUM
      EXTERNAL  LENGTH
!
!     initialize all values
!
      P%LREAL=.FALSE.
      P%LEXCH=0
      P%ZVALF=0
      P%POMASS=0
      P%RWIGS=0
      P%EATOM=0
      P%ENMAXA=0
      P%ENMINA=0
      P%QOPT1=0
      P%QOPT2=0
      P%EAUG=0
      P%EKECUT=0
      P%EKEERR=0
      P%EGGA=0
      P%DEXCCORE=0
!-MM- added for use in MAGCAR
      P%ELEMENT='  '
!-MM- end of addition
      

      I=1
!  next line
  120 CONTINUE

      ITEM=1
!  next item in line
  110 CONTINUE

      IF (ITEM==1) THEN
        TAG=STRING(I)(4:9)
        VALUE=STRING(I)(13:80)
      ELSE
        TAG=STRING(I)(23:28)
        VALUE=STRING(I)(32:80)
      ENDIF

      EDUML=0
      EDUM=0
      IDUM=0
      LDUM=.FALSE.

      READ(VALUE,'(G10.4)',IOSTAT=IERR) EDUML
      READ(VALUE,'(G8.3)',IOSTAT=IERR) EDUM
      READ(VALUE,'(L8)',IOSTAT=IERR) LDUM
      READ(VALUE,'(I8)',IOSTAT=IERR) IDUM
!
!     set values
!
      L=LENGTH(TAG)
      IF (TAG(1:L)=='QCUT')   THEN
        P%LREAL=.TRUE.
        P%QOPT1=EDUM
        P%QOPT2=EDUM*2
      ENDIF
      IF (TAG(1:L)=='QGAM')   THEN
        P%LREAL=.TRUE.
        P%QOPT2=EDUM
      ENDIF
      IF (TAG(1:L)=='ZVAL')   THEN
         P%ZVALF  =EDUM
         P%ZVALF_ORIG=EDUM
      ENDIF
      IF (TAG(1:L)=='POMASS') P%POMASS=EDUM
      IF (TAG(1:L)=='RWIGS')  P%RWIGS =EDUM
      IF (TAG(1:L)=='ENMAX')  THEN
        P%ENMAXA=EDUM
      ENDIF
      IF (TAG(1:L)=='ENMIN')  P%ENMINA=EDUM
      IF (TAG(1:L)=='EAUG')   P%EAUG  =EDUM
      IF (TAG(1:L)=='EATOM')  P%EATOM =EDUML
      IF (TAG(1:L)=='DEXC')   P%DEXCCORE=EDUM
      IF (TAG(1:L)=='LEXCH')  THEN
        CALL EXTYP(VALUE(1:2),P%LEXCH)
      ENDIF
      IF (TAG(1:L)=='GGA')  THEN
        READ(STRING(I)(13:80),'(4G10.4)',IOSTAT=IERR) P%EGGA
      ENDIF
!  kinetic energy error lines
      IF (TAG(1:L)=='Error')  THEN
        NLINE=8
        I=I+1
        READ(STRING(I)(13:80),'(I8)',IOSTAT=IERR)   NDATA
        I=I+1
        READ(STRING(I)(13:80),'(2G8.3)',IOSTAT=IERR) START,STEP
        IF (NDATA/=NEKERR) GOTO 130
        DO I1=0,NDATA-1,NLINE
        I=I+1
        IM=MIN(NDATA-I1,NLINE)
        READ(STRING(I),'(8G10.4)',IOSTAT=IERR) (P%EKEERR(I1+I2),I2=1,IM)
        DO I2=1,IM
          P%EKECUT(I1+I2)=STEP**(I1+I2-1)*START
        ENDDO
        ENDDO
!  next item
      ENDIF
!-MM- for use in MAGCAR
      IF (TAG(1:L)=='VRHFIN')  THEN
         P%ELEMENT(1:1)=STRING(I)(12:12)
         IF (STRING(I)(13:13).NE.':') P%ELEMENT(2:2)=STRING(I)(13:13)
      ENDIF
!-MM- end of addition
      ITEM=ITEM+1
      IF (ITEM<=2)  GOTO 110
      I=I+1
      IF (I<=IREAD) GOTO 120
      RETURN

!   error in kinetic energy lines
  130 CONTINUE
      IF (IU6>=0) &
      WRITE(*,*)'WARNING: can not read POTCAR file (kinetic energy)'

      END SUBROUTINE

!******************* SUBROUTINE EXTYP *********************************
!
!  this subroutine interprets the string CEXCH
!  which determines the type of exchange correlation
!  and sets the integer LEXCH accordingly
!
!**********************************************************************

      SUBROUTINE EXTYP(CEXCH,LEXCH)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER*2 CEXCH
      LEXCH=-1
        IF (CEXCH=='  ') THEN
          LEXCH=0
        ELSE IF (CEXCH=='HL') THEN
          LEXCH=1
        ELSE IF (CEXCH=='PZ') THEN
          LEXCH=2
        ELSE IF (CEXCH=='CA') THEN
          LEXCH=2
        ELSE IF (CEXCH=='WI') THEN
          LEXCH=3
        ELSE IF (CEXCH=='PB') THEN
          LEXCH=4
        ELSE IF (CEXCH=='PW') THEN
          LEXCH=5
        ELSE IF (CEXCH=='LM') THEN
          LEXCH=6
        ELSE IF (CEXCH=='91') THEN
          LEXCH=7
        ELSE IF (CEXCH=='PE') THEN
          LEXCH=8
        ELSE IF (CEXCH=='RP') THEN
          LEXCH=9
        ENDIF
      RETURN
      END SUBROUTINE


!******************* SUBROUTINE PSPOST ********************************
!
!  postprocessing of pseudopotential reader
!  checks
!
!**********************************************************************

      SUBROUTINE POST_PSEUDO(NTYPD,NTYP_PP,NTYP,NIONS,NITYP,P,INFO, &
     &        LREALD,ROPT,IDIOT,IU6,IU0,LMAX_CALC,L_NO_US)
      USE prec
      USE base
      USE  ini
      USE constant

      IMPLICIT NONE

      INTEGER NTYPD
      TYPE (potcar) P(NTYPD)
      TYPE (info_struct) INFO
      INTEGER NIONS         ! number of ions
      INTEGER NTYP,NTYP_PP  ! number of spec. on INCAR  / POTCAR
      INTEGER NITYP(NTYPD)  ! numer of ions of each species
      LOGICAL LREALD        ! default for LREAL from INCAR
      REAL(q)    ROPT(NTYPD)   ! cutoff for automatic real space opt
      LOGICAL L_NO_US       ! no US PP
      INTEGER LMAX_CALC     ! maximum L quantum number in PAW method
! arryas for tutor call
      INTEGER IU6,IU0,IDIOT,ITUT(3),CDUM,LDUM
      REAL(q)    RTUT(3)
! temporary variables
      LOGICAL LREALT
      INTEGER NT,LOXCH
      REAL(q) NVALEL
      REAL(q) QMAXL,QMAXNL,QMAXNL2,EERROR,ENAUG2,DUMMY
      REAL(q) SPLFIT(NEKERR,5)

!-----------------------------------------------------------------------
! check number of species
!-----------------------------------------------------------------------
      IF (NTYP/=NTYP_PP) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'ERROR: number of potentials on File POTCAR', &
     &              ' incompatible with number of species', &
     &             'INCAR :',NTYP,'POTCAR: ',NTYP_PP
        STOP
      ENDIF
!-----------------------------------------------------------------------
! all PP generated with the same XC-type
! also find out whether all PP are real space optimized
!-----------------------------------------------------------------------
      INFO%LEXCH=P(1)%LEXCH
      LREALT=P(1)%LREAL
      DO NT=2,NTYP
         IF (INFO%LEXCH /= P(NT)%LEXCH) THEN
         ITUT(1)=P(NT)%LEXCH
         ITUT(2)=INFO%LEXCH
         ITUT(3)=NT
         CALL VTUTOR('E','DIFFERENT LDA-XC TYPES',RTUT,1, &
     &                ITUT,3,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('S','DIFFERENT LDA-XC TYPES',RTUT,1, &
     &               ITUT,3,CDUM,1,LDUM,1,IU0,IDIOT)
        ENDIF
        LREALT=LREALT.AND.P(1)%LREAL
      ENDDO
!-----------------------------------------------------------------------
!  set default value for LREAL if not given in INCAR
!-----------------------------------------------------------------------
      IF (LREALD) INFO%LREAL=.FALSE.
!---- complain about using non-optimized projectors if LREAL=.T.
      IF (INFO%LREAL.AND.(.NOT.LREALT)) THEN
         CALL VTUTOR('W','REAL-SPACE WITHOUT OPTIMIZATION', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('W','REAL-SPACE WITHOUT OPTIMIZATION', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!---- 'very large' cell and no real-space-projection scheme ???
      IF ((.NOT.INFO%LREAL).AND.(NIONS>16)) THEN
         IF (LREALT) THEN
            CALL VTUTOR('A','NO REAL-SPACE AND YOU COULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
            CALL VTUTOR('A','NO REAL-SPACE AND YOU COULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
         ELSE
            CALL VTUTOR('A','NO REAL-SPACE AND YOU SHOULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
            CALL VTUTOR('A','NO REAL-SPACE AND YOU SHOULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
         ENDIF
      ENDIF
!---- 'very small' cell and still real-space projection scheme ???
      IF (INFO%LREAL.AND.(NIONS<=8)) THEN
         CALL VTUTOR('A','REAL-SPACE NOMORE RECOMMENDED', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','REAL-SPACE NOMORE RECOMMENDED', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!-----------------------------------------------------------------------
!  set energy cutoff and
!  check whether range of potentials is sufficient
!-----------------------------------------------------------------------
      IF (INFO%ENMAX==-1) THEN
      DO NT=1,NTYP
        IF (P(NT)%ENMAXA==0 &
     &   .OR.INFO%SZPREC(1:1)=='l'.AND.P(NT)%ENMINA==0 ) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Fatal error! Could not find entry for ENMAX'// &
     &               ' on file INCAR. MUST be specified'
         STOP
        ENDIF
        IF  (INFO%SZPREC(1:1)=='l') THEN
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMINA)
        ELSE IF (INFO%SZPREC(1:1)=='h') THEN
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMAXA*1.25_q)
          IF (IU0>0) THEN
             WRITE(0,*) 'WARNING: for PREC=h ENMAX is automatically increase by 25 %'
             WRITE(0,*) '       this was not the case for versions prior to vasp.4.4'
          ENDIF
        ELSE
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMAXA)
        ENDIF
      ENDDO
      ENDIF
      IF (INFO%ENINI==-1) INFO%ENINI=INFO%ENMAX

      QMAXNL=   SQRT(INFO%ENMAX /RYTOEV)/AUTOA
      QMAXL = 2*SQRT(INFO%ENMAX /RYTOEV)/AUTOA

      DO NT=1,NTYP
      IF (ROPT(NT)/=0) THEN
        QMAXNL2    =QMAXNL*2.0
        IF  (INFO%SZPREC(1:1)=='h') THEN
           QMAXNL2    =QMAXNL*2.5
        ENDIF

        IF ( ABS(ROPT(NT)) > 0.1) THEN
          P(NT)%PSRMAX=12.5_q/ABS(QMAXNL)*ABS(ROPT(NT))**(1._q/3._q)
          IF  (INFO%SZPREC(1:1)=='h') &
          P(NT)%PSRMAX=10/ABS(QMAXNL)*ABS(ROPT(NT))**(1._q/3._q)
        ELSE
          P(NT)%PSRMAX=12.5_q/ABS(QMAXNL)*(0.5)**(1._q/3._q)
          IF  (INFO%SZPREC(1:1)=='h') &
          P(NT)%PSRMAX=10/ABS(QMAXNL)*(0.5)**(1._q/3._q)
        ENDIF

	IF (ROPT(NT)>0) THEN
        CALL OPTREAL(IU6,P(NT)%LMAX,P(NT)%LPS(1),NPSRNL, &
     &         P(NT)%PSMAXN/NPSRNL, &
     &         P(NT)%PSPRNL (1,1,1),P(NT)%PSPNL(0,1), &
     &         QMAXNL2,QMAXNL,P(NT)%PSRMAX,ROPT(NT)<0.1,ROPT(NT))
        ELSE
        CALL OPTREAL_NEW(IU6,P(NT)%LMAX,P(NT)%LPS(1),NPSRNL, &
     &         P(NT)%PSMAXN/NPSRNL, &
     &         P(NT)%PSPRNL (1,1,1),P(NT)%PSPNL(0,1), &
     &         QMAXNL2,QMAXNL,P(NT)%PSRMAX,ABS(ROPT(NT))<0.1,ABS(ROPT(NT)))

        ENDIF
        P(NT)%QOPT1=QMAXNL *AUTOA
        P(NT)%QOPT2=QMAXNL2*AUTOA
      ENDIF

      IF (P(NT)%PSGMAX< QMAXL ) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WARNING: PSGMAX for local potential too small'
        IF (IU6>=0) THEN
        WRITE(IU6,*)'WARNING: PSGMAX for local potential too small'
        WRITE(IU6,*)'         PSGMAX should be >',QMAXL,' NTYP=',NT
        ENDIF
      ENDIF
      IF (P(NT)%PSMAXN<  QMAXNL ) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WARNING: PSMAXN for non-local potential too small'
        IF (IU6>=0) THEN
        WRITE(IU6,*)'WARNING: PSMAXN for non-local potential too small'
        WRITE(IU6,*)'         PSMAXN should be >', QMAXNL,' NTYP=',NT
        ENDIF
      ENDIF
      IF (INFO%LREAL.AND. P(NT)%QOPT1/=0) THEN
        IF ((P(NT)%QOPT1**2*RYTOEV-INFO%ENMAX)/INFO%ENMAX &
     &      > 0.1_q) THEN
         RTUT(1)=P(NT)%QOPT1**2*RYTOEV
         RTUT(2)=INFO%ENMAX
         RTUT(3)=SQRT(INFO%ENMAX/RYTOEV)
         CALL VTUTOR('A','WRONG OPTIMZATION REAL-SPACE', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','WRONG OPTIMZATION REAL-SPACE', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
        ENDIF
      ENDIF
      ENDDO
!-----------------------------------------------------------------------
! switch to specific LDA/GGA if tag was found on INCAR file
!-----------------------------------------------------------------------
      LOXCH=INFO%LEXCH
      CALL EXTYP(INFO%SZGGA,INFO%LEXCH)
      IF (INFO%LEXCH<0) INFO%LEXCH=LOXCH

      INFO%LEXCHG = 0

      IF (INFO%LEXCH>=4) THEN
         INFO%LEXCHG=INFO%LEXCH-3
!  add GGA corrections to EATOM
         DO NT=1,NTYP
           P(NT)%EATOM=P(NT)%EATOM-P(NT)%EGGA(MOD(INFO%LEXCHG-1,4)+1)
         ENDDO
      ENDIF
      IF (INFO%LEXCH/=LOXCH) THEN
         CALL VTUTOR('A','ENFORCED LDA', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','ENFORCED LDA', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!-----------------------------------------------------------------------
!  set number of electrons and total atomic energy
!-----------------------------------------------------------------------
      INFO%EALLAT=0
      NVALEL=0
      ENAUG2=0
      DO NT=1,NTYP
       ENAUG2=MAX(P(NT)%EAUG,ENAUG2)
       EERROR=0
       IF (IU6>=0) &
       WRITE(IU6,449) P(NT)%SZNAMP,NT,-P(NT)%EATOM

  449  FORMAT( A40,':'/ &
     & ' energy of atom ',I2,'       EATOM=',F10.4)
       IF ((P(NT)%EKEERR(1)/=0 .OR. P(NT)%EKECUT(1)/=0)) &
     & THEN
       CALL SPLCPY(P(NT)%EKECUT,P(NT)%EKEERR,SPLFIT,NEKERR,NEKERR,10E30_q)
       CALL SPLVAL(INFO%ENMAX,EERROR,DUMMY,SPLFIT,NEKERR,NEKERR)
       IF (IU6>=0) &
       WRITE(IU6,448) EERROR

  448  FORMAT( &
     & ' kinetic energy error for atom=',F10.4, &
     & ' (will be added to EATOM!!)')
       ENDIF
       INFO%EALLAT=INFO%EALLAT+NITYP(NT)*(P(NT)%EATOM-EERROR)
       NVALEL=NVALEL+NITYP(NT)*P(NT)%ZVALF
      ENDDO
      IF (INFO%ENAUG==-1) INFO%ENAUG=ENAUG2

      IF (INFO%NELECT==0) INFO%NELECT= NVALEL
      IF (IU6>=0) WRITE(IU6,*)
!-----------------------------------------------------------------------
!  set truncation for LMAX for PAW
!-----------------------------------------------------------------------
      IF (LMAX_CALC >= 0) THEN
         P%LMAX_CALC=LMAX_CALC
      ENDIF

      L_NO_US=.TRUE.
      DO NT=1,NTYP
         IF (ASSOCIATED( P(NT)%QDEP )) THEN
            L_NO_US=.FALSE.
         ENDIF
      ENDDO

      END SUBROUTINE

!**********************************************************************
!
!  postprocessing of pseudopotential reader
!  checks
!
!**********************************************************************

      SUBROUTINE LDIM_PSEUDO(LORBIT, NTYPD, P, LDIMP, LMDIMP)

      USE prec
      USE base
      USE  ini
      USE constant

      IMPLICIT NONE

      INTEGER LORBIT       ! how do calculate partial dos
      INTEGER NTYPD        ! how many types are there
      INTEGER LDIMP,LMDIMP ! required dimension of arrays
      TYPE (potcar) P(NTYPD)
! local
      INTEGER MAXIMUM_L,NTYP,I

      MAXIMUM_L=0

      DO NTYP=1,NTYPD
         DO I=1,P(NTYP)%LMAX
            MAXIMUM_L = MAX(MAXIMUM_L,P(NTYP)%LPS(I))
         ENDDO
      ENDDO
      
      ! if the maximum L is larger than 2 we have to set LDIMP

      IF ( MAXIMUM_L > 2 ) THEN
         LDIMP= MAXIMUM_L+1
      ELSE
         LDIMP=3
      ENDIF
      LMDIMP=LDIMP*LDIMP

      END SUBROUTINE LDIM_PSEUDO

      END MODULE
