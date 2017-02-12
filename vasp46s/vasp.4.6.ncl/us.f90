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





      MODULE us
      USE prec
!***********************************************************************
! RCS:  $Id: us.F,v 1.6 2003/06/27 13:22:23 kresse Exp kresse $
!
! a few comments concerning implementation:
! It would be simplest to have assumed shape arrays for CDIJ...
! But I have not done this up do now.
! One problem (which I really hate) is that VASP uses complex
! to real FFT. This requires that the same array is in some
! places addressed as a real, and in other places as a complex
! array. Because no cast exists in f90, I have to use the
! strange INTERFACE construct below to get at least some error
! checking.
! This fact also keeps me from using assumed shape arrays (gK)
!***********************************************************************
!***********************************************************************
!
! specify the interface for SETDIJ and DEPLE
! the dummy argument CVTOT differs here from the actual definition
! in the subroutine used below
!
!***********************************************************************

      CONTAINS


!************************ SUBROUTINE DEPSUM ****************************
!
! this subroutine calculates  the total "occupancy" of each
! ll'mm'augmentation charge from the  FERMI weights, the weights of
! each k-point WTKPT and  projection operators acting on all bands
! result is stored in CRHODE (10.32)
!
!***********************************************************************

      SUBROUTINE DEPSUM(W,WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      COMPLEX(q) CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)


      IF (.NOT.LOVERL) RETURN
!=======================================================================
! initialise to (0._q,0._q)
!=======================================================================
      CRHODE=0
!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      band:   DO N=1,WDES%NBANDS

      WEIGHT=WDES%RSPIN*W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *WDES%NPRO/2
      LMBASE_=ISPINOR_*WDES%NPRO/2

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
              WEIGHT*W%CPROJ(L+LMBASE,N,NK,ISP)*CONJG(W%CPROJ(LP+LMBASE_,N,NK,ISP))
        ENDDO
        ENDDO
	
      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
      


      RETURN
      END SUBROUTINE



!************************ SUBROUTINE DEPATO ****************************
!
! set the "occupancy" of each llp channel CRHODE to the values during
! the pseudopotential generation
!
! in the spinpolarized collinear case the CRHODE(:,2) is set
! to CRHODE(:,1) / ZVAL * ATOMOM
! in the non collinear case
! the magnetization in each direction x,y and z is set in the way
!
!***********************************************************************

      SUBROUTINE DEPATO(WDES, LMDIM, CRHODE, LOVERL, P, T_INFO)
      USE prec
      USE pseudo
      USE wave
      USE poscar
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      LOGICAL LOVERL
      INTEGER LMDIM
      COMPLEX(q) CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! local variables
      INTEGER NI,NIP,LOW,LM,NT,LL,LHI,MMAX,L,LP,M,ISP
      REAL(q) :: VALUE(WDES%NCDIJ)

      IF (.NOT.LOVERL) RETURN
      CRHODE=0

  ion: DO NI=1,T_INFO%NIONS
      NIP=NI_LOCAL(NI, WDES%COMM_INB)
      IF (NIP==0) CYCLE ion ! not on local node

      LOW=1
      LM =1
      NT=T_INFO%ITYP(NI)

      IF (.NOT.ASSOCIATED(P(NT)%QATO)) CYCLE ion
      
      !
      ! this little trick does push the PAW occupancies in the right
      ! direction for core level shift calculations
      !
      VALUE(1) = P(NT)%ZVALF/P(NT)%ZVALF_ORIG
      IF (WDES%NCDIJ==2) THEN
         VALUE(2)= T_INFO%ATOMOM(NI)/P(NT)%ZVALF
      ELSE IF (WDES%NCDIJ == 4) THEN
         VALUE(2)= T_INFO%ATOMOM(1+(NI-1)*3)/P(NT)%ZVALF
         VALUE(3)= T_INFO%ATOMOM(2+(NI-1)*3)/P(NT)%ZVALF
         VALUE(4)= T_INFO%ATOMOM(3+(NI-1)*3)/P(NT)%ZVALF
      ENDIF

      block: DO
         LL=P(NT)%LPS(LOW)
         ! search block with same L
         DO LHI=LOW,P(NT)%LMAX
            IF (LL/=P(NT)%LPS(LHI)) EXIT
         ENDDO
         LHI=LHI-1
         MMAX=2*LL+1

         DO ISP=1,WDES%NCDIJ

         DO L =LOW,LHI
         DO LP=LOW,LHI
         DO M =0,MMAX-1
            CRHODE(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)= &
                    P(NT)%QATO(L,LP)*VALUE(ISP)
         ENDDO
         ENDDO
         ENDDO

         ENDDO

         ! set new LOW value and LM value and go on
         LM=LM+(LHI-LOW+1)*MMAX
         LOW=LHI+1
         IF (LOW > P(NT)%LMAX) EXIT block
      ENDDO block

      ENDDO ion
      RETURN
      END SUBROUTINE


!************************* SUBROUTINE FORDEP ***************************
!
! this subroutine calculates the forces, which are related to the
! the change in the position of the augmentation charges in a fixed
! local potential
! finite central differences are usually used
! (only if noHighPrec is specified the old version with 1. order
!  finite differences is used)
! the accuracy of the routine is at least 7 significant digits
! previously it was 4-5 digits
!***********************************************************************

      SUBROUTINE FORDEP(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
       LMDIM,CDIJ,CQIJ,CRHODE, CVTOT, IRDMAX, FORNL)
      USE prec

      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      COMPLEX(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      REAL(q)     FORNL(3,T_INFO%NIONS)
      LOGICAL  LOVERL
! work array
      DIMENSION DISPL(3)
      DIMENSION EDEP(T_INFO%NIONS)
      REAL(q)   FORAUG(3,T_INFO%NIONS)

      IF (.NOT. LOVERL) RETURN

      ! we need the spinor representation of CRHODE at this point
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .TRUE.)

      DIS=1E-5_q
!TEST
!     WRITE(*,'(4E14.7)') DIS
!     DIS=1E-2
!1000 DIS=DIS/2
!TEST
      FORAUG=0
!=======================================================================
! calculate the contribution to the energy from the augmentation-hole
! for displacement 0
! (1st order finite difference formula)
!=======================================================================
      DISPL=0
!=======================================================================
! calculate the contribution to the energy from the augmentation charges
! for displacement X
! (and displacement -X if central differences are used)
!=======================================================================
      dir: DO IDIR=1,3
      DISPL=0
      DISPL(IDIR)=-DIS/2
      CALL SETDIJ(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX,DISPL(1),DISPL(2),DISPL(3))

      NIS=1
      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
          CALL FORDE1(LMDIM,WDES%NIONS,WDES%NCDIJ,LMMAXC,NI,CDIJ,CRHODE,ADD)
          EDEP(NI)=ADD
        ENDDO
      NIS = NIS+WDES%NITYP(NT)
      ENDDO

      DISPL=0
      DISPL(IDIR)=DIS/2

      CALL SETDIJ(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX,DISPL(1),DISPL(2),DISPL(3))

      NIS=1
      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
          CALL FORDE1(LMDIM,WDES%NIONS,WDES%NCDIJ,LMMAXC,NI,CDIJ,CRHODE,ADD)
          NIP=NI_GLOBAL(NI, WDES%COMM_INB)     !  local storage index
          FORAUG(IDIR,NIP)=FORAUG(IDIR,NIP)-(ADD-EDEP(NI))/DIS
        ENDDO
        NIS = NIS+WDES%NITYP(NT)
      ENDDO

      ENDDO dir

      

!TEST
!     WRITE(*,'(4E19.12)') DIS,FORAUG(1:3,1:1)
!     GOTO 1000
!TEST
      FORNL=FORNL+FORAUG

      ! we need the spinor representation of CRHODE
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .FALSE.)

      RETURN
      END SUBROUTINE

!
! small helper routine
!
      SUBROUTINE FORDE1(LMDIM,NIONS,ISPIN,LMMAXC,NI,CDIJ,CRHODE,ADD)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      COMPLEX(q) CRHODE(LMDIM,LMDIM,NIONS,ISPIN)
      COMPLEX(q) CDIJ(LMDIM,LMDIM,NIONS,ISPIN)

      ADD=0

      DO ISP=1,ISPIN
      DO L=1,LMMAXC
      DO LP=1,LMMAXC
        ADD=ADD+ CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP))
      ENDDO; ENDDO; ENDDO

      RETURN

      END SUBROUTINE


!************************* SUBROUTINE STRDEP ***************************
!
!  subroutine for calculating the contributions to the stress,
!  related to the augmentation charges
!  uncomment CTEST-lines if you want to test finite differences
!
!***********************************************************************

      SUBROUTINE STRDEP(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
       LMDIM,CDIJ,CQIJ,CRHODE, CVTOT, IRDMAX, ISIF,AUGSIF)
      USE prec

      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR,LATT_FIN

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      COMPLEX(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      LOGICAL  LOVERL
      REAL(q)     AUGSIF(3,3)

      IF (.NOT. LOVERL) RETURN

      ! we need the spinor representation of CRHODE at this point
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .TRUE.)

      DIS=1E-7_q
!TEST
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
!=======================================================================
! calculate the contribution to the energy from the augmentation-hole
! for undistorted lattice
!=======================================================================
      AUGSIF=0
      CALL LATTIC(LATT_CUR)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX,.0_q,.0_q,.0_q)

      EAUG=0
      NIS=1

      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
        DO ISP=1,WDES%NCDIJ
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
          EAUG=EAUG+ REAL( CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP)),KIND=q)
        ENDDO; ENDDO; ENDDO; ENDDO

        NIS = NIS+WDES%NITYP(NT)
      ENDDO
!=======================================================================
! calculate the contribution to the energy from the augmentation charges
! for distortion X
! 1. order finite differences are used
! (for stress high precision is not really required)
!=======================================================================
      DO IDIR=1,3
      DO JDIR=1,3

      LATT_FIN=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3; DO J=1,3
          LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
        ENDDO; ENDDO
      ELSE
!  all directions
        DO I=1,3
          LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
        ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_FIN,P,T_INFO,LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX,.0_q,.0_q,.0_q)

      EAUGD=0
      NIS=1
      DO  NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
        DO ISP=1,WDES%NCDIJ
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
          EAUGD=EAUGD+ REAL( CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP)),KIND=q)
        ENDDO; ENDDO; ENDDO; ENDDO

        NIS = NIS+WDES%NITYP(NT)
      ENDDO

      AUGSIF(IDIR,JDIR)=EAUG-EAUGD
!
!  only isotrop pressure terminate loop (IDIR=1 and JDIR=1)
!
      IF (ISIF==1) THEN
        AUGSIF(2,2)= AUGSIF(1,1)
        AUGSIF(3,3)= AUGSIF(1,1)
        GOTO 400 ! terminate (not very clean but who cares)
      ENDIF

      ENDDO
      ENDDO
!=======================================================================
! calculation finished  scale stress
!=======================================================================
  400 CONTINUE
      

      AUGSIF=AUGSIF/DIS

!TEST
!      WRITE(*,'(E10.3,"  ",7E14.7)')DIS,AUGSIF
!      IF (DIS>1E-10) GOTO 1000
!TEST
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .FALSE.)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE CHARGE     ************************
!
! calculate the charge density on the soft and fine grid
! results are returned in the convention
! (total, magnetization (x,y,z))
!
!***********************************************************************

      SUBROUTINE SET_CHARGE(W, WUP, WDW, WDES, LOVERL, &
                  GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
                  LATT_CUR, P, SYMM, T_INFO, &
                  CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      USE paw
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave

      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT,GRIDUS
      TYPE (transit)     C_TO_US,SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM
      TYPE (wavefun)     WUP,WDW
      TYPE (wavespin)    W

      INTEGER LMDIM
      INTEGER IRDMAX
      INTEGER N_MIX_PAW
      REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
! on return the following arrays are set
      COMPLEX(q)     CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) ! occupancy of augmentation channels
      COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)         ! soft pseudo charge
      COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ)             ! total charge
      LOGICAL LOVERL
      INTEGER ISP

!---  update of charge if necessary
      CALL SOFT_CHARGE(GRID,GRID_SOFT,W, WDES, CHDEN)
!     change storage convention to (total, magnetization)
      CALL RC_FLIP(CHDEN,GRID_SOFT,WDES%NCDIJ,.FALSE.)

      IF (SYMM%ISYM ==2) THEN
         IF (WDES%LNONCOLLINEAR) THEN
            CALL RHOSYM(CHDEN(1,1),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
! Marsman: Insert symmetrization of vector field
            IF (.NOT.WDES%LSPIRAL) &
           &   CALL SYMFIELD(CHDEN(1,2),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
         ELSE
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(CHDEN(1,ISP),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
            ENDDO
         ENDIF
      ENDIF

      CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                 LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
                 LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)

      CALL SET_RHO_PAW(WDES, P, T_INFO, LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)

      IF (SYMM%ISYM ==1) THEN
         IF (WDES%LNONCOLLINEAR) THEN
            CALL RHOSYM(CHTOT(1,1),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
! Marsman: Insert symmetrization of vector field
            IF (.NOT.WDES%LSPIRAL) &
           &   CALL SYMFIELD(CHTOT(1,2),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
         ELSE
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(CHTOT(1,ISP),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
            ENDDO
         ENDIF
      ENDIF

      END SUBROUTINE SET_CHARGE

END MODULE

!************************ SUBROUTINE SETYLM_AUG ************************
!
! this subroutine performes the following tasks
! ) finds the points, which are within a certain cutoff around (1._q,0._q) ion
! ) calculates the distance of each of this points from the ion
! ) calculates the spherical harmonics Y_lm(Omega(r-R(ion))
! DISX,Y,Z are additional displacements of the ions
!
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
!
!***********************************************************************

      SUBROUTINE SETYLM_AUG(GRID,LATT_CUR,POSION,PSDMAX,NPSRNL, &
     &        LMYDIM,LYDIM,YLM,IRMAX,INDMAX,DISX,DISY,DISZ,DIST,NLI, &
     &         XS,YS,ZS)
      USE prec

      USE mpimy
      USE mgrid
      USE lattice
      USE asa
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR

      DIMENSION POSION(3)

      DIMENSION DIST(IRMAX)
      DIMENSION YLM(IRMAX,LMYDIM)
! work-arrays
      DIMENSION NLI(IRMAX)
!-MM- changes to accomodate spin spirals
!     REAL(q), ALLOCATABLE :: XS(:),YS(:),ZS(:)
!     ALLOCATE(XS(IRMAX),YS(IRMAX),ZS(IRMAX))
      REAL(q) XS(IRMAX),YS(IRMAX),ZS(IRMAX)
      XS=0;YS=0;ZS=0
!-MM- end of alterations
      IF ((LYDIM+1)**2 > LMYDIM) THEN
         WRITE(0,*)'internal error: LMYDIM is too small',LYDIM,LMYDIM
         STOP
      ENDIF
!=======================================================================
! find lattice points contained within the cutoff-sphere
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
! which is a somewhat strange convention
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

      ARGSC=NPSRNL/PSDMAX
!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= PSDMAX*LATT_CUR%BNORM(1)*GRID%NGX
      D2= PSDMAX*LATT_CUR%BNORM(2)*GRID%NGY
      D3= PSDMAX*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(POSION(3)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
      N2LOW= INT(POSION(2)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
      N1LOW= INT(POSION(1)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

      N3HI = INT(POSION(3)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
      N2HI = INT(POSION(2)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
      N1HI = INT(POSION(1)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX

!-----------------------------------------------------------------------
! MPI version z ist the fast index
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! conventional version x is fast index
!-----------------------------------------------------------------------
      IND=1

      DO N3=N3LOW,N3HI
      X3=(N3*F3-POSION(3))
      N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

      DO N2=N2LOW,N2HI
      X2=(N2*F2-POSION(2))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      NCOL=GRID%RL%INDEX(N2P,N3P)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-POSION(1))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL) THEN
        N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
        NLI (IND) = N1P+(NCOL-1)*GRID%NGX+1
        IF (NLI(IND) /= 1+N1P+ GRID%NGX*( N2P + GRID%NGY*N3P )) THEN
          WRITE(*,*)'SETYLM ERROR:',N1P,N2P,N3P, NCOL
          STOP
        ENDIF
        ZZ=Z-DISZ
        YY=Y-DISY
        XX=X-DISX
        ! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
        ! is done using the well known formula  | R+d | = | R | + d . R/|R|
        ! this improves the stability of finite differences considerable
        IF (D<1E-4_q) THEN
          DIST(IND)=1E-4_q
        ELSE
          DIST(IND)=MAX(D-(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
        ENDIF

        XS(IND)  =XX/DIST(IND)
        YS(IND)  =YY/DIST(IND)
        ZS(IND)  =ZZ/DIST(IND)

        IND=IND+1
      ENDIF
      ENDDO; ENDDO; ENDDO
!-----------------------------------------------------------------------
!  compare maximum index with INDMAX
!-----------------------------------------------------------------------
      INDMAX=IND-1
      IF (INDMAX>IRMAX) THEN
        WRITE(*,*) &
     &  'internal ERROR: DEPLE:  IRDMAX must be increased to',INT(INDMAX*1.1)
        STOP
      ENDIF
      CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)
!-MM- changes to accomodate spin spirals
!     DEALLOCATE(XS,YS,ZS)
! XS, XY, and XZ, are now returned to calculate the necessary
! phase shifts in SETDIJ
      DO IND=1,INDMAX
         XS(IND)=XS(IND)*DIST(IND)
         YS(IND)=YS(IND)*DIST(IND)
         ZS(IND)=ZS(IND)*DIST(IND)
      ENDDO
!-MM- end of alterations

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE MAXL ******************************
!
! calculate the maximum L quantum number for the augmentation charges
! (required to allocate the tables for e.g. the spherical harmonics)
!
!***********************************************************************

   FUNCTION MAXL_AUG(NTYP,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL_AUG,NTYP
      TYPE (potcar) P(NTYP)
   ! local varibale
      INTEGER I,NT,LTMP,CHANNELS

      MAXL_AUG=0

      DO NT=1,NTYP
         CHANNELS=P(NT)%LMAX
         LTMP=0
         DO I=1,CHANNELS
            LTMP=MAX( P(NT)%LPS(I),LTMP )
         ENDDO
         ! paw requires 2*L for the augmentation charges
         IF ( ASSOCIATED( P(NT)%QPAW) ) THEN
            LTMP=LTMP*2
         ENDIF
         MAXL_AUG=MAX( MAXL_AUG, LTMP)
      END DO

    END FUNCTION MAXL_AUG

!************************ SUBROUTINE MAXL1 *****************************
!
! calculate the maximum L quantum number for (1._q,0._q) particular type
!
!***********************************************************************

   FUNCTION MAXL1(P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL1
      TYPE (potcar) P
   ! local varibale
      INTEGER I,LTMP,CHANNELS

      MAXL1=0

      CHANNELS=P%LMAX
      LTMP=0
      DO I=1,CHANNELS
         LTMP=MAX( P%LPS(I),LTMP )
      ENDDO
      MAXL1=LTMP

    END FUNCTION MAXL1

!************************ SUBROUTINE MAXL ******************************
!
! calculate the maximum L quantum number found in all pseudopotential
! arrays
!
!***********************************************************************

   FUNCTION MAXL(NTYP,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL,NTYP
      TYPE (potcar) P(NTYP)
   ! local varibale
      INTEGER I,NT,LTMP,CHANNELS

      MAXL=0

      DO NT=1,NTYP
         CHANNELS=P(NT)%LMAX
         LTMP=0
         DO I=1,CHANNELS
            LTMP=MAX( P(NT)%LPS(I),LTMP )
         ENDDO
         MAXL=MAX( MAXL, LTMP)
      END DO

    END FUNCTION MAXL

!************************ SUBROUTINE SETDEP ****************************
!
! this subroutine interpolates the augmentation charge on the grid
! around (1._q,0._q) ion using a cubic spline interpolation
! result is returned in DEP
!***********************************************************************

      SUBROUTINE SETDEP(QDEP,PSDMAX,NPSRNL,OMEGA,INDMAX,DIST,DEP)
      USE prec
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)


      DIMENSION QDEP(NPSRNL,5)
      DIMENSION DIST(INDMAX)
      DIMENSION DEP(INDMAX)

      FAKT= OMEGA

      ARGSC=NPSRNL/PSDMAX
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        I  =MIN(INT(DIST(IND)*ARGSC)+1,NPSRNL-1)
        REM=DIST(IND)-QDEP(I,1)
        DEP(IND)=(QDEP(I,2)+REM*(QDEP(I,3)+ &
     &               REM*(QDEP(I,4)+REM*QDEP(I,5))))*FAKT
      ENDDO

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE SETDIJ  ****************************
!
! this subroutine calculates the corrections DION(I,J) corresponding
! to the integral over the augmentation-holes * total local potential
! as input it requires the total potential CVTOT (6.16) (6.12)
! the subroutine also sets up CQIJ(I,J) (i.e. Q(i,j) in thesis)
! DISX,DISY,DISZ are additional displacements of the ions
! small changes required
!
!***********************************************************************

      SUBROUTINE SETDIJ(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,CQIJ, CVTOT_, IRDMAA,IRDMAX, DISX,DISY,DISZ)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      COMPLEX(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),TARGET :: CVTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      LOGICAL  LOVERL,LADDITIONAL
!  work arrays
      REAL(q)   DLM(256)
      REAL(q)   ,ALLOCATABLE ::   DIST(:),DEP(:),POT(:),YLM(:,:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      COMPLEX(q),POINTER :: CVTOT(:),CWORK(:)
!-MM- spin spiral stuff
      REAL(q) QVEC(3),QR
      REAL(q),ALLOCATABLE :: XS(:),YS(:),ZS(:)
!-MM- end of addition
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)
      CDIJ  =0
      IRDMAA=0
      
 overl: IF (LOVERL) THEN
      ! storage convention to (total,magnetization)
      CALL RL_FLIP(CVTOT_(1,1),GRIDC_,WDES%NCDIJ,.FALSE.)

      LYDIM=MAXL_AUG(T_INFO%NTYP,P)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),POT(IRDMAX), &
     &          YLM(IRDMAX,LMYDIM),NLI(IRDMAX))
      IF (LADDITIONAL) THEN
         ALLOCATE(CVTOT(GRIDUS%MPLWV),CWORK(GRIDC_%MPLWV))
      ENDIF
      
!-MM- spin spiral stuff
      ALLOCATE(XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (WDES%LSPIRAL) THEN
      ! Take QSPIRAL from direct to cartesian coordinates
         QVEC(1)=WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3)
         QVEC(2)=WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3)
         QVEC(3)=WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3)
      ENDIF
!-MM- end of addition      
 
 spin:DO ISP=1,WDES%NCDIJ

      IF (LADDITIONAL) THEN
         RINPL=1._q/GRIDC_%NPLWV
         CALL RL_ADD(CVTOT_(1,ISP),RINPL,CWORK,0.0_q,CWORK,GRIDC_)
         CALL FFT3RC(CWORK(1),GRIDC_,-1)

         CVTOT=0
         CALL CPB_GRID(GRIDUS,GRIDC_,C_TO_US,CWORK(1),CVTOT(1))
         CALL FFT3RC(CVTOT(1),GRIDUS,1)
         GRIDC => GRIDUS
      ELSE
         CVTOT => CVTOT_(:,ISP)
         GRIDC => GRIDC_
      ENDIF
      RINPL=1._q/GRIDC%NPLWV

!=======================================================================
! loop over all ions
!=======================================================================

      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      ! for this ion (this type of ion) no depletion charge
      IF (P(NT)%PSDMAX==0) CYCLE
!-----------------------------------------------------------------------
! calculate the spherical harmonics YLM and the distance DIST  between
! grid-points and central atom (DEP and POT are work-arrays)
!-----------------------------------------------------------------------
      LYMAX=MAXL1(P(NT))
      IF ( ASSOCIATED(P(NT)%QPAW) ) THEN
      ! in paw method we truncate the augmentation charge at L=2
         LYMAX=MIN(4,LYMAX*2)
      ENDIF
      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),P(NT)%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISX,DISY,DISZ,DIST(1),NLI(1),XS(1),YS(1),ZS(1))

      IRDMAA=MAX(IRDMAA,INDMAX)

      DO N=1,INDMAX
         POT(N)=CVTOT(NLI(N))
      ENDDO
!=======================================================================
! US-PP
!=======================================================================
  lpaw: IF ( .NOT. ASSOCIATED(P(NT)%QPAW) ) THEN
      LDEP_INDEX=1

    ! loop over all channels (l,epsilon)
      LM =1
      l_loop:  DO L =1,P(NT)%LMAX
      LMP=LM
      lp_loop: DO LP=L,P(NT)%LMAX
      IF (P(NT)%NDEP(L,LP)==0) GOTO 510

      CALL SETDEP(P(NT)%QDEP(1,1,LDEP_INDEX),P(NT)%PSDMAX,NPSRNL, &
     &            LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
      LDEP_INDEX=LDEP_INDEX+ABS(P(NT)%NDEP(L,LP))

    ! quantum numbers l and lp of these two channels
      LL =P(NT)%LPS(L )
      LLP=P(NT)%LPS(LP)

    ! loop over all m mp
      m_loop:  DO M=1,2*LL+1
      MPLOW=1 ; IF (L==LP) MPLOW=M
      mp_loop: DO MP=MPLOW,2*LLP+1

      !   calculate the indices into the array containing the spherical
      !   harmonics
         INDYLM =LL**2   +M
         INDPYL =LLP**2  +MP

         SUM=0

      !   sum over all r
!DIR$ IVDEP
!OCL NOVREC
         DO IND=1,INDMAX
            SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)*YLM(IND,INDPYL)
         ENDDO

         !   add to array CDIJ and make symmetric
         CDIJ(LM+M-1,LMP+MP-1,NI,ISP)=SUM*RINPL
         CDIJ(LMP+MP-1,LM+M-1,NI,ISP)=SUM*RINPL

      ENDDO mp_loop
      ENDDO m_loop

 510  LMP=LMP+2*LLP+1
      ENDDO lp_loop
      LM =LM +2*LL +1
      ENDDO l_loop
   ELSE lpaw
!=======================================================================
! PAW approach
! calculate first the integral int V Y(L,M) Q(L)
! (Q(L) are the L dependent compensation charges in the PAW method)
! around (1._q,0._q) atom
! then transform to  L,L',M,M' use Clebsch-Gordan coefficients
!=======================================================================
!-MM- changes to accommodate spin spirals
      IF (WDES%LSPIRAL .AND. ISP==2) THEN
! V_x -> cos(qr)V_x - sin(qr)V_y
! corresponds to the multiplication of V_12 and V_21 in the spinor 
! representation of the potential by exp(-iqr) and exp(+iqr), respectively
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
            ! calculate q dot r, both given in cartesian coordinates
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM=SUM+DEP(IND)*YLM(IND,INDYLM)* &
              &     (CVTOT_(NLI(IND),2)*COS(QR)-CVTOT_(NLI(IND),3)*SIN(QR))
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
      !   WRITE(0,'("DLM",I2,10F7.4)') L,(DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CDIJ(:,:,NI,ISP), DLM, P(NT))
      ENDIF
      
      IF (WDES%LSPIRAL .AND. ISP==3) THEN
! V_y -> cos(qr)V_y + sin(qr)V_x
! corresponds to the multiplication of V_12 and V_21 in the spinor 
! representation of the potential by exp(-iqr) and exp(+iqr), respectively
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
            ! calculate q dot r, both given in cartesian coordinates
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM=SUM+DEP(IND)*YLM(IND,INDYLM)* &
              &     (CVTOT_(NLI(IND),3)*COS(QR)+CVTOT_(NLI(IND),2)*SIN(QR))
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
      !   WRITE(0,'("DLM",I2,10F7.4)') L,(DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CDIJ(:,:,NI,ISP), DLM, P(NT))
      ENDIF
      
      IF (.NOT.WDES%LSPIRAL .OR. ISP==1 .OR. ISP==4) THEN
! no phase factor for total charge and m_z or if LSPIRAL=.FALSE.
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
            SUMN=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)
               SUMN=SUMN+DEP(IND)*YLM(IND,INDYLM)
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
      !   WRITE(0,'("DLM",I2,10F7.4)') L,(DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CDIJ(:,:,NI,ISP), DLM, P(NT))
      ENDIF
!-MM- end of alterations
      ENDIF lpaw
!=======================================================================
      ENDDO ion
!-----------------------------------------------------------------------
      ENDDO spin

!-MM- changes to accomodate spin spirals
!     DEALLOCATE(DIST,DEP,POT,YLM,NLI)
      DEALLOCATE(DIST,DEP,POT,YLM,NLI,XS,YS,ZS)
!-MM- end of alterations
      IF (LADDITIONAL) DEALLOCATE(CVTOT,CWORK)
! reduce CDIJ
      
      ! back to spinor representation
      CALL RL_FLIP(CVTOT_(1,1),GRIDC_,WDES%NCDIJ,.TRUE.)

      ENDIF overl
!-----------------------------------------------------------------------
! now set up CQIJ and add diagonal part to CDIJ
! find blocks with same quantum number L
! the routine is somewhat complicated
! only terms with same quantum numbers L L' and M M' are non-(0._q,0._q)
!-----------------------------------------------------------------------
      ion2: DO NI=1,T_INFO%NIONS
      NIP=NI_LOCAL(NI, WDES%COMM_INB)
      IF (NIP==0) CYCLE ion2

        LOW=1
        LM =1 
        NT=T_INFO%ITYP(NI)

        DO ISP=1,WDES%NCDIJ
           CQIJ(:,:,NIP,ISP)=0
        ENDDO

        FAKT=WDES%ISPIN

        IF (WDES%LNONCOLLINEAR) THEN 
           ISP_INC=3
           FAKT   =2
        ELSE
           ISP_INC=1
        ENDIF

        block: DO
           LL=P(NT)%LPS(LOW)
           ! search block with same L
           DO LHI=LOW,P(NT)%LMAX
              IF (LL/=P(NT)%LPS(LHI)) EXIT
           ENDDO
           LHI=LHI-1
           MMAX=2*LL+1

           DO ISP=1,WDES%NCDIJ,ISP_INC
           DO L =LOW,LHI
           DO LP=LOW,LHI
           DO M =0,MMAX-1
              CQIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)=P(NT)%QION(L,LP)
           ENDDO
           ENDDO
           ENDDO
           ENDDO

           ISP=1
           DO L =LOW,LHI
           DO LP=LOW,LHI
           DO M =0,MMAX-1
              CDIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)=P(NT)%DION(L,LP)*FAKT &
                   &   +CDIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)
           ENDDO
           ENDDO
           ENDDO

        ! set new LOW value and LM value and goon
           LM=LM+(LHI-LOW+1)*MMAX
           LOW=LHI+1
           IF (LOW > P(NT)%LMAX) EXIT block
        ENDDO block

!       IF (GRIDC%COMM%NODE_ME == GRIDC%COMM%IONODE) THEN
!       DO ISP=1,WDES%NCDIJ
!       WRITE(*,*)'ion',NI,NIP,ISP
!       DO LP=1,P(1)%LMMAX
!         WRITE(*,'(16(F10.3,1X))') (CDIJ(L,LP,NIP,ISP),L=1,MIN(8,P(1)%LMMAX))
!       ENDDO
!       WRITE(*,*)
!       ENDDO
!       ENDIF

      ENDDO ion2

      CALL US_FLIP(WDES, LMDIM, CDIJ, .TRUE., .TRUE.)


      RETURN
      END SUBROUTINE


!************************ SUBROUTINE US_FLIP ***************************
!
! rearranges the storage mode for spin components of array CRHODE:
! given crhode_up and crhode_down on input the quantities
! (crhode_up+crhode_down) and (crhode_up-crhode_down) = total charge
! and magnetization are returned
! also the reverse operation is possible if setting LBACK=.TRUE.
!
!***********************************************************************

      SUBROUTINE US_FLIP(WDES, LMDIM, CRHODE, LOVERL, LBACK)
      USE prec
      USE wave

      IMPLICIT NONE

      TYPE (wavedes)     WDES
      LOGICAL LOVERL, LBACK
      INTEGER LMDIM,LMMAXC,NT,NI,NIS,L,LP,LMBASE
      REAL(q) FAC
      COMPLEX(q) :: CQU,CQD,C01,C10
      COMPLEX(q) :: C11,C00,CX,CY,CZ
      COMPLEX(q) :: CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)


      IF (.NOT.LOVERL) RETURN

      IF (WDES%NCDIJ==2 ) THEN
!=======================================================================
         FAC=1._q
         IF (LBACK) FAC=0.5_q
      
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 100

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CQU=CRHODE(L,LP,NI,1)
                  CQD=CRHODE(L,LP,NI,2)
                  CRHODE(L,LP,NI,1)=FAC*(CQU+CQD)
                  CRHODE(L,LP,NI,2)=FAC*(CQU-CQD)
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 100     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. .NOT. LBACK) THEN
!=======================================================================
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 200

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  C01=CRHODE(L,LP,NI,2)
                  C10=CRHODE(L,LP,NI,3)
                  C11=CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= C00+C11
                  CRHODE(L,LP,NI,2)= C01+C10
                  CRHODE(L,LP,NI,3)=(C01-C10)*(0._q,1._q)
                  CRHODE(L,LP,NI,4)= C00-C11             
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 200     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. LBACK) THEN
!=======================================================================
         FAC=0.5_q
         NIS=1
         LMBASE=0
         typ:  DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 300

            ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  CX =CRHODE(L,LP,NI,2)
                  CY =CRHODE(L,LP,NI,3)
                  CZ =CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= (C00+CZ)*FAC
                  CRHODE(L,LP,NI,2)= (CX-CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,3)= (CX+CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,4)= (C00-CZ)*FAC
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO ion

 300     NIS = NIS+WDES%NITYP(NT)
         ENDDO typ
      ENDIF

      END SUBROUTINE


!************************ SUBROUTINE DEPLE  ****************************
!
! this subroutine calculates  the augmentation charge-density-
! distribution in real space
! as input it requires  CRHODE(LM,LMP,ION,ISP)
! at the end it calculates the total charge-density CHTOT
!
!***********************************************************************

      SUBROUTINE DEPLE(WDES, GRID_SOFT,GRIDC_,GRIDUS,C_TO_US, &
        LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
        LMDIM,CRHODE, CHTOT_,CHDEN, IRDMAX )
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET :: GRID_SOFT,GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER   IRDMAX,ISP      ! allocation required for augmentation
      COMPLEX(q)   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),TARGET  :: CHTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      COMPLEX(q),POINTER :: CHTOT(:)
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      LOGICAL   LOVERL,LADDITIONAL
!  work arrays
      REAL(q)   RHOLM(256)
      REAL(q),ALLOCATABLE ::   DIST(:),DEP(:),SUM(:),YLM(:,:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      LOGICAL L_SYM
!-MM- spin spiral stuff
      REAL(q)   QVEC(3),QR
      REAL(q),ALLOCATABLE :: XS(:),YS(:),ZS(:)
      REAL(q)   RHOLMX(256),RHOLMY(256)
!-MM- end of addition
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)

!-MM- spin spiral stuff
      ALLOCATE(XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (WDES%LSPIRAL) THEN
      ! Take QSPIRAL from direct to cartesian coordinates
         QVEC(1)=WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3)
         QVEC(2)=WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3)
         QVEC(3)=WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3)
      ENDIF
!-MM- end of addition 
!=======================================================================
! if no overlap copy CHDEN to CHTOT and that s it
!=======================================================================
      overl: IF (.NOT.LOVERL) THEN
         DO ISP=1,WDES%NCDIJ
            CALL RC_ADD(CHDEN(1,ISP),1.0_q,CHDEN(1,ISP),0.0_q,CHTOT_(1,ISP),GRID_SOFT)
         ENDDO
      ELSE overl

! find the maximum L for augmentation charge (usually just 2 l)
      LYDIM=MAXL_AUG(T_INFO%NTYP,P)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( &
     &          DIST(IRDMAX),DEP(IRDMAX),SUM(IRDMAX),YLM(IRDMAX,LMYDIM), &
     &          NLI(IRDMAX))

      IF (LADDITIONAL) THEN
         ALLOCATE(CHTOT(GRIDUS%MPLWV))
      ENDIF

!-----------------------------------------------------------------------
! do  symmetrization of the CRHODE
! (in MPI version this is the only position where I can do that
!  without additional communication)
!-----------------------------------------------------------------------
! if PAW is selected, do symmetrization in any case
      L_SYM=.FALSE.
      DO NT=1,T_INFO%NTYP
        IF ( ASSOCIATED(P(NT)%QPAW) ) L_SYM=.TRUE.
      ENDDO
! no symmetry used, well switch it off
      IF (SYMM%ISYM<=0) L_SYM=.FALSE.
! CHDEN is symmetrized and not CHTOT, in that case do symmetrization in any case
      IF (SYMM%ISYM==2) L_SYM=.TRUE.
! now do the symmetrization
      IF (L_SYM) THEN
         IF (WDES%LNONCOLLINEAR) THEN
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP(1), &
              CRHODE(1,1,1,1), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, 1)
! Marsman insert symmetrization here
! Symmetrize the vectors (DX,DY,DZ)
           IF (.NOT.WDES%LSPIRAL) &
          &   CALL AUGSYM_NONCOL_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                CRHODE(1,1,1,2), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), WDES%SAXIS, LATT_CUR%A, LATT_CUR%B)
         ELSE
         DO ISP=1,WDES%NCDIJ
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP(1), &
              CRHODE(1,1,1,ISP), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, ISP)
         ENDDO
        ENDIF
     ENDIF

!=======================================================================
! now the actual work starts
!=======================================================================
     spin:DO ISP=1,WDES%NCDIJ
      IF (LADDITIONAL) THEN
         GRIDC => GRIDUS
      ELSE
         CHTOT => CHTOT_(:,ISP)
         GRIDC => GRIDC_
      ENDIF
!=======================================================================
! loop over all ions
!=======================================================================
      CHTOT=0
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
!-----------------------------------------------------------------------
! for this ion (this type of ion) no depletion charge
!-----------------------------------------------------------------------
      IF (P(NT)%PSDMAX==0) CYCLE
!-----------------------------------------------------------------------
! calulate the spherical harmonics (DEP is Work-arrays)
!-----------------------------------------------------------------------
      DISX=0
      DISY=0
      DISZ=0

      LYMAX=MAXL1(P(NT))
      IF ( ASSOCIATED(P(NT)%QPAW) ) THEN
         LYMAX=MIN(4,LYMAX*2)
      ENDIF
      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),P(NT)%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISX,DISY,DISZ,DIST(1),NLI(1),XS,YS,ZS)

      SUM=0
!=======================================================================
! US-PP
! now loop over pseudopotential indexes L and LP
!=======================================================================
  lpaw: IF ( .NOT. ASSOCIATED(P(NT)%QPAW) ) THEN
      LDEP_INDEX=1

    ! loop over all channels (l,epsilon)
      LM=1
      l_loop:  DO L =1,P(NT)%LMAX
      LMP=LM
      lp_loop: DO LP=L,P(NT)%LMAX
      IF (P(NT)%NDEP(L,LP)==0) GOTO 510

      CALL SETDEP(P(NT)%QDEP(1,1,LDEP_INDEX),P(NT)%PSDMAX,NPSRNL, &
           LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
      LDEP_INDEX=LDEP_INDEX+ABS(P(NT)%NDEP(L,LP))

    ! quantum numbers l and lp of these two channels
      LL =P(NT)%LPS(L )
      LLP=P(NT)%LPS(LP)

    ! loop over all m mp
      m_loop:  DO M=1,2*LL+1
      MPLOW=1
      IF (L==LP) MPLOW=M
      mp_loop: DO MP=MPLOW,2*LLP+1
         FAKT=1
         IF (LMP+MP/=LM+M) FAKT=2

      !   calculate the indexes into the array containing the spherical
      !   harmonics
         INDYLM =LL **2  +M
         INDPYL =LLP**2  +MP

         TFAKT=CRHODE(LM+M-1,LMP+MP-1,NI,ISP)*FAKT

      !   add augmentation charge (augmentation charge is real)
!DIR$ IVDEP
!OCL NOVREC
         DO IND=1,INDMAX
            SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*YLM(IND,INDPYL)*TFAKT
         ENDDO

      ENDDO mp_loop
      ENDDO m_loop
  510 LMP=LMP+2*LLP+1
      ENDDO lp_loop
      LM =LM +2*LL +1
      ENDDO l_loop
   ELSE lpaw
!=======================================================================
! PAW
! transform CRHODE from the basis L,L',M,M'  to
! the basis L,LP,Lmain,Mmain using Clebsch-Gordan coefficients
! then add the compensation charge to grid
!=======================================================================
!-MM- changes to accommodate spin spirals
      IF (WDES%LSPIRAL .AND. ISP==2) THEN
! calculate cell periodic part of the augmentation magnetization density
! by rotating rho_x and rho_y against the spiral
! rho_x -> cos(qr)rho_x + sin(qr)rho_y
      CALL CALC_RHOLM( LYMAX, CRHODE(:,:,NI,2) , RHOLMX, P(NT))
      CALL CALC_RHOLM( LYMAX, CRHODE(:,:,NI,3) , RHOLMY, P(NT))

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*(RHOLMX(INDYLM)*COS(QR)+RHOLMY(INDYLM)*SIN(QR))
            ENDDO

         ENDDO
      ENDDO
      ENDIF
      
      IF (WDES%LSPIRAL .AND. ISP==3) THEN
! calculate cell periodic part of the augmentation magnetization density
! by rotating rho_x and rho_y against the spiral
! rho_y -> cos(qr)rho_y - sin(qr)rho_x
      CALL CALC_RHOLM( LYMAX, CRHODE(:,:,NI,2) , RHOLMX, P(NT))
      CALL CALC_RHOLM( LYMAX, CRHODE(:,:,NI,3) , RHOLMY, P(NT))

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*(-RHOLMX(INDYLM)*SIN(QR)+RHOLMY(INDYLM)*COS(QR))
            ENDDO

         ENDDO
      ENDDO
      ENDIF
      
      IF (.NOT.WDES%LSPIRAL .OR. ISP==1 .OR. ISP==4) THEN
! no phase factor for total charge and m_z or if LSPIRAL=.FALSE.
      CALL CALC_RHOLM( LYMAX, CRHODE(:,:,NI,ISP) , RHOLM, P(NT))

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(P(NT)%QDEP(1,1,L),P(NT)%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*TFAKT
            ENDDO

         ENDDO
      ENDDO
      ENDIF
    ENDIF lpaw
!=======================================================================
! add the calculated augmentation charge to the total charge
!=======================================================================
      SUMN=0
      DO IND=1,INDMAX
        CHTOT(NLI(IND))=CHTOT(NLI(IND))+SUM(IND)
        SUMN=SUMN+SUM(IND)
      ENDDO
!-----------------------------------------------------------------------
      ENDDO ion
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! transform the charge-density to reciprocal space
! and add valenz-charge-density
!-----------------------------------------------------------------------
      CALL FFT_RC_SCALE(CHTOT(1),CHTOT(1),GRIDC)
      IF (LADDITIONAL) CALL CP_GRID(GRIDUS,GRIDC_,C_TO_US,CHTOT(1),CHTOT_(1,ISP))

      CALL ADD_GRID(GRIDC_,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CHTOT_(1,ISP))

      CALL SETUNB_COMPAT(CHTOT_(1,ISP),GRIDC_)
      ENDDO spin

      
!-MM- changes to accommodate spin spirals
! hard set of m_z to (0._q,0._q); this should not be necessary
      IF (WDES%LSPIRAL.AND.WDES%LZEROZ) CHTOT_(:,4)=0
!-MM- end of addition

      IF (LADDITIONAL) DEALLOCATE(CHTOT)
      DEALLOCATE(DIST,DEP,SUM,YLM,NLI,XS,YS,ZS)

      ENDIF overl



      RETURN
      END SUBROUTINE

!************************************************************************
!
! (non explicit) interface for AUGSYM
!
!************************************************************************
      SUBROUTINE AUGSYM_(P,LMDIM,NIONS,NIOND,NTYP,NITYP,MAT,ROTMAP,MAGROT,A,B,ISP)
      USE prec
      USE paw
      USE pseudo
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P(NTYP)
      INTEGER LMDIM,NIONS,NIOND,NTYP,NITYP(NTYP)
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS)
      REAL(q) MAGROT(48,NPCELL)
      INTEGER ROTMAP(NIOND,48,NIOND)
      REAL(q) A(3,3),B(3,3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      CALL AUGSYM(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
            NIOND,NROTK,NPCELL,ROTMAP,MAGROT,ISYMOP,A,B,ISP)

      END SUBROUTINE
      
!************************************************************************
!
! (non explicit) interface for AUGSYM_NONCOL
!
!************************************************************************
      SUBROUTINE AUGSYM_NONCOL_(P,LMDIM,NIONS,NIOND,NTYP,NITYP,MAT,ROTMAP,MAGROT,SAXIS,A,B)
      USE prec
      USE paw
      USE pseudo
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P(NTYP)
      INTEGER LMDIM,NIONS,NIOND,NTYP,NITYP(NTYP)
      COMPLEX(q) MAT(LMDIM,LMDIM,NIONS,3)
      REAL(q) MAGROT(48,NPCELL),SAXIS(3)
      INTEGER ROTMAP(NIOND,48,NIOND)
      REAL(q) A(3,3),B(3,3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      CALL AUGSYM_NONCOL(P,LMDIM,NIONS,NTYP,NITYP,MAT, &
     &      NIOND,NROTK,NPCELL,ROTMAP,MAGROT,SAXIS,ISYMOP,INVMAP,A,B)

      END SUBROUTINE
