#include "symbol.inc"
!***********************************************************************
! RCS:  $Id: fileio.F,v 1.12 2003/06/27 13:22:18 kresse Exp kresse $
!
!  collection of subroutines to performe File IO
!  not yet realy supported heavily
!  but we want to switch all IO form the main program to this modul
!
!***********************************************************************
      MODULE fileio
      USE prec

!*************************SUBROUTINE OUTPOT ****************************
!   write potential     to a specified unit
!   HEADER is currently created in the main program
!***********************************************************************

      INTERFACE
      SUBROUTINE OUTPOT(GRIDC, IU,LLONG,CVTOT)
      USE prec
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRIDC
      COMPLEX(q) CVTOT(GRIDC%RC%NP)
      LOGICAL LLONG
      END SUBROUTINE
      END INTERFACE

      CONTAINS

!*************************SUBROUTINE INWAV  ****************************
!
!   read wavefunctions header from file
!   just parse the old cell shape old cutoff etc.
!   this is a new version that uses only real quantities and allows
!   the WAVECAR file to be exchanged between different machines
!   more easily (IEEE standard)
!
!***********************************************************************

      SUBROUTINE INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IU0)
      USE prec
      USE wave
      USE lattice
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (wavedes)  WDES
      TYPE (latt)     LATT_INI,LATT_CUR
      LOGICAL     LDIFF
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif

      READ(12,REC=2,ERR=100) RKPTSF,RBANDF,ENMAXI, &
     &                       ((LATT_INI%A(I,J),I=1,3),J=1,3)

      IF (IU0 >= 0) WRITE(IU0,*) 'found WAVECAR, reading the header'

      NKPTSF=NINT(RKPTSF)
      NBANDF=NINT(RBANDF)

      CALL LATTIC(LATT_INI)

      IF (ISTART==2 .AND. ENMAXI /= WDES%ENMAX) THEN
        IF (IU0>=0) WRITE(IU0,*) 'ERROR: ENMAX changed please set ISTART to 1'
        STOP
      ENDIF

      IF ((NBANDF/=WDES%NB_TOT .AND. NBANDF*WDES%NRSPINORS/=WDES%NB_TOT)) THEN
!           .OR. NKPTSF/=WDES%NKPTS) THEN
        IF (NBANDF/=WDES%NB_TOT .AND. NBANDF*WDES%NRSPINORS/=WDES%NB_TOT) THEN
           IF (IU0 >= 0) WRITE(IU0,'(2X,A,I4,A,I4)') &
              'nup: number of bands has changed, file:',NBANDF,' present:',WDES%NB_TOT
        ENDIF
        GOTO 100
      ENDIF
      IF (NKPTSF/=WDES%NKPTS) THEN
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I4,A,I4)') &
              'number of k-points has changed, file:',NKPTSF,' present:',WDES%NKPTS
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I4,A,I4)') &
              'trying to continue reading WAVECAR, but it might fail'
      ENDIF
      IF (ISTART ==1) THEN

      LDIFF=.FALSE.
      DO I=1,3
        DO J=1,3
          IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
        ENDDO
      ENDDO
      IF (ENMAXI /= WDES%ENMAX) LDIFF=.TRUE.
      IF (LDIFF) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
      ENDIF

      ENDIF

      RETURN
!
! WAVECAR does not exist or can not be read at all 
!
  100 CONTINUE
      LATT_INI%A=LATT_CUR%A
      CALL LATTIC(LATT_INI)

      IF (ISTART==3) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)"ERROR: can't restart (ISTART=3) with wavefunctions on file"
        STOP
      ELSE
        ISTART=0
      ENDIF
      RETURN

      END SUBROUTINE


!*************************SUBROUTINE INWAV  ****************************
!
!   read wavefunctions from file and repad them according to the new
!   cutoff
!   the repadding is at the moment only supported on the
!   full grid (i.e. column layout is not supported)
!   it is possible to restart a spin polarised calculation
!   from a non spin polarised wavefunction file
!   it also possible to restart a non collinear calculation from
!   a collinear calculation
!
!   this routine is now replaced by the similar routine INWAV_FAST
!
!***********************************************************************

      SUBROUTINE INWAV(WDES,W,GRID,IO,LATT_INI, ISTART)
      USE prec
      USE wave
      USE base
      USE mpimy
      USE mgrid
      USE lattice

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (wavedes)  WDES,WDESI
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID,GRIDI
      TYPE (in_struct) IO
      TYPE (latt)      LATT_INI

      COMPLEX(q) , ALLOCATABLE :: CWORK(:),CW(:),EIG(:)
#ifdef WAVECAR_double
      COMPLEX(q), ALLOCATABLE :: CRD(:)
#else
      COMPLEX(qs), ALLOCATABLE :: CRD(:)
#endif

      INTEGER TIU6, TIU0
      TIU6 = IO%IU6       ! some SGI compilers fail for WRITE(IO%IU6,...)
      TIU0 = IO%IU0

      NODE_ME=0
      IONODE =0
!
! parse the header
!
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      IF (WDES%COMM%NODE_ME == WDES%COMM%IONODE) &
          WRITE(*,*)'internal error: INWAV not supported in MPI version'
      STOP
#endif
        WDESI=WDES
        GRIDI=GRID

        WDESI%NKPTS =1
        WDESI%NKDIM =1
        ALLOCATE(WDESI%VKPT (3,1))

        IF (IO%IU6>=0) THEN
           WRITE(TIU6,*)
           WRITE(TIU6,*)'WAVECAR: Reading wavefunctions'
        ENDIF

        READ(12,REC=1,ERR=200) RDUM,RISPIN
        WDESI%ISPIN=NINT(RISPIN)

        READ(12,REC=2,ERR=200) RKPTSF,RBANDS,WDESI%ENMAX, &
     &                         ((LATT_INI%A(I,J),I=1,3),J=1,3)
        IREC=2
        NBANDF      =NINT(RBANDS)
        NKPTSF      =NINT(RKPTSF)
        WDESI%NBANDS=NINT(RBANDS)
        IF (ISTART==2 .AND. WDESI%ENMAX /= WDES%ENMAX) THEN
          IF (IO%IU0>=0) &
          WRITE(TIU0,*)'ERROR: ENMAX changed please set ISTART to 1'
          STOP
        ENDIF
        CALL LATTIC(LATT_INI)
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
      IF ( WDESI%NBANDS == WDES%NBANDS ) THEN

        ALLOCATE(EIG(WDESI%NBANDS))

        spin:    DO ISP=1,MIN(WDESI%ISPIN,WDES%ISPIN)
        kpoints: DO K=1,WDES%NKPTS

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        IREC=IREC+1
        READ(12,REC=IREC,ERR=200) RNPL,WDESI%VKPT(1:3,1), &
                               (EIG(J),W%FERWE(J,K,ISP),J=1,WDES%NBANDS)
        IF (ISTART==2 .AND. &
         (    ABS(WDESI%VKPT(1,1)- WDES%VKPT(1,K)) > 1E-10_q .OR. &
              ABS(WDESI%VKPT(2,1)- WDES%VKPT(2,K)) > 1E-10_q .OR. &
              ABS(WDESI%VKPT(3,1)- WDES%VKPT(3,K)) > 1E-10_q  )) THEN
          IF (IO%IU0>=0) &
          WRITE(TIU0,*)'ERROR: k-point changed set ISTART to 1'
          STOP
        ENDIF
        NPLREAD=NINT(RNPL)
! calculate indexing scheme (.F. means no allocation for kinetic energy arrays)
        CALL GEN_LAYOUT(GRIDI,WDESI, LATT_INI%B,LATT_INI%B,-1, .FALSE.)
        CALL GEN_INDEX(GRIDI,WDESI, LATT_INI%B,LATT_INI%B,-1,-1, .FALSE.)

        NPL=WDESI%NPLWKP(1)
        NGVECTOR=WDESI%NGVECTOR(1)

        IF (NPL /= NPLREAD) GOTO 200

        MALLOC=MAX(GRIDI%MPLWV,GRID%MPLWV)
        ALLOCATE(CWORK(MALLOC),CW(NPL),CRD(NPL))
!
! read band and repad it if required
!
        band: DO J=1,WDES%NBANDS
          W%CELEN(J,K,ISP)=EIG(J)
          IREC=IREC+1
          READ(12,REC=IREC,ERR=200) (CRD(I),I=1,NPL)
          CW(1:NPL)=CRD(1:NPL)
          CWORK=0
          DO ISPINOR=0,WDES%NRSPINORS-1
          DO I=1,NGVECTOR
             II=I+ISPINOR*NGVECTOR
             CWORK(WDESI%NINDPW(I,1))= CW(II)
          ENDDO
          DO I=1,WDES%NGVECTOR(K)
             II=I+ISPINOR*WDES%NGVECTOR(K)
             W%CPTWFP(II,J,K,ISP)=CWORK(WDES%NINDPW(I,K))
          ENDDO
          ENDDO
        ENDDO band

        DEALLOCATE(CWORK,CW,CRD)
        CALL DEALLOCWDES(WDESI,.FALSE.)
        CALL DEALLOC_GRD(GRIDI)

        ENDDO kpoints
        ENDDO spin
        DEALLOCATE(EIG)

        IF (TIU0 >= 0) WRITE(TIU0,*) 'the WAVECAR file was read sucessfully'

        IF (WDES%ISPIN<=WDESI%ISPIN) RETURN
!
!  spin down is missing
!
        IF (IO%IU0>=0) &
        WRITE(TIU0,*) 'No down-spin wavefunctions found', &
     &             ' --> setting down-spin equal up-spin ...'

        DO K=1,WDES%NKPTS
          W%CELEN(1:WDES%NBANDS,K,2)=W%CELEN(1:WDES%NBANDS,K,1)
          W%FERWE(1:WDES%NBANDS,K,2)=W%FERWE(1:WDES%NBANDS,K,1)
          NPL=WDES%NPLWKP(K)
          W%CPTWFP(1:NPL,1:WDES%NBANDS,K,2)=W%CPTWFP(1:NPL,1:WDES%NBANDS,K,1)
        ENDDO
        RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
      ELSE IF (WDESI%NBANDS*WDES%NRSPINORS == WDES%NBANDS ) THEN

        WDESI%NRSPINORS=1
        ALLOCATE(EIG(WDESI%NBANDS))

        spin2:    DO ISP=1,WDESI%ISPIN
        IF (IO%IU0>=0.AND. ISP==1) &
        WRITE(TIU0,*) 'reading wavefunctions of collinear run, up'
        IF (IO%IU0>=0.AND. ISP==2) &
        WRITE(TIU0,*) 'reading wavefunctions of collinear run, down'
        kpoints2: DO K=1,NKPTSF

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        IREC=IREC+1
        READ(12,REC=IREC,ERR=200) RNPL,WDESI%VKPT(1:3,1), &
               (EIG(J),W%FERWE(J+WDESI%NBANDS*(ISP-1),K,1),J=1,WDESI%NBANDS)
        IF (ISTART==2 .AND. &
         (    ABS(WDESI%VKPT(1,1)- WDES%VKPT(1,K)) > 1E-10_q .OR. &
              ABS(WDESI%VKPT(2,1)- WDES%VKPT(2,K)) > 1E-10_q .OR. &
              ABS(WDESI%VKPT(3,1)- WDES%VKPT(3,K)) > 1E-10_q  )) THEN
          IF (IO%IU0>=0) &
          WRITE(TIU0,*)'ERROR: k-point changed set ISTART to 1'
          STOP
        ENDIF
        NPLREAD=NINT(RNPL)
! calculate indexing scheme (.F. means no allocation for kinetic energy arrays)
        CALL GEN_LAYOUT(GRIDI,WDESI, LATT_INI%B,LATT_INI%B,-1, .FALSE.)
        CALL GEN_INDEX(GRIDI,WDESI, LATT_INI%B,LATT_INI%B,-1,-1, .FALSE.)

        NPL=WDESI%NPLWKP(1)
        NGVECTOR=WDESI%NGVECTOR(1)

        IF (NPL /= NPLREAD) GOTO 200

        MALLOC=MAX(GRIDI%MPLWV,GRID%MPLWV)
        ALLOCATE(CWORK(MALLOC),CW(NPL),CRD(NPL))
!
! read band and repad it if required
! store the up wavefunctions in the upper spinor for 1...WDESI%NBANDS
! and the down wavefunctions in the lower spinor for WDESI%NBANDS+1...2*WDESI%NBANDS
! 
        IF (ISP==1) THEN
        DO J=1,WDESI%NBANDS
          W%CELEN(J,K,1)=EIG(J)

          ! copy immedeatly to second panel
          W%CELEN(J+WDESI%NBANDS,K,1)=EIG(J)
          W%FERWE(J+WDESI%NBANDS,K,1)=W%FERWE(J,K,1)

          IREC=IREC+1
          READ(12,REC=IREC,ERR=200) (CRD(I),I=1,NPL)
          CW(1:NPL)=CRD(1:NPL)
          CWORK=0
          DO I=1,NGVECTOR
             CWORK(WDESI%NINDPW(I,1))= CW(I)
          ENDDO
          ISPINOR=0
          W%CPTWFP(:,J,K,ISP)=0
          DO I=1,WDES%NGVECTOR(K)
             II=I+ISPINOR*WDES%NGVECTOR(K)
             W%CPTWFP(II,J,K,1)=CWORK(WDES%NINDPW(I,K))
          ENDDO

          ISPINOR=1
          W%CPTWFP(:,J+WDESI%NBANDS,K,1)=0
          DO I=1,WDES%NGVECTOR(K)
             II=I+ISPINOR*WDES%NGVECTOR(K)
             W%CPTWFP(II,J+WDESI%NBANDS,K,1)=CWORK(WDES%NINDPW(I,K))
          ENDDO
        ENDDO
        ELSE
        DO J=1,WDESI%NBANDS
          W%CELEN(J+WDESI%NBANDS,K,1)=EIG(J)
          IREC=IREC+1
          READ(12,REC=IREC,ERR=200) (CRD(I),I=1,NPL)
          CW(1:NPL)=CRD(1:NPL)
          CWORK=0
          DO I=1,NGVECTOR
             CWORK(WDESI%NINDPW(I,1))= CW(I)
          ENDDO
          ISPINOR=1
          DO I=1,WDES%NGVECTOR(K)
             II=I+ISPINOR*WDES%NGVECTOR(K)
             W%CPTWFP(II,J+WDESI%NBANDS,K,1)=CWORK(WDES%NINDPW(I,K))
          ENDDO
        ENDDO
        ENDIF

        DEALLOCATE(CWORK,CW,CRD)
        CALL DEALLOCWDES(WDESI,.FALSE.)
        CALL DEALLOC_GRD(GRIDI)

        ENDDO kpoints2
        ENDDO spin2
        DEALLOCATE(EIG)

        IF (TIU0 >= 0) WRITE(TIU0,*) 'the WAVECAR file was read sucessfully'
        RETURN
      ENDIF

!
!  number of bands changed 
!
      IF (IO%IU0>=0) THEN
         WRITE(TIU0,*)'ERROR: while reading wavefunctions, file is incompatible'
         IF  (WDESI%NBANDS /= WDES%NBANDS) &
         WRITE(*,*)'the number of bands has changed (new,old) ',WDESI%NBANDS,WDES%NBANDS
         IF (NKPTSF /= WDES%NKPTS) &
         WRITE(*,*)'the number of k-points has changed (new,old) ',WDES%NKPTS,NKPTSF
      ENDIF
      STOP
!=======================================================================
! can not do anything with WAVECAR 
!=======================================================================
 200  CONTINUE
      IF (IO%IU0>=0) &
      WRITE(TIU0,*)'ERROR: while reading wavefunctions file is possibly corrupt'

      STOP

      END SUBROUTINE


!*************************SUBROUTINE INWAV_FAST ************************
!
!   read wavefunctions from file
!
!   the routine can read WAVECAR files generated for different
!   lattice vectors and cutoff (both in the serial and parallel case)
!
!   it is also possible to restart a spin polarised calculation
!   from a non spin polarised wavefunction file or
!   to restart a non collinear calculation from a collinear calculation
!
!***********************************************************************

      SUBROUTINE INWAV_FAST(WDES, W, GRID, LATT_CUR, LATT_INI, ISTART, IU0)
      USE prec
      USE wave
      USE mgrid
      USE lattice

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR,LATT_INI

! local work arrays
      REAL(q) VKPT(3)
      INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
      COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
      LOGICAL ::  SINGLE_PREC
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)

      NALLOC = MAXVAL(WDES%NPLWKP_TOT)

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif
!
! parse the header
!
        RTAG=0
        READ(12,REC=1,ERR=200) RDUM,RISPIN,RTAG
        IF (RTAG==45200) THEN
           SINGLE_PREC=.TRUE.
        ELSE
           IF (IU0>=0) WRITE(IU0,*) "double precession WAVECAR encountered, converting it"
           SINGLE_PREC=.FALSE.
        ENDIF

        ISPINREAD=NINT(RISPIN)
        READ(12,REC=2,ERR=200) RKPTSF,RBANDF,ENMAXF, &
     &                         ((LATT_INI%A(I,J),I=1,3),J=1,3)

        IREC=2
        NKPTSF=NINT(RKPTSF)
        NBANDF=NINT(RBANDF)

        CALL LATTIC(LATT_INI)
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
      IF ( NBANDF == WDES%NB_TOT ) THEN

        spin:    DO ISP=1, MIN(WDES%ISPIN, ISPINREAD)
        kpoints: DO K=1,WDES%NKPTS

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)
        IREC=IREC+1

        io_begin
        READ(12,REC=IREC,ERR=230) RNPL,VKPT, &
                       (W%CELTOT(J,K,ISP),W%FERTOT(J,K,ISP),J=1,WDES%NB_TOT)

        NPLREAD=NINT(RNPL)
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
        NPL=WDES%NPLWKP_TOT(K)

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(MALLOC),CW2(MALLOC),CRD(MALLOC),IND(MALLOC),INDI(MALLOC))

        io_begin
        CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 220

        band: DO J=1,WDES%NB_TOT
          IREC=IREC+1
          io_begin
            IF (SINGLE_PREC) THEN
               READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
               CW2(1:NPLREAD)=CRD(1:NPLREAD)
            ELSE
               READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
            ENDIF

            CW1=0
            DO IS=1,WDES%NRSPINORS
               CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                  CW1((IS-1)*NPL/WDES%NRSPINORS+1), CW2((IS-1)*NPLREAD/WDES%NRSPINORS+1))
            ENDDO

          io_end
          CALL DIS_PW_BAND(WDES1, J, CW1, W%CW(1,1,K,ISP))
        ENDDO band

        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(WDES%NB_TOT+1)
        ENDIF

        ! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALLMPI( M_bcast_z(WDES%COMM, W%CELTOT(1,1,ISP),NCOMM ))
        CALLMPI( M_bcast_d(WDES%COMM, W%FERTOT(1,1,ISP),NCOMM ))

        ENDDO spin

        IF (ISPINREAD > WDES%ISPIN .AND. IU0>=0 ) THEN
           WRITE(IU0,*) 'down-spin wavefunctions not read'
        ENDIF
        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read sucessfully'
        IF (WDES%ISPIN<=ISPINREAD) RETURN
!
!  spin down is missing
!
        IF (IU0>=0) &
        WRITE(IU0,*) 'No down-spin wavefunctions found', &
     &             ' --> setting down-spin equal up-spin ...'

        DO K=1,WDES%NKPTS
          W%CELTOT(1:WDES%NB_TOT,K,2)=W%CELTOT(1:WDES%NB_TOT,K,1)
          W%FERTOT(1:WDES%NB_TOT,K,2)=W%FERTOT(1:WDES%NB_TOT,K,1)
          NPL=WDES%NPLWKP(K)
          W%CW(1:NPL,1:WDES%NBANDS,K,2)=W%CW(1:NPL,1:WDES%NBANDS,K,1)
        ENDDO

        RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
      ELSE IF (NBANDF*WDES%NRSPINORS == WDES%NB_TOT ) THEN
        W%CW=0

        spin2:    DO ISP=1,ISPINREAD
        IF (IU0>=0.AND. ISP==1) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
        IF (IU0>=0.AND. ISP==2) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
        kpoints2: DO K=1,NKPTSF

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)
        IREC=IREC+1

        io_begin
        READ(12,REC=IREC,ERR=230) RNPL, VKPT(1:3), &
                       (W%CELTOT(J+NBANDF*(ISP-1),K,1),W%FERTOT(J+NBANDF*(ISP-1),K,1),J=1,NBANDF)

        NPLREAD=NINT(RNPL)
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
        NPL=WDES%NPLWKP_TOT(K)/2

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),CRD(2*MALLOC),IND(MALLOC),INDI(MALLOC))

        io_begin
        CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 220

        IF (ISP==1) THEN
           DO J=1,NBANDF
              IREC=IREC+1
              io_begin
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF
                   
                      
                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
              io_end
              CALL DIS_PW_BAND(WDES1, J, CW1, W%CW(1,1,K,1))

              ! copy immedeatly to second panel
              W%CELTOT(J+NBANDF,K,1)=W%CELTOT(J,K,1)
              W%FERTOT(J+NBANDF,K,1)=W%FERTOT(J,K,1)

              CW1(NPL+1:2*NPL)=CW1(1:NPL)
              CW1(1:NPL)=0

              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CW(1,1,K,1))
           ENDDO
        ELSE
           DO J=1,NBANDF
              IREC=IREC+1
              io_begin
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF
                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
              io_end
              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CW(1,1,K,1))
           ENDDO

        ENDIF

        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints2

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(NBANDF+1)
        ENDIF

        ENDDO spin2

        ! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALLMPI( M_bcast_z(WDES%COMM, W%CELTOT(1,1,1),NCOMM ))
        CALLMPI( M_bcast_d(WDES%COMM, W%FERTOT(1,1,1),NCOMM ))

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read sucessfully'
        RETURN
      ENDIF

 210  CONTINUE
      IF (IU0>=0) THEN
         WRITE(IU0,*)'ERROR: while reading WAVECAR, file is incompatible'
         IF (WDES%ENMAX /= ENMAXF) &
            WRITE(*,*)'the energy cutoff has changed (new,old) ',WDES%ENMAX,ENMAXF
         IF  (NBANDF /= WDES%NB_TOT) &
            WRITE(*,*)'the number of bands has changed (new,old) ',WDES%NB_TOT,NBANDF
         IF (NKPTSF /= WDES%NKPTS) &
            WRITE(*,*)'the number of k-points has changed (new,old) ',WDES%NKPTS,NKPTSF
      ENDIF
      STOP
!=======================================================================
! can not do anything with WAVECAR 
! hard stop, pull all breaks
!=======================================================================
  200 CONTINUE

      IF (IU0>=0) &
      WRITE(IU0,*)'ERROR: while reading WAVECAR, header is corrupt'
      STOP

  220 CONTINUE

      IF (IU0>=0) &
      WRITE(IU0,*)'ERROR: while reading WAVECAR, plane wave coefficients changed', &
          NPL,NPLREAD
      STOP

  230 CONTINUE

      IF (IU0>=0) &
      WRITE(IU0,*)'ERROR: while reading eigenvalues from WAVECAR',K,ISP
      STOP

  240 CONTINUE

      IF (IU0>=0) &
      WRITE(IU0,*)'ERROR: while reading plane wave coef. from WAVECAR',K,ISP,J
      STOP
      
      END SUBROUTINE


!*************************SUBROUTINE OUTWAV    ************************
!
!   write wavefunctions to file layout must be correct
!
!***********************************************************************

      SUBROUTINE OUTWAV(IRECLW,WDES,W,LATT_INI,IU0)
      USE prec
      USE wave
      USE lattice
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER IRECLW    ! record length
      TYPE (latt)     LATT_INI
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
! local work arrays
      COMPLEX(q),ALLOCATABLE :: CW(:),EIG(:)
#ifdef WAVECAR_double
      COMPLEX(q), ALLOCATABLE :: CRD(:)
#else
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
#endif

      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
      ALLOCATE(CW(NPL_TOT),CRD(NPL_TOT),EIG(WDES%NB_TOT))

      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif
!
! header
!
      IF (IU0>=0) WRITE(IU0,*)'writing wavefunctions'

      RISPIN=WDES%ISPIN
      RDUM  =IRECLW

#ifdef WAVECAR_double
      RTAG  =45210
#else
      RTAG  =45200
#endif
      do_io WRITE(12,REC=1) RDUM,RISPIN,RTAG

      IREC=2
      io_begin
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
      RNKPTS =WDES%NKPTS
      RNB_TOT=WDES%NB_TOT
      WRITE(12,REC=2) RNKPTS,RNB_TOT,WDES%ENMAX,((LATT_INI%A(I,J),I=1,3),J=1,3)
      io_end
!
! loop over spin, kpoints, bands
!
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO K=1,WDES%NKPTS

       CALL SETWDES(WDES,WDES1,K)
       NPL=WDES%NPLWKP_TOT(K)
       RNPL=NPL
       IREC=IREC+1
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
       io_begin
! write eigenvalues in real format
         DO J=1,WDES%NB_TOT
           EIG(J)=REAL(W%CELTOT(J,K,ISP),KIND=q)
         ENDDO
         WRITE(12,REC=IREC) RNPL,WDES%VKPT(1,K),WDES%VKPT(2,K), &
                       WDES%VKPT(3,K),(EIG(J),W%FERTOT(J,K,ISP),J=1,WDES%NB_TOT)
       io_end

       DO J=1,WDES%NB_TOT
         CALL MRG_PW_BAND(WDES1, J, CW, W%CW(1,1,K,ISP))
         IREC=IREC+1
         do_io CRD(1:NPL)=CW(1:NPL)
         do_io WRITE(12,REC=IREC) (CRD(I),I=1,NPL)
       ENDDO
      ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CW,CRD,EIG)

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OUTCHG ****************************
!
!   write chargedensity to a specified unit
!   HEADER is currently created in the main program
!
!***********************************************************************

      SUBROUTINE OUTCHG(GRIDC, IU, LLONG,CHTOT)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRIDC

      COMPLEX(q)  CHTOT(GRIDC%RC%NP)
      LOGICAL LLONG        ! long or short format
! work arrays
      RGRID,  ALLOCATABLE ::  CWORK(:)
      REAL(q),ALLOCATABLE ::  WORK (:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      CHARACTER*40 FORM
      INTEGER ISTAT

      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(CWORK(DIMREAL(GRIDC%MPLWV)),WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the charge immediate exit

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        WRITE(*,*)'internal ERROR: OUTCHG NPLWV is not compatibel', &
     &   ' with  NGX,NGY,NGZ'
        WRITE(*,*)'   ',GRIDC%NPLWV,GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ
        STOP
      ENDIF
      ! transfer to CWORK and FFT
      CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
      CALL FFT3RC(CWORK,GRIDC,1)

      IF (LLONG) THEN
        FORM='(1(1X,E17.11))'
        NWRITE=5
      ELSE
        FORM='(1(1X,G11.5))'
        NWRITE=10
      ENDIF

      do_io WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CWORK, NZ)
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
         io_end
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(CWORK,WORK)

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE  INCHG ****************************
!   read chargedensity form a specified unit
!   HEADER must have been read previously
!***********************************************************************

      SUBROUTINE INCHG(IU,CHTOT,GRID,IERR)
      USE prec
      USE charge
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q)   CHTOT(GRID%RC%NP)
! work array
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      INTEGER ISTAT

      NALLOC=GRID%NGX*GRID%NGY
      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) THEN
         IERR=2
         RETURN                 ! can not read potential, immediate exit
      ENDIF

      IERR=0
      IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
        WRITE(*,*)'internal ERROR: INCHG GRID%NPLWV is not compatible', &
     &   ' with  GRID%NGX,GRID%NGY,NGZC'
        WRITE(*,*)'   ',GRID%NPLWV,GRID%NGX,GRID%NGY,GRID%NGZ
        STOP
      ENDIF

      IERR=0
      NGXFIL=0; NGYFIL=0; NGZFIL=0

#ifdef MPI
      IF (NODE_ME==IONODE) THEN
#endif
        READ(IU,*,ERR=120,END=120) NGXFIL,NGYFIL,NGZFIL
 120    CONTINUE
        IF (NGXFIL==0 .AND. NGYFIL==0 .AND. NGZFIL==0) THEN
           IERR=2
        ELSE IF (NGXFIL/=GRID%NGX .OR. NGYFIL/=GRID%NGY .OR. NGZFIL/=GRID%NGZ ) THEN
          IERR=1
        ENDIF
#ifdef MPI
      ENDIF
#endif

      CALLMPI( M_bcast_i( GRID%COMM, IERR, 1))
      IF (IERR /=0 ) GOTO 200

      NWRITE=5
      NWRITTEN=0
      DO NZ=1,GRID%NGZ
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            ! non advancing read is not the same on all machine
            ! this is the only version that works on all tested platforms
            IF ( MOD(NWRITTEN,NWRITE)==0 .OR. (N==NALLOC .AND. NZ==GRID%NGZ)) THEN
               READ(IU,'(1(1X,E17.11))',ERR=100,END=100) WORK(N)
            ELSE
               READ(IU,'(1(1X,E17.11))',ADVANCE='NO',ERR=100,END=100) WORK(N)
            ENDIF
         ENDDO
         io_end
         CALL DIS_GRID_RL_PLANE(GRID, WORK, CHTOT, .TRUE., NZ )
      ENDDO
 
      CALL FFT_RC_SCALE(CHTOT,CHTOT,GRID)
      CALL SETUNB_COMPAT(CHTOT,GRID)

  200 CONTINUE
      DEALLOCATE(WORK)
      RETURN

  100 CONTINUE
      IERR=2
      DEALLOCATE(WORK)
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READNI ****************************
!   reads total number of ions from a line consisting of numbers
!   for several species
!***********************************************************************

      SUBROUTINE READNI(ZPARSE,NIONSF)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER*80 ZPARSE
      CHARACTER*80 ZITEM
      CHARACTER*80 ZWORK
      CHARACTER*15 ZFORM
      CHARACTER*1  ZCHAR

!  get number of data items
      NDATA=NITEMS(ZPARSE,ZWORK,.TRUE.,'I')

!  parse and read list:
      NIONSF=0

      DO 100 IDATA=1,NDATA
         CALL SUBWRD(ZPARSE,ZITEM,IDATA,1)
         CALL CHKINT(ZITEM,ZWORK,ZCHAR,ZFORM)
         IF (ZCHAR=='Y') THEN
            ZWORK='('//ZFORM//')'
            CALL STRIP(ZWORK,LFORM,'A')
            READ(ZITEM,ZWORK(1:LFORM)) NI
         ELSE
            WRITE(*,*) 'Fatal error in READNI: Invalid data found ...'
            STOP
         ENDIF

         NIONSF=NIONSF+NI

  100 CONTINUE

      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READCH ****************************
!
!  This routine reads in a chargedensity from file IU
!  when the chargedensity of the file is incompatible
!  i.e. NGX,Y,Z differ from that used in the programm a warning is
!  reported and ICHARG is set to 0
!  if the magnetisation charge density could not be read in 
!  ICHARG is set to -1
!  on entry to the routine:
!     CHTOT  may be uninitialised
!     RHOLM  must be already correctly initialised according to MAGMOM
!  on exit:
!     CHTOT  as read from file
!     RHOLM  as read from file
!  the routine is complicated by several special cases:
!  o if the magnetisation charge density could not be read in 
!     ICHARG is set to -1, CHTOT(:,2:4)=0
!  o non collinear case: if only mz could be read
!     the x and y components are set to 0 in CHTOT and RHOLM
!
!***********************************************************************

      SUBROUTINE READCH(GRIDC, LOVERL, T_INFO, CHTOT, RHOLM, ICHARG, ISPIN, &
           LATT_CUR, P, CSTRF, IU, IU0)
      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE paw

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      INTEGER ISPIN, IU, IU0, ICHARG
      REAL(q)         :: RHOLM(:,:)
! local work arrays
      TYPE (type_info)   T_INFO_OLD
      REAL(q)         :: TMP(T_INFO%NIONS)
      COMPLEX(q) CHTOT(GRIDC%MPLWV,ISPIN)
      CHARACTER*40 TITEL
      CHARACTER*80 ZPARSE
      COMPLEX(q), ALLOCATABLE :: CD(:),CSTRF_OLD(:,:)
      INTEGER I,NIONSF,IERR,N
      LOGICAL LOVERL

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

!=======================================================================
! read in old ionic-positions and header
!=======================================================================
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IU,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IU,*,ERR=1000,END=1000)
      ENDDO
      READ(IU,'(A80)') ZPARSE
      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: number of atoms are different on CHGCAR file'
        ICHARG=0
        RETURN
      ENDIF
      READ(IU,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IU,*) T_INFO_OLD%POSION(:,I)
      ENDDO
!=======================================================================
! read in the charge-density
!=======================================================================
      CALL INCHG(IU,CHTOT(1,1),GRIDC,IERR)
      IF (IERR==1) THEN
        IF (IU0>=0) &
        WRITE(*,*) 'WARNING: dimensions on CHGCAR file are different'
        ICHARG=0
        RETURN
      ENDIF
      IF (IERR==2) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: chargedensity file is incomplete'
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! now subtract charge density according to old position and add
! that one according to new positions
!=======================================================================
!     GOTO 2000
      ALLOCATE(CD(GRIDC%MPLWV),CSTRF_OLD(GRIDC%MPLWV,T_INFO%NTYP))

!---- subtract atomic charge-denisty for old positions
      CALL STUFAK(GRIDC,T_INFO_OLD,CSTRF_OLD)
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO_OLD,LATT_CUR%B,P,CSTRF_OLD,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT(N,1)= CHTOT(N,1)- CD(N)
      ENDDO

!---- add atomic charge-denisty for new positions
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT (N,1)= CHTOT(N,1)+ CD(N)
      ENDDO

      DEALLOCATE(CD,CSTRF_OLD,T_INFO_OLD%POSION)
 2000 CONTINUE

      ! just in case initialize the magnetization density to 0
      DO ISP=2,ISPIN
        CALL RC_ADD(CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),GRIDC)
      ENDDO

      CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,1), GRIDC%COMM, IU, IERR )
      IF (IERR/=0) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: PAW occupancies are missing on CHGCAR '
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! set/read spin density
!=======================================================================
      IF (IU0>=0) &
        WRITE(IU0,*) 'charge-density read from file: ',TITEL

      DO ISP=2,ISPIN
! read in the spin-density
        IERR=0
        ! in vasp.4.4 the ATOMOM array was read from the CHGCAR file
        ! reading the magnetic moments from the file really makes no sense
        ! and spoiles non collinear calculations
        do_io READ(IU,*,IOSTAT=IERR) (TMP(I),I=1,T_INFO%NIONS)
        CALLMPI( M_bcast_i( GRIDC%COMM, IERR, 1))

        IF (IERR==0) CALL INCHG(IU,CHTOT(1,ISP),GRIDC,IERR)
        IF (IERR==0) CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,ISP), GRIDC%COMM, IU, IERR )

        ! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
           ! we could not read any entry in the magnetisation density return ICHARG=-1
           ! (magnetisation according to overlapping atoms)
              ICHARG=-1
              RETURN
           ELSE
           ! non collinear case we could read at least m_z(r)
           ! copy that to the right place and initialise everything else to 0
              CALL RC_ADD(CHTOT(1,2),1.0_q,CHTOT(1,2),0.0_q,CHTOT(1,4),GRIDC) ! copy
              CALL RC_ADD(CHTOT(1,2),0.0_q,CHTOT(1,2),0.0_q,CHTOT(1,2),GRIDC) ! =0
              CALL RC_ADD(CHTOT(1,3),0.0_q,CHTOT(1,3),0.0_q,CHTOT(1,3),GRIDC) ! =0
              
              RHOLM(:,4)=RHOLM(:,2)
              RHOLM(:,2)=0
              RHOLM(:,3)=0

              RETURN
           ENDIF
        ENDIF

        IF (IU0>=0) WRITE(IU0,'(A,I2)') ' magnetization density read from file',ISP-1

      ENDDO
      RETURN

 40   CONTINUE
      ICHARG=-1
      RETURN

 1000 CONTINUE
      IF (IU0>=0) &
       WRITE(IU0,*) 'WARNING: chargedensity file is incomplete'
      ICHARG=0
      RETURN
      END SUBROUTINE
      END MODULE

!*************************SUBROUTINE OUTPOT ****************************
!   write potential     to a specified unit
!   HEADER is currently created in the main program
!***********************************************************************

      SUBROUTINE OUTPOT(GRIDC, IU,LLONG,CVTOT)
      USE prec
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRIDC
      RGRID CVTOT(GRIDC%RL%NP)
      LOGICAL LLONG
      CHARACTER*40 FORM
! local variables
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER ISTAT

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the potential immediate exit

      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        WRITE(*,*)'internal ERROR: OUTPOT GRIDC%NPLWV is not compatibel', &
     &   ' with  GRIDC%NGX,GRIDC%NGY,NGZC'
        WRITE(*,*)'   ',GRIDC%NPLWV,GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ
        STOP
      ENDIF

      IF (LLONG) THEN
        FORM='(1(1X,E17.11))'
        NWRITE=5
      ELSE
        FORM='(1G11.5)'
        NWRITE=10
      ENDIF

      do_io WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CVTOT, NZ)
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
         io_end
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(WORK)

      RETURN
      END SUBROUTINE



!***********************************************************************
!  write out initial header for PCDAT
!***********************************************************************

      SUBROUTINE PCDAT_HEAD(IU, T_INFO, LATT_CUR, DYN, PACO, SZNAM1)
      USE prec
      USE lattice
      USE poscar
      USE base
      IMPLICIT NONE

      INTEGER IU
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      TYPE (dynamics)  :: DYN
      TYPE (paco_struct)  PACO
      CHARACTER*40 SZNAM1
! local variables
      INTEGER I
      REAL(q) AOMEGA
      AOMEGA=LATT_CUR%OMEGA/T_INFO%NIONS

      WRITE(IU,'(4I4,2E15.7)')1,T_INFO%NIONS,1,0,AOMEGA,DYN%TEMP
      WRITE(IU,*) ' CAR '
      WRITE(IU,*) SZNAM1
      WRITE(IU,'(3I4)') 0,0,0
      WRITE(IU,'(2I4)') 1,DYN%NBLOCK
      WRITE(IU,'(3I4)') PACO%NPACO,PACO%NPACO,PACO%NPACO
      WRITE(IU,'(1I4)') PACO%NPACO
      WRITE(IU,'(1E15.7)') 1E-10
      WRITE(IU,'(1E15.7)') PACO%APACO*1E-10/PACO%NPACO
      WRITE(IU,'(1I4)') DYN%NSW/DYN%NBLOCK/DYN%KBLOCK
      WRITE(IU,'(4E15.7)') DYN%POTIM*1E-15,((LATT_CUR%ANORM(I)*1E-10),I=1,3)

      RETURN
      END SUBROUTINE
