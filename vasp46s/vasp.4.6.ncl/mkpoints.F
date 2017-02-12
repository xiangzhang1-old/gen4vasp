#include "symbol.inc"
!************************************************************************
! RCS:  $Id: mkpoints.F,v 1.3 2003/06/27 13:22:20 kresse Exp kresse $
!
!  this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE MKPOINTS
      USE prec
      INCLUDE "mkpoints.inc"
      CONTAINS
!************************************************************************
!  Read UNIT=14: KPOINTS
!  number of k-points and k-points in reciprocal lattice
!  than call the symmetry package to generate the k-points for
!  us
!
!************************************************************************

      SUBROUTINE RD_KPOINTS(KPOINTS,LATT_CUR, NKDIM,IKPTD,NTETD,LINVERSION,IU6,IU0)
      USE prec
      USE lattice

      IMPLICIT NONE

      TYPE (kpoints_struct) KPOINTS
      TYPE (latt)        LATT_CUR

      INTEGER   NKDIM   ! maximal number of k-points
      INTEGER   IKPTD   ! maximal number of division
      INTEGER   NTETD   ! maximal number of tetrahedrons
      INTEGER   IU0,IU6 ! units for output
      CHARACTER*1   CSEL,CLINE
      REAL(q)    BK(3,3),SHIFT(3),SUPL_SHIFT(3),BK_REC(3,3)
! required for reallocation
      REAL(q),POINTER   :: VKPT(:,:),WTKPT(:),VKPT2(:,:)
      INTEGER,POINTER:: IDTET(:,:)
      INTEGER :: IERR,INDEX,NINTER,N

! local variables
      INTEGER KTH,NKPX,NKPY,NKPZ,NKP,ITET,I,NT,NK
      REAL(q) RKLEN,WSUM
      LOGICAL LINVERSION
! warnings from tutor
      INTEGER,PARAMETER :: NTUTOR=1
      REAL(q)     RTUT(NTUTOR),RDUM
      INTEGER  ITUT(NTUTOR),IDUM
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM

      OPEN(UNIT=14,FILE='KPOINTS',STATUS='OLD')

      KPOINTS%NKDIM=NKDIM
      KPOINTS%NTET=0

      ALLOCATE(VKPT(3,NKDIM),WTKPT(NKDIM),IDTET(0:4,NTETD))
      ALLOCATE(KPOINTS%IKPT(IKPTD,IKPTD,IKPTD))

      IF (IU6>=0) WRITE(IU6,*)
!-----K-points
      ITUT(1)=1
      READ(14,'(A40)',ERR=70111,END=70111) KPOINTS%SZNAMK
      IF (IU6>=0) WRITE(IU6,*)'KPOINTS: ',KPOINTS%SZNAMK

      ITUT(1)=ITUT(1)+1
      READ(14,*,ERR=70111,END=70111) KPOINTS%NKPTS
      IF (KPOINTS%NKPTS>KPOINTS%NKDIM) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'ERROR: MAIN: increase NKDIM'
        STOP
      ENDIF
      ITUT(1)=ITUT(1)+1
      READ(14,'(A1)',ERR=70111,END=70111) CSEL
!
! if CSEL is starting with l for line, k-points are interpolated
! between the points read from  KPOINTS
!
      IF (CSEL=='L'.OR.CSEL=='l') THEN
         CLINE='L'
         ITUT(1)=ITUT(1)+1
         READ(14,'(A1)',ERR=70111,END=70111) CSEL
         KPOINTS%NKPTS=MAX( KPOINTS%NKPTS,2)
        IF (IU6 >=0)  &
     &      WRITE(IU6,*)' interpolating k-points between supplied coordinates'

      ELSE
         CLINE=" "
      ENDIF
      
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        IF (IU6 >=0 .AND. KPOINTS%NKPTS>0)  &
     &      WRITE(IU6,*)' k-points in cartesian coordinates'
      ELSE
        IF (IU6 >= 0 .AND. KPOINTS%NKPTS>0)  &
     &     WRITE(IU6,*)' k-points in reciprocal lattice'
      ENDIF

!=======================================================================
! read in a set of k-points and interpolate NKPTS between each
!=======================================================================

      IF (KPOINTS%NKPTS>0) THEN

      kr: IF (CLINE=='L') THEN
         IF (KPOINTS%LTET) THEN
           CALL VTUTOR('E','LINTET',RTUT,1, &
     &          ITUT,3,CDUM,1,LDUM,1,IU6,1)
           CALL VTUTOR('E','LINTET',RTUT,1, &
     &          ITUT,3,CDUM,1,LDUM,1,IU0,1)
           STOP
         ENDIF
         ALLOCATE(VKPT2(3,NKDIM))

         NINTER=KPOINTS%NKPTS
         NKP=0  ! counter for the number of k-points already read in
         DO 
            NKP=NKP+1
            IF (NKP>KPOINTS%NKDIM) THEN
               IF (IU0>=0) &
                    WRITE(IU0,*)'ERROR: MAIN: increase NKDIM'
               STOP
            ENDIF
            ITUT(1)=ITUT(1)+1
            READ(14,*,IOSTAT=IERR) &
     &           VKPT2(1,NKP),VKPT2(2,NKP),VKPT2(3,NKP)
            IF (IERR/=0) EXIT
         ENDDO

         KPOINTS%NKPTS=NKP-1
         IF (CSEL=='K') THEN
            VKPT2(:,1:KPOINTS%NKPTS)=  &
     &           VKPT2(:,1:KPOINTS%NKPTS)/LATT_CUR%SCALE

            CALL KARDIR(KPOINTS%NKPTS,VKPT2,LATT_CUR%A)
         ENDIF

         INDEX=0
         ! make NKPTS even
         KPOINTS%NKPTS=(KPOINTS%NKPTS/2)*2
         DO NKP=1,KPOINTS%NKPTS-1,2
            SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
            DO N=0,NINTER-1
               INDEX=INDEX+1
               IF (INDEX>KPOINTS%NKDIM) THEN
                  IF (IU0>=0) WRITE(IU0,*)'ERROR: MAIN: increase NKDIM'
                  STOP
               ENDIF
               VKPT(:,INDEX)=VKPT2(:,NKP)+SHIFT*N
               WTKPT(INDEX)=1._q/(KPOINTS%NKPTS/2*NINTER)
            ENDDO
         ENDDO
         KPOINTS%NKPTS=(KPOINTS%NKPTS/2)*NINTER

         ! k-point lines
         CALL XML_KPOINTS_3(KPOINTS%NKPTS, VKPT, WTKPT, NINTER)

      ELSE kr
!=======================================================================
! Read in a given set of arbitrary k-points:
!=======================================================================
         WSUM=0
         DO NKP=1,KPOINTS%NKPTS
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) &
     &           VKPT(1,NKP),VKPT(2,NKP),VKPT(3,NKP), &
     &           WTKPT(NKP)
            WSUM=WSUM+WTKPT(NKP)
         ENDDO

         IF (WSUM==0) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'ERROR: sum of weights is zero'
            STOP
         ENDIF

         IF (CSEL=='K') THEN

            VKPT(:,1:KPOINTS%NKPTS)=  &
     &           VKPT(:,1:KPOINTS%NKPTS)/LATT_CUR%SCALE

            CALL KARDIR(KPOINTS%NKPTS,VKPT,LATT_CUR%A)
         ENDIF

         WTKPT(1:KPOINTS%NKPTS)=WTKPT(1:KPOINTS%NKPTS)/WSUM

         IF (KPOINTS%LTET) THEN
! Read in tetrahedra if you want to use tetrahedron method:
            ITUT(1)=ITUT(1)+1
            READ(14,'(A)',ERR=70119,END=70119) CSEL
            IF (CSEL/='T' .AND. CSEL/='t') GOTO 70119
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) KPOINTS%NTET,KPOINTS%VOLWGT
            DO ITET=1,KPOINTS%NTET
               ITUT(1)=ITUT(1)+1
               READ(14,*,ERR=70119,END=70119) (IDTET(KTH,ITET),KTH=0,4)
            ENDDO
         ENDIF
         CALL XML_KPOINTS_1(KPOINTS%NKPTS, VKPT, WTKPT,&
              KPOINTS%NTET, IDTET, KPOINTS%VOLWGT)

       ENDIF kr

      ELSE
!=======================================================================
! Automatic generation of a mesh if KPOINTS%NKPTS<=0:
!=======================================================================
         IF (IU6>=0 ) WRITE(IU6,'(/A)') 'Automatic generation of k-mesh.'
         SHIFT(1)=0._q
         SHIFT(2)=0._q
         SHIFT(3)=0._q
         SUPL_SHIFT=SHIFT
! k-lattice basis vectors in cartesian or reciprocal coordinates?
         IF ((CSEL/='M').AND.(CSEL/='m').AND. &
     &       (CSEL/='G').AND.(CSEL/='g').AND. &
     &       (CSEL/='A').AND.(CSEL/='a')) THEN
! Here give a basis and probably some shift (warning this shift is
! always with respect to the point (0,0,0) ... !
            IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &          CSEL=='C'.OR.CSEL=='c') THEN
               CSEL='K'
               IF (IU6>=0 ) WRITE(IU6,*)' k-lattice basis in cartesian coordinates'
            ELSE
               IF (IU6>=0 )WRITE(IU6,*)' k-lattice basis in reciprocal lattice'
            ENDIF
! Read in the basis vectors for the k-lattice (unscaled!!):
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) BK(1,1),BK(2,1),BK(3,1)
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) BK(1,2),BK(2,2),BK(3,2)
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) BK(1,3),BK(2,3),BK(3,3)
! Correct scaling with LATT_CUR%SCALE ('lattice constant'):
            IF (CSEL=='K') BK=BK/LATT_CUR%SCALE
! Routine IBZKPT needs cartesian coordinates:
            IF (CSEL/='K') THEN
               CALL DIRKAR(3,BK,LATT_CUR%B)
            ENDIF
! Read in the shift of the k-mesh: these values must be given in
! k-lattice basis coordinates (usually 0 or 1/2 ...):
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70112,END=70112) SHIFT(1),SHIFT(2),SHIFT(3)
            SUPL_SHIFT=SHIFT
70112       CONTINUE
         ELSE IF ((CSEL=='A').OR.(CSEL=='a')) THEN
            ITUT(1)=ITUT(1)+1
            READ(14,*) RKLEN
            NKPX =MAX(1._q,RKLEN*LATT_CUR%BNORM(1)+0.5_q)
            NKPY =MAX(1._q,RKLEN*LATT_CUR%BNORM(2)+0.5_q)
            NKPZ =MAX(1._q,RKLEN*LATT_CUR%BNORM(3)+0.5_q)
            IF (IU6 >= 0 ) THEN
              IF (IU0>=0) &
              WRITE(IU0,99502) NKPX,NKPY,NKPZ
              WRITE(IU6,99502) NKPX,NKPY,NKPZ
            ENDIF
99502       FORMAT( ' generate k-points for:',3I4)
            DO 99501 I=1,3
               BK(I,1)=LATT_CUR%B(I,1)/FLOAT(NKPX)
               BK(I,2)=LATT_CUR%B(I,2)/FLOAT(NKPY)
               BK(I,3)=LATT_CUR%B(I,3)/FLOAT(NKPZ)
99501       CONTINUE
         ELSE
! Here we give the Monkhorst-Pack conventions ... :
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) NKPX,NKPY,NKPZ
! Shift (always in units of the reciprocal lattice vectors!):
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70113,END=70113) SHIFT(1),SHIFT(2),SHIFT(3)
            SUPL_SHIFT=SHIFT
70113       CONTINUE
! Internal rescaling and centering according to Monkhorst-Pack:
            IF ((CSEL=='M').OR.(CSEL=='m')) THEN
               SHIFT(1)=SHIFT(1)+0.5_q*MOD(NKPX+1,2)
               SHIFT(2)=SHIFT(2)+0.5_q*MOD(NKPY+1,2)
               SHIFT(3)=SHIFT(3)+0.5_q*MOD(NKPZ+1,2)
            ENDIF
            DO I=1,3
               BK(I,1)=LATT_CUR%B(I,1)/FLOAT(NKPX)
               BK(I,2)=LATT_CUR%B(I,2)/FLOAT(NKPY)
               BK(I,3)=LATT_CUR%B(I,3)/FLOAT(NKPZ)
            ENDDO
! At this point hope that routine IBZKPT accepts all ... !
         ENDIF
! Find all irreducible points in the first Brillouin zone ... :
         IF (.NOT.LINVERSION) THEN
            CALL VTUTOR('W','LNOINVERSION',RTUT,1, &
     &           ITUT,1,CDUM,1,LDUM,1,IU6,2)
            CALL VTUTOR('W','LNOINVERSION',RTUT,1, &
     &           ITUT,1,CDUM,1,LDUM,1,IU0,2)
         ENDIF
         CALL IBZKPT(LATT_CUR%B(1,1),BK,SHIFT,KPOINTS%NKPTS, &
              VKPT(1,1),WTKPT(1),KPOINTS%NKDIM, &
              KPOINTS%LTET,KPOINTS%NTET,IDTET(0,1),NTETD,KPOINTS%VOLWGT, &
              KPOINTS%IKPT(1,1,1),IKPTD,LATT_CUR%SCALE,LINVERSION,IU6)

         ! k-point lines
         BK_REC=BK
         CALL KARDIR(3,BK_REC,LATT_CUR%A)
         
         CALL XML_KPOINTS_2(KPOINTS%NKPTS, VKPT, WTKPT,&
              KPOINTS%NTET, IDTET, KPOINTS%VOLWGT, &
              CSEL, RKLEN, NKPX, NKPY, NKPZ, SUPL_SHIFT, SHIFT, BK_REC )

      ENDIF

! set old k-points
      GOTO 70222
70111 CONTINUE
      CALL VTUTOR('E','KPOINTS',RTUT,1, &
     &     ITUT,1,CDUM,1,LDUM,1,IU6,1)
      CALL VTUTOR('E','KPOINTS',RTUT,1, &
     &     ITUT,1,CDUM,1,LDUM,1,IU0,1)
      STOP
70119 CONTINUE
      CALL VTUTOR('E','KPOINTSTET',RTUT,1, &
     &     ITUT,1,CDUM,1,LDUM,1,IU6,1)
      CALL VTUTOR('E','KPOINTSTET',RTUT,1, &
     &     ITUT,1,CDUM,1,LDUM,1,IU0,1)
      STOP
70222 CONTINUE
#ifdef gammareal
      IF (KPOINTS%NKPTS/=1._q .OR. &
     &    VKPT(1,1)**2+VKPT(2,1)**2 &
     &   +VKPT(3,1)**2>1E-10_q) THEN

         CALL VTUTOR('E','GAMMAK',RTUT,1, &
     &        ITUT,1,CDUM,1,LDUM,1,IU6,1)
         CALL VTUTOR('E','GAMMAK',RTUT,1, &
     &        ITUT,1,CDUM,1,LDUM,1,IU0,1)
        STOP
      ENDIF
#endif
      DEALLOCATE(KPOINTS%IKPT)
!
!  reallocate everthing with the minimum number of kpoints set
!
      NK=KPOINTS%NKPTS
      NT=MAX(KPOINTS%NTET,1)

      ALLOCATE(KPOINTS%VKPT(3,NK),KPOINTS%WTKPT(NK),KPOINTS%IDTET(0:4,NT))
      KPOINTS%VKPT = VKPT(1:3,1:NK)
      KPOINTS%WTKPT= WTKPT(1:NK)
      KPOINTS%IDTET= IDTET(0:4,1:NT)
      DEALLOCATE(VKPT,WTKPT,IDTET)
      KPOINTS%NKDIM=NK
      CLOSE(UNIT=14)

      RETURN
      END SUBROUTINE
      END MODULE
