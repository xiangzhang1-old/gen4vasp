#include "symbol.inc"

MODULE finite_differences
  USE prec
  IMPLICIT NONE

CONTAINS

!*********************************************************************
! RCS:  $Id: finite_diff.F,v 1.4 2003/06/27 13:22:18 kresse Exp kresse $
!
! calculate second derivatives using finite differences
! implemented by Orest Dubay hopefully alone :)
!
!*********************************************************************

  SUBROUTINE FINITE_DIFF( LSTOP, STEP, NIONS, NTYP, NITYP,MASSES, POS, FORCE, NDISPL, &
       LSDYN, LSFOR, A, B, IU6, IU0)

    USE lattice

    IMPLICIT NONE

    LOGICAL :: LSTOP              ! on return: true to stop main code
    LOGICAL :: LSDYN              ! selective dynamics (yes/ no)
    INTEGER :: NIONS              ! number ions
    INTEGER :: NTYP               ! number of types of ions
    INTEGER :: NITYP(NTYP)        ! number of species
    REAL(q) :: MASSES(NTYP)       ! masses of species
    REAL(q) :: STEP               ! step size
    REAL(q) :: POS(3,NIONS)       ! positions in terms of direct lattice
    REAL(q) :: FORCE(3,NIONS)     ! forces in cartesian coordinates
    REAL(q) :: A(3,3)             ! lattice vectors
    REAL(q) :: B(3,3)             ! reciprocal lattice vectors
    LOGICAL :: LSFOR(3,NIONS)     ! selective 
    INTEGER :: IU6                ! OUTCAR file
    INTEGER :: IU0                ! stdout
    INTEGER :: IUDYNMAT           ! DYNMAT file
    INTEGER :: NDISPL             ! number of displacement

! local variables
    REAL(q) :: X             
    REAL(q),ALLOCATABLE,SAVE :: INITIAL_POSITIONS(:,:)
    REAL(q),ALLOCATABLE,SAVE :: INITIAL_FORCE(:,:)
    REAL(q),ALLOCATABLE,SAVE :: DISPL_FORCES(:,:,:,:)
    REAL(q),ALLOCATABLE      :: SUM_FORCES(:,:,:)
    REAL(q),ALLOCATABLE      :: SECOND_DERIV(:,:)
    INTEGER,SAVE             :: DOF
    LOGICAL,SAVE             :: INIT=.FALSE.
    INTEGER,SAVE             :: PROCESSED_DOF
    INTEGER,SAVE             :: PROCESSED_DISPL
    INTEGER                  :: I,J,K,M,N
    REAL(q),ALLOCATABLE      :: WORK(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVECTORS(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVALUES(:)
    INTEGER                  :: IERROR

    IF (.NOT.INIT) THEN
      CALL COUNT_DOF(NIONS, LSFOR, LSDYN, DOF)
      ALLOCATE(INITIAL_POSITIONS(3,NIONS))
      ALLOCATE(INITIAL_FORCE(3,NIONS))
      ALLOCATE(DISPL_FORCES(NDISPL,DOF,3,NIONS))
      INITIAL_POSITIONS             = POS
      INITIAL_FORCE                 = FORCE
      PROCESSED_DISPL               = NDISPL
      PROCESSED_DOF                 = 0
      IF (IU0>=0) WRITE (IU0,*) 'Finite differences POTIM=',STEP,' DOF=',DOF
      IF (IU6>=0) THEN
        WRITE (IU6,*) 'Finite differences:'
        WRITE (IU6,*) '  Step               POTIM = ',STEP
        WRITE (IU6,*) '  Degrees of freedom DOF   = ',DOF
      END IF
      INIT                          = .TRUE.
!      PRINT *,"NTYP=",NTYP
!      DO I=1,NTYP
!        PRINT *,"spec. ",I,"nr:",NITYP(I),"mass:",MASSES(I)
!      END DO
    END IF

    IF (PROCESSED_DOF>0) THEN
      IF (IU0>=0) THEN
        WRITE (IU0,*)'Finite differences progress:'
        WRITE (IU0,'(A,I3,A,I3)') &
           '  Degree of freedom: ',PROCESSED_DOF,'/',DOF
        WRITE (IU0,'(A,I3,A,I3)') &
           '  Displacement:      ',PROCESSED_DISPL,"/",NDISPL
        WRITE (IU0,'(A,I3,A,I3)') &
           '  Total:             ',(PROCESSED_DOF-1)*NDISPL+PROCESSED_DISPL,&
           '/',DOF*NDISPL
      END IF
     
      IF (IU6>=0) THEN
        WRITE (IU6,*)'Finite differences progress:'
        WRITE (IU6,'(A,I3,A,I3)') &
           '  Degree of freedom: ',PROCESSED_DOF,'/',DOF
        WRITE (IU6,'(A,I3,A,I3)') &
           '  Displacement:      ',PROCESSED_DISPL,"/",NDISPL
        WRITE (IU6,'(A,I3,A,I3)') &
           '  Total:             ',(PROCESSED_DOF-1)*NDISPL+PROCESSED_DISPL,&
           '/',DOF*NDISPL
      END IF
     
      ! CALL PRINT_FORCE(NIONS,FORCE,IU6)
      DISPL_FORCES(PROCESSED_DISPL,PROCESSED_DOF,:,:)=FORCE
    END IF

    PROCESSED_DISPL   = PROCESSED_DISPL+1
    IF (PROCESSED_DISPL>NDISPL) THEN
      PROCESSED_DISPL = 1
      PROCESSED_DOF   = PROCESSED_DOF+1
    END IF

    IF (PROCESSED_DOF.LE.DOF) THEN
      LSTOP=.FALSE.
      ! IF (IU0>=0) WRITE (IU0,*) 'Finite differences DOF             = ',DOF
      ! IF (IU0>=0) WRITE (IU0,*) 'Finite differences PROCESSED_DOF   = ',PROCESSED_DOF
      ! IF (IU0>=0) WRITE (IU0,*) 'Finite differences PROCESSED_DISPL = ',PROCESSED_DISPL
      ! IF (IU0>=0) WRITE (IU0,*) 'Finite differences POSITIONS BEFORE DISPL'

      POS=INITIAL_POSITIONS
      ! CALL PRINT_POSITIONS(NIONS,POS,IU0,A,B)

    ! PROCESSED_DISPL ordering:
    !
    ! \_           _/   
    !   \__     __/
    !      \---/
    !  |  |  |  |  |
    !       (0) 1     NDISPL=1
    !     2     1     NDISPL=2
    !  4  3     2  1  NDISPL=4
    
      IF(NDISPL<=2) THEN
         SELECT CASE(PROCESSED_DISPL)
         CASE(1)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF,   STEP,A,B)
         CASE(2)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF,  -STEP,A,B)
         END SELECT
      ELSE
         SELECT CASE(PROCESSED_DISPL)
         CASE(1)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF, 2*STEP,A,B)
         CASE(2)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF,   STEP,A,B)
         CASE(3)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF,  -STEP,A,B)
         CASE(4)
            CALL MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,PROCESSED_DOF,-2*STEP,A,B)
         END SELECT
      END IF
      ! IF (IU0>=0) WRITE (IU0,*) 'Finite differences POSITIONS AFTER DISPL'
      ! CALL PRINT_POSITIONS(NIONS,POS,IU0,A,B)

      RETURN
    ELSE

    !
    ! Final processing + output is here:
    !

      IF (IU6>=0) THEN
        WRITE (IU6,*)
        WRITE (IU6,*) 'FORCES'
        WRITE (IU6,*) '------'
        WRITE (IU6,*)
        WRITE (IU6,*) 'INITIAL FORCE'
        CALL PRINT_FORCE(NIONS,INITIAL_FORCE,IU6)

        DO N=1,DOF
          CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
          DO M=1,NDISPL
            WRITE (IU6,'("DOF:",I4," ATOM:",I4," AXIS:",I4," DISPLACEMENT:",I4)')N,J,I,M
            CALL PRINT_FORCE(NIONS,DISPL_FORCES(M,N,:,:),IU6)
            WRITE (IU6,*)
          END DO
        END DO
        WRITE (IU6,*) ' '
      END IF

      ALLOCATE(SUM_FORCES(DOF,3,NIONS))
      SELECT CASE(NDISPL)
      CASE(1)
        DO N=1,DOF
          SUM_FORCES(N,:,:)=(DISPL_FORCES(1,N,:,:)-INITIAL_FORCE)/STEP
        END DO
      CASE(2)
        SUM_FORCES = (1._q/(2._q*STEP))*(DISPL_FORCES(1,:,:,:)-DISPL_FORCES(2,:,:,:))
      CASE(4)
        SUM_FORCES = (1._q/(12._q*STEP))* &
                                (8._q*DISPL_FORCES(2,:,:,:)-8._q*DISPL_FORCES(3,:,:,:) &
                                     -DISPL_FORCES(1,:,:,:)+     DISPL_FORCES(4,:,:,:))
      END SELECT

      IF (IU6>=0) THEN 
        WRITE (IU6,*) 'DYNMAT'
        WRITE (IU6,*) '------'
        CALL PRINT_DYNMAT(NIONS,DOF,STEP,NTYP,NITYP,MASSES,STEP*SUM_FORCES,LSDYN,LSFOR,IU6)

       ! DO N=1,DOF
       !   CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
       !   WRITE (IU6,'(I4,I4,I4,I4)') N,I,J
       !   CALL PRINT_FORCE(NIONS,SUM_FORCES(N,:,:),IU6)
       ! END DO
        WRITE (IU6,*) ' '
      END IF      
     
      ALLOCATE(SECOND_DERIV(DOF,DOF))
      DO N=1,DOF
        CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)       
        DO M=1,DOF
          SECOND_DERIV(M,N)=SUM_FORCES(M,I,J)
        END DO
      END DO
     
      IF (IU6>=0) THEN
        WRITE (IU6,*) 
        WRITE (IU6,*) 'SECOND DERIVATIVES (NOT SYMMETRIZED)'
        WRITE (IU6,*) '------------------------------------'
        CALL PRINT_SECOND_DERIV(NIONS,DOF,SECOND_DERIV,LSFOR,LSDYN,IU6)
      END IF

      DO N=1,DOF
        DO M=N+1,DOF
          X=0.5_q*(SECOND_DERIV(N,M)+SECOND_DERIV(M,N))
          !WRITE(0,*) N,M,SECOND_DERIV(N,M),SECOND_DERIV(M,N),X
          SECOND_DERIV(N,M)=X
          SECOND_DERIV(M,N)=X
          !WRITE(0,*) N,M,SECOND_DERIV(N,M),SECOND_DERIV(M,N),X
        END DO
      END DO

      IF (IU6>=0) THEN
        WRITE (IU6,*)  
        WRITE (IU6,*) 'SECOND DERIVATIVES (SYMMETRYZED)'
        WRITE (IU6,*) '--------------------------------'
        CALL PRINT_SECOND_DERIV(NIONS,DOF,SECOND_DERIV,LSFOR,LSDYN,IU6)
      END IF

      N=1
      DO I=1,NTYP
        DO J=1,NITYP(I)
          DO K=1,3
            CALL FIND_DOF_INDEX(NIONS,LSFOR,LSDYN,K,N,M)
            IF (M>0) SECOND_DERIV(:,M)=SECOND_DERIV(:,M)/SQRT(MASSES(I))
            IF (M>0) SECOND_DERIV(M,:)=SECOND_DERIV(M,:)/SQRT(MASSES(I))
          END DO
          N=N+1
        END DO
      END DO

      IF (IU6>=0) THEN
        WRITE (IU6,*)  
        WRITE (IU6,*) 'MASS-WEIGHTED SECOND DERIVATIVES'
        WRITE (IU6,*) '--------------------------------'
        CALL PRINT_SECOND_DERIV(NIONS,DOF,SECOND_DERIV,LSFOR,LSDYN,IU6)
      END IF
      
      ALLOCATE(WORK(DOF,32),EIGENVECTORS(DOF,DOF),EIGENVALUES(DOF))
      EIGENVECTORS=SECOND_DERIV

      CALL DSYEV &
              ('V','U',DOF,EIGENVECTORS,DOF, &
              EIGENVALUES,WORK,32*DOF, IERROR)

      IF (IERROR/=0) THEN
        IF (IU6>=0) THEN
          WRITE(IU6,*) "Error while diagonalisation DSYEV INFO=",IERROR
          WRITE(IU6,*) "Some of (or all) eigenvectors and eigenvalues are not correct !"
        END IF
      END IF

      CALL PRINT_EIGENVECTORS(NIONS,DOF,INITIAL_POSITIONS,A, &
                                        EIGENVECTORS,      &
                                        EIGENVALUES,       &
                                        LSFOR,LSDYN,IU6)
      IF (IU6>=0) THEN
        WRITE(IU6,*) "Eigenvectors after division by SQRT(mass)"
      END IF
     
      N=1
      DO I=1,NTYP
        DO J=1,NITYP(I)
          DO K=1,3
            CALL FIND_DOF_INDEX(NIONS,LSFOR,LSDYN,K,N,M)
            IF (M>0) EIGENVECTORS(M,:)=EIGENVECTORS(M,:)/SQRT(MASSES(I))
          END DO
          N=N+1
        END DO
      END DO

      CALL PRINT_EIGENVECTORS(NIONS,DOF,INITIAL_POSITIONS,A, &
                                        EIGENVECTORS,      &
                                        EIGENVALUES,       &
                                        LSFOR,LSDYN,IU6)


      DEALLOCATE(WORK)
      DEALLOCATE(INITIAL_POSITIONS)
      DEALLOCATE(INITIAL_FORCE)
      DEALLOCATE(DISPL_FORCES)
      DEALLOCATE(SECOND_DERIV,SUM_FORCES)
      LSTOP=.TRUE.
    END IF


    IF (IU0>=0) WRITE (IU0,*) 'Finite differences POTIM=',STEP
    IF (IU6>=0) WRITE (IU6,*) 'Finite differences POTIM=',STEP

    ! converts a vector with three entries from direct to cartesian
    !CALL DIRKAR(1,X,A)
    
    ! converts a vector from cartesian to direct coordinates
    !CALL KARDIR(1,X,B)


!         CALL DSYEV &
!              ('V','U',NDIM,MATRIX,NDIM, &
!              EIGENVALUES,CWORK,*NDIM, IERROR)

    LSTOP=.TRUE.
  END SUBROUTINE FINITE_DIFF

  SUBROUTINE COUNT_DOF(NIONS, LSFOR, LSDYN, N)
    INTEGER :: NIONS,N
    LOGICAL :: LSFOR(3,NIONS)  ! selective 
    LOGICAL :: LSDYN
    ! local    
    INTEGER :: I,J

    IF (.NOT.LSDYN) THEN
      N=NIONS*3
      RETURN
    END IF

    N=0
    
    DO J=1,NIONS
      DO I=1,3
        IF (LSFOR(I,J)) N=N+1
      END DO
    END DO
  END SUBROUTINE COUNT_DOF

  ! Find the indexes of N-th degree of freedom
  SUBROUTINE FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
    INTEGER :: NIONS,N
    LOGICAL :: LSFOR(3,NIONS)
    LOGICAL :: LSDYN
    INTEGER :: I,J
    !local
    INTEGER :: M

    IF (.NOT.LSDYN) THEN
      I=1+MOD(N-1,3)
      J=1+(N-1)/3
      RETURN
    END IF

    M=0
ijloop: DO J=1,NIONS
      DO I=1,3
        IF (LSFOR(I,J)) M=M+1
        IF (M==N) EXIT ijloop
      END DO
    END DO ijloop
  END SUBROUTINE FIND_IJ

  SUBROUTINE FIND_DOF_INDEX(NIONS,LSFOR,LSDYN,AXE,N,DOFI)
    INTEGER     :: NIONS,N,AXE
    INTEGER     :: DOFI
    LOGICAL     :: LSFOR(3,NIONS)
    LOGICAL     :: LSDYN
    ! local    
    INTEGER :: I,J

    IF (.NOT.LSDYN) THEN
      DOFI=3*N+AXE-3
      RETURN
    END IF

    DOFI=0
    
    DO J=1,NIONS
      DO I=1,3
        IF (LSFOR(I,J)) THEN
          DOFI=DOFI+1
          IF ( (J==N).AND.(I==AXE) ) RETURN
        END IF
      END DO
    END DO
    DOFI=-1
    RETURN
  END SUBROUTINE FIND_DOF_INDEX

  SUBROUTINE MAKE_DISPLACED(NIONS,LSFOR,LSDYN,POS,N,STEP,A,B)
    USE lattice
    INTEGER :: NIONS,N
    REAL(q) :: POS(3,NIONS),STEP
    LOGICAL :: LSFOR(3,NIONS)
    LOGICAL :: LSDYN
    REAL(q) :: A(3,3)          ! lattice vectors
    REAL(q) :: B(3,3)          ! reciprocal lattice vectors

    !local
    INTEGER::I,J

    CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
    CALL DIRKAR(1,POS(:,J),A)
    POS(I,J) = POS(I,J)+STEP
    CALL KARDIR(1,POS(:,J),B)
  END SUBROUTINE MAKE_DISPLACED

  SUBROUTINE PRINT_POSITIONS(NIONS,POS,OUT,A,B)
    USE lattice
    INTEGER :: NIONS,OUT
    REAL(q) :: POS(3,NIONS)
    REAL(q) :: A(3,3),B(3,3)
    !local
    INTEGER I
    REAL(q) :: X(3)
    IF (OUT>=0) THEN
      WRITE (OUT,*) 'Positions of atoms (Carthesian coordinates)'
      WRITE (OUT,*) '-------------------------------------------'
      DO I=1,NIONS
        X=POS(:,I)
        CALL DIRKAR(1,X,A)
        WRITE(OUT,'(3F11.7)') X
      END DO
    END IF
  END SUBROUTINE PRINT_POSITIONS

  SUBROUTINE PRINT_FORCE(NIONS,FORCE,OUT)
    USE lattice
    INTEGER :: NIONS,OUT
    REAL(q) :: FORCE(3,NIONS)
    !local
    INTEGER I
    IF (OUT>=0) THEN
      DO I=1,NIONS
        WRITE(OUT,'(F10.6," ",F10.6," ",F10.6)') FORCE(:,I)
      END DO
    END IF
  END SUBROUTINE PRINT_FORCE

  SUBROUTINE PRINT_DYNMAT(NIONS,DOF,STEP,NTYP,NITYP,MASSES,FORCES,LSDYN,LSFOR,OUT)
    USE lattice
    INTEGER :: NIONS,OUT,DOF,NTYP
    INTEGER :: NITYP(NTYP)
    REAL(q) :: MASSES(NTYP)
    REAL(q) :: STEP
    REAL(q) :: FORCES(DOF,3,NIONS)
    LOGICAL :: LSFOR(3,NIONS)
    LOGICAL :: LSDYN    
    !local
    INTEGER :: N,I,J,OLDJ,DISPL

    IF (OUT>=0) THEN
       WRITE(OUT,'(I5)',ADVANCE='NO') NTYP
       DO N=1,NTYP
          WRITE(OUT,'(I5)',ADVANCE='NO') NITYP(N) 
       END DO
       WRITE(OUT,'(I5)') DOF
       DO N=1,NTYP
          WRITE(OUT,'(F7.3)',ADVANCE='NO') MASSES(N) 
       END DO
       WRITE(OUT,*)

       OLDJ  = 0
       DISPL = 0
       DO N=1,DOF
          CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
          IF (J==OLDJ) THEN
            DISPL=DISPL+1
          ELSE
            DISPL=1
            OLDJ=J
          ENDIF
          WRITE (OUT,'(I5,I5)',ADVANCE='NO') J,DISPL
          SELECT CASE(I)
          CASE(1)
             WRITE (OUT,'(F8.4,F8.4,F8.4)') STEP,0._q,0._q
          CASE(2)
             WRITE (OUT,'(F8.4,F8.4,F8.4)') 0._q,STEP,0._q
          CASE(3)
             WRITE (OUT,'(F8.4,F8.4,F8.4)') 0._q,0._q,STEP
          CASE DEFAULT
             WRITE (OUT,*) "?"
          END SELECT                            
          CALL PRINT_FORCE(NIONS,FORCES(N,:,:),OUT)
       END DO
    END IF
  END SUBROUTINE PRINT_DYNMAT

  SUBROUTINE PRINT_SECOND_DERIV(NIONS,DOF,SD,LSFOR,LSDYN,OUT)
    INTEGER::NIONS,DOF,OUT
    LOGICAL::LSFOR(3,NIONS)
    LOGICAL::LSDYN
    REAL(q)::SD(DOF,DOF)
    !local
    INTEGER   :: I,J,K,M,N
    CHARACTER :: C

    IF (OUT>=0) THEN
      WRITE(OUT,'(A)',ADVANCE='NO') "      "
      DO N=1,DOF
        CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
        SELECT CASE(I)
        CASE(1)
           C="X"
        CASE(2)
           C="Y"
        CASE(3)
           C="Z"
        CASE DEFAULT
           C="?"
        END SELECT
        WRITE(OUT,'(I10,A,A)',ADVANCE='NO') J,C," "
      END DO
      WRITE(OUT,*)

      DO N=1,DOF
        CALL FIND_IJ(NIONS,LSFOR,LSDYN,N,I,J)
        SELECT CASE(I)
        CASE(1)
           C="X"
        CASE(2)
           C="Y"
        CASE(3)
           C="Z"
        CASE DEFAULT
           C="?"
        END SELECT
        WRITE(OUT,'(I3,A,A)',ADVANCE='NO') J,C," "
     
        DO M=1,DOF
          WRITE(OUT,'(F12.6)',ADVANCE='NO') SD(N,M)
        END DO
        WRITE(OUT,*)
      END DO
      WRITE(OUT,*)
    END IF
  END SUBROUTINE PRINT_SECOND_DERIV

  SUBROUTINE PRINT_EIGENVECTORS (NIONS,DOF,POS,A,EVEC,EVAL,LSFOR,LSDYN,OUT)
    USE constant
    USE lattice
    INTEGER::NIONS,DOF,OUT
    LOGICAL::LSFOR(3,NIONS)
    LOGICAL::LSDYN
    REAL(q)::EVEC(DOF,DOF),EVAL(DOF)
    REAL(q)::POS(3,NIONS)
    REAL(q)::A(3,3)
    !local
    INTEGER   :: N,NI,AXE,DI
    REAL(q),PARAMETER :: PLANK=6.626075E-34
    REAL(q),PARAMETER :: C= 2.99792458E10
    REAL(q)::X(3)
!    REAL(q)   :: ELECT  = 1.602199E-19
!    REAL(q)   :: M0     = 1.6725E-27
    REAL(q)   :: FACTOR,W
                                               
!---- frequenz Sqrt(d E / M) d x (cgs System)
    FACTOR=SQRT(EVTOJ/((1E-10)**2)/AMTOKG)

    IF (OUT>=0) THEN
      WRITE(OUT,*)
      WRITE(OUT,*)'Eigenvectors and eigenvalues of the dynamical matrix'
      WRITE(OUT,*)'----------------------------------------------------'
      WRITE(OUT,*)
      DO N=1,DOF
        WRITE(OUT,*)
        W=FACTOR*SQRT(ABS(EVAL(N)))
        IF (EVAL(N).GT.0) THEN
          WRITE(OUT,'(I4," f/i=",F12.6," THz ",F12.6," 2PiTHz",F12.6," cm-1 ",F12.6," meV")') &
               N,W/(1E12*2*PI),W/(1E12),W/(C*PI*2),W*1000*PLANK/EVTOJ/2/PI
        ELSE
          WRITE(OUT,'(I4," f  =",F12.6," THz ",F12.6," 2PiTHz",F12.6," cm-1 ",F12.6," meV")') &
               N,W/(1E12*2*PI),W/(1E12),W/(C*PI*2),W*1000*PLANK/EVTOJ/2/PI
        END IF
        WRITE(OUT,'("             X         Y         Z","           dx          dy          dz")')
        DO NI=1,NIONS
          X=POS(:,NI)
          CALL  DIRKAR(1,X,A)
          WRITE (OUT,'(A,3F10.6,A)',ADVANCE='NO')'    ',X,'   '
          DO AXE=1,3
            CALL FIND_DOF_INDEX(NIONS,LSFOR,LSDYN,AXE,NI,DI)
            IF (DI>0) THEN
              WRITE(OUT,'(F10.6,A)',ADVANCE='NO') EVEC(DI,N),'  '
            ELSE
              WRITE(OUT,'(A)',ADVANCE='NO')'         0  '
            END IF
          END DO
          WRITE(OUT,*)
        END DO
      END DO
      WRITE(OUT,*)
    END IF
  END SUBROUTINE PRINT_EIGENVECTORS
  
END MODULE finite_differences
