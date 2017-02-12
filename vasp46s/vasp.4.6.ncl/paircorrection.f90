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





!*********************************************************************
!   SUBROUTINE PAIR_CORRECTION
!
! this subroutine can be used to calculate the effect of an additional
! pairwise additative potential
! the resulting energy is returned in ENERGY and the forces are
! added to the array FORCES
! the original routine was written by Gilles de Wijs
! it was cleaned up a little bit by gK
!
!*********************************************************************

      SUBROUTINE PAIR_CORRECTION(NIONS, NTYP, NI_TYP, &
           A,BNORM,XC,FORCES,ENERGY,IU5,IU6)

      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)

      PARAMETER (NDIM=100)        ! maximal number of points used for the presentation 
                                  ! of the pair potential
      PARAMETER (VOLF=1.005)      ! finite differences are used to evaluate the stress tensor
                                  ! set the magnitude of the volume change 
      INTEGER NIONS               ! number of ions
      INTEGER NTYP                ! number of species
      INTEGER NI_TYP(NTYP)        ! number of atoms per species
 
      REAL(q) XC(3,NIONS)         ! position of each ion

      REAL(q) A(3,3),BNORM(3)
      REAL(q) FORCES(3,NIONS)
      REAL(q),ALLOCATABLE :: CHECKF(:,:)
      REAL(q) V(3),TMP(NDIM)
      ! local saved variables
      REAL(q),SAVE :: PSP(NDIM,5),APACO,XMIN,XMAX,ARGSC
      INTEGER,SAVE :: NPOINTS
      ! status of routine (0=uninitialized, 1=used/initialized 2=unused)
      INTEGER :: STATUS=0
      CHARACTER*1 CHARAC
      LOGICAL LDUM


      ENERGY=0
!
! status== 2 immedate return
!
      IF ( STATUS == 2) RETURN
!
! status== 0 initialize (try reading PAIRPOT from INCAR file)
!
      IF ( STATUS == 0 ) THEN 
         CALL RDATAB(.TRUE.,'INCAR',IU5,'PAIRPOT','=','#',';','F', &
     &     IDUM,TMP,CDUM,LDUM,CHARAC,N,NDIM,IERR)
         IF ((IERR/=0) .AND. (IERR/=3)) THEN
            WRITE(*,*)'Error reading item ''PAIRPOT'' from INCAR.'
            STOP
         ENDIF
         IF (IERR/=0 .OR. N==0) THEN
            STATUS=2
            RETURN
         ENDIF
         STATUS=1

         NPOINTS=N/2
         IF (IU6>0) THEN
            WRITE(IU6,"('found', I4,' pair potential data on INCAR')") NPOINTS
         ENDIF

         PSP(1:NPOINTS,1)=TMP(1:NPOINTS)
         PSP(1:NPOINTS,2)=TMP(NPOINTS+1:NPOINTS*2)
         
         APACO =PSP(NPOINTS,1)
         XMIN  =PSP(1,1)

         DO I=1,NPOINTS
            PSP(I,2)=PSP(I,2)-PSP(NPOINTS,2)
            PSP(I,1)=PSP(I,1)-XMIN
         ENDDO

         XMAX=PSP(NPOINTS,1)

         CALL SPLCOF(PSP(1,1),NPOINTS,NDIM,1E30_q)

         ARGSC=(NPOINTS-1)/XMAX

      ENDIF
!-----------------------------------------------------------------------
!  pressure is currently only correct for cubic cells
!-----------------------------------------------------------------------
      VOL= SQRT(A(1,1)**2+A(2,1)**2+A(3,1)**2)* &
     &     SQRT(A(1,2)**2+A(2,2)**2+A(3,2)**2)* &
     &     SQRT(A(1,3)**2+A(2,3)**2+A(3,3)**2)
      VOLP=VOLF
      VOLM=2.-VOLF
      VOLD=VOL*(VOLP-VOLM)
      VOLP=EXP(LOG(VOLP)/3.)
      VOLM=EXP(LOG(VOLM)/3.)
      AA1=ABS(A(1,1)*A(1,2)+A(2,1)*A(2,2)+A(3,1)*A(3,2))
      AA2=ABS(A(1,3)*A(1,2)+A(2,3)*A(2,2)+A(3,3)*A(3,2))
      AA3=ABS(A(1,1)*A(1,3)+A(2,1)*A(2,3)+A(3,1)*A(3,3))
      IF ( &
     &    (AA1.GT.0.001) .OR. &
     &    (AA2.GT.0.001) .OR. &
     &    (AA3.GT.0.001)) THEN
        IF (IU6.GT.0) WRITE(IU6,*) 'CELL IS NOT SC, PRESSURE CALCULATION INCORRECT'
      ENDIF

      ENERGY =0.
      ENERGYP=0.
      ENERGYM=0.
      
      ALLOCATE(CHECKF(3,NIONS))
      CHECKF=0
!-----------------------------------------------------------------------
!  IXMAX,Y,Z defines the maximum number of cells which
!  have to be used for the summation
!-----------------------------------------------------------------------
      I1MAX=AINT(APACO*BNORM(1)-.001)
      I2MAX=AINT(APACO*BNORM(2)-.001)
      I3MAX=AINT(APACO*BNORM(3)-.001)

      APACO2=APACO**2
!-----------------------------------------------------------------------
      DO I1=-I1MAX-1,I1MAX
      DO I2=-I2MAX-1,I2MAX
      DO I3=-I3MAX-1,I3MAX

      DO NI=1,NIONS
      DO NII=1,NIONS
         IF (NI==NII) CYCLE
         D1= I1+MOD(XC(1,NI)-XC(1,NII)+1,1._q)
         D2= I2+MOD(XC(2,NI)-XC(2,NII)+1,1._q)
         D3= I3+MOD(XC(3,NI)-XC(3,NII)+1,1._q)

         V(1)= D1*A(1,1)+D2*A(1,2)+D3*A(1,3)
         V(2)= D1*A(2,1)+D2*A(2,2)+D3*A(2,3)
         V(3)= D1*A(3,1)+D2*A(3,2)+D3*A(3,3)
         R2= V(1)**2 + V(2)**2 + V(3)**2

         IF (R2 > APACO2) CYCLE

         R2=SQRT(R2)
         R2P=R2*VOLP
         R2M=R2*VOLM

         DIST=R2
         DO J=1,3
            V(J)=V(J)/R2
         ENDDO

         CALL INTERP(R2 -XMIN,E ,DE ,PSP(1,1),NPOINTS,NDIM,ARGSC,IU6)
         CALL INTERP(R2P-XMIN,EP,DEP,PSP(1,1),NPOINTS,NDIM,ARGSC,IU6)
         CALL INTERP(R2M-XMIN,EM,DEM,PSP(1,1),NPOINTS,NDIM,ARGSC,IU6)
         ! double counting
         E =E /2.
         EP=EP/2.
         EM=EM/2.
!       WRITE (*,*) 'OO',R2,APACO2,ENERGY,I
!       WRITE (*,*) 'NI',NI,NII,' DE',DE,' DIST',DIST,' V1',V

         ENERGY =ENERGY +E
         ENERGYP=ENERGYP+EP
         ENERGYM=ENERGYM+EM
       
         DO J=1,3
            FORCES(J,NI) =FORCES(J,NI) -V(J)*DE
            CHECKF(J,NI) =CHECKF(J,NI) -V(J)*DE
         ENDDO

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      PRESS=(ENERGYM-ENERGYP)/VOLD
!     EV/A**3 = 1.6022E-19 M^2 KG / S^2 / A**3 =
!               1.6022E-19 M^3 PA / A**3 =
!               1.6022E-19 M^3 BAR / A**3 / 100 000 =
!               1.6022E-19 BAR 10^30 / 100 000 = 1.6022E6 BAR

      IF (IU6>=0) THEN
         WRITE(IU6,*) 'PAIR: correction to pressure:',PRESS * 1.6022E3,'KB'
         WRITE(IU6,*) 'energy correction',ENERGY
         WRITE(IU6,*) 'force corrections'
         DO J=1,NIONS
            WRITE(IU6,*) CHECKF(1,J),CHECKF(2,J),CHECKF(3,J)
         ENDDO
      ENDIF

      DEALLOCATE (CHECKF)

      RETURN
      END

!=======================================================================
!
!  interpolate the array PSP using evaluated spline coefficients
!
!=======================================================================

	SUBROUTINE INTERP(RIN,E,DE,PSP,NPOINTS,NDIM,ARGSC,IU6)
        USE prec
        IMPLICIT REAL(q) (A-H,O-Z)
        DIMENSION PSP(NDIM,5)

        R2=RIN
        IF (R2 >= PSP(NPOINTS,1)) THEN
          E=0.
          DE=0.
          RETURN
        ENDIF
        
        IF (R2.LE.0) THEN
          IF (IU6.GT.0) WRITE(IU6,*) 'PAIR, WARNING: R2 TOO SMALL'
          R2=0.
        ENDIF
        I  =MAX(INT(R2*ARGSC)+1,1)
        REM=R2-PSP(I,1)
        E  =(PSP(I,2)+REM*(PSP(I,3)+ &
     &                REM*(PSP(I,4)  +REM*PSP(I,5))))
        DE = PSP(I,3)+REM*(PSP(I,4)*2+REM*PSP(I,5)*3)
        RETURN
        END
