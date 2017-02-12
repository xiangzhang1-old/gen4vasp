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





      MODULE ini
      USE prec
      END MODULE
!**************** SUBROUTINE SPLCOF, SPLCOF_N0 *************************
! RCS:  $Id: ini.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!  Subroutine for calculating Spline-Koeffizients
!  using the routines of the book 'numerical  recipes'
!  on input P(1,N) must contain x-values
!           P(2,N) must contain function-values
!  YP is the first derivatives at the first point
!  if >= 10^30 natural boundary-contitions (y''=0) are used
!
!  for point N always natural boundary-conditions are used in
!  SPLCOF, whereas SPLCOF_N0 assume 0 derivative at N
!
!***********************************************************************

      SUBROUTINE SPLCOF(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO

      P(N,4)=0.0_q
      P(N,3)=0.0_q
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE



      SUBROUTINE SPLCOF_N0(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO
      YNP=0
      IF (YNP> .99E30_q) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5_q
        UN=(3._q/(P(N,1)-P(N-1,1)))*(YNP-(P(N,2)-P(N-1,2))/ &
     &             (P(N,1)-P(N-1,1)))
      ENDIF
      P(N,4)=(UN-QN*P(N-1,3))/(QN*P(N-1,4)+1.)
      P(N,3)=0  ! never used
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE

!
!  helper routine, which copies X and Y arrays to P
!  and than performes the fit on the array Y
!
      SUBROUTINE SPLCPY(X,Y,P,NAC,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
      DIMENSION X(NAC)
      DIMENSION Y(NAC)
      DO 100 N=1,NAC
        P(N,1)=X(N)
        P(N,2)=Y(N)
  100 CONTINUE
      CALL SPLCOF(P,NAC,NDIM,Y1P)
      RETURN
      END SUBROUTINE
!
!  helper routine, which evaluates the spline fit at a specific
!  position
!
      SUBROUTINE SPLVAL(X,F,FDER,P,NAC,NDIM)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!  interval bisectioning
      I=1
      J=NAC
      IF (X   <P(I,1)) GO TO 60
      IF (X   <P(J,1)) GO TO 70
      K=J-1
      GOTO 90
   60 K=1
      GOTO 90
   70 K=(I+J)/2
      IF(I==K) GOTO 90
      IF (X   <P(K,1)) GO TO 80
      I=K
      GOTO 70
   80 J=K
      GOTO 70
!
   90 DX=X   -P(K,1)
      F   =((P(K,5)*DX+P(K,4))*DX+P(K,3))*DX+P(K,2)
      FDER=(3.0_q*P(K,5)*DX+2.0_q*P(K,4))*DX+P(K,3)
  100 RETURN
      END SUBROUTINE


!***********************************************************************
!
! system name date and time
!
!***********************************************************************

      SUBROUTINE MY_DATE_AND_TIME(IU6)
      USE prec
      USE vaspxml
      IMPLICIT NONE
      INTEGER IU6
      CHARACTER*(8)  DATE
      CHARACTER*(10) TIME

      CALL DATE_AND_TIME( DATE,TIME)
      IF (IU6>=0) &
      WRITE(IU6,"(' executed on ',A20,' date ',A4,'.',A2,'.',A2,'  ',A2,':',A2,':',A2 )") & 
           "LinuxIFC",DATE(1:4),DATE(5:6),DATE(7:8),TIME(1:2),TIME(3:4),TIME(5:6)

      CALL XML_TAG_STRING("platform" , "LinuxIFC")
      CALL XML_TAG_STRING("date" , DATE(1:4)//" "//DATE(5:6)//" "//DATE(7:8))
      CALL XML_TAG_STRING("time" , TIME(1:2)//":"//TIME(3:4)//":"//TIME(5:6))

      END SUBROUTINE

!***********************************************************************
!
! write out current memory requirements
!
!***********************************************************************

      SUBROUTINE MEMORY_CHECK(LOOP,STR)
      USE prec
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      CHARACTER*(*) STR
      REAL*8 SUM

      CALL TIMING(0,UTIME,STIME,DAYTIM,MINPGF,MAJPGF, &
     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,'(A)') &
     &   ' General timing and accounting informations for this job:'
      WRITE(*,'(A)') &
     &   ' ========================================================'
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,F12.0)') '  Maximum memory used (kb): ',RSIZM
      WRITE(*,'(17X,A,F12.0)') '  Average memory used (kb): ',AVSIZ
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,I12)')   '         Minor page faults: ',MINPGF
      WRITE(*,'(17X,A,I12)')   '         Major page faults: ',MAJPGF
      WRITE(*,'(17X,A,I12)')   'Voluntary context switches: ',IVCSW

      END SUBROUTINE
